#!/usr/bin/perl
#use strict;
# read in HMM results and determines the start and end of enriched regions.
# changed by JC at 08/11/2008

$minregionlength = 2; #3 * 25 = 75 bp, 2 * 25 = 50 bp.

$argc = @ARGV;
$argc == 7 || die "Provide probability threshold, species, significance level, window size, input filename, parameter file and result file name.\n";

$sp2 = $ARGV[1];
if ($sp2 eq "human" || "hg19") 
{
        $nchr2 = 24;
        $fold2 = "hg19";
	$gcover = 0.9*3100000000;
}
elsif ($sp2 eq "hg18") 
{
        $nchr2 = 24;
        $fold2 = "hg18";
	$gcover = 0.9*3100000000;
}
elsif ($sp2 eq "mouse" || "mm9") 
{
        $nchr2 = 21;
        $fold2 = "mm9";
	$gcover = 0.7*2655000000;
}       
else
{
        print "species error!\n";
        exit(0); 
}

$threshold = $ARGV[0];
$siglevel = $ARGV[2];
$windowSize = $ARGV[3];
$inputfilename = $ARGV[4];
$parasfilename = $ARGV[5];
$outputfilename = $ARGV[6];
my @chr;
my @startpos;
my @endpos;
my @regionwid;
my @supheight;
my @totalmaxproba;
$totalcount = 0;
for($j=1;$j<($nchr2 + 1);$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	my $reached = 0;
	while($datline = <INPUTFILE>)
	{
		chomp($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
		$height = $datas[3];
                $logproba = $datas[4];
		if($proba >= $threshold) ##### changed 08/10/08
		{
			$reached = 1;
 			if($first == 1)
			{
				$start = $order;
				$end = $order;
				$maxheight = $height;
				$first = 0;
				$length = 1;
				$count = 1;
				$sumlogproba = $logproba;
                                if($proba >=1)       
                                {
                                        $maxprobacount =1;
                                }
                                else
                                {
                                        $maxprobacount =0;
                                }
			}
			else
			{
				if($order == ($end +1))
				{# continue
					$end ++;
					$length ++;
					if($height >$maxheight)
					{
						$maxheight = $height;
					}
					$sumlogproba = $sumlogproba + $logproba;
                                        if($proba >=1)
                                        {
                                                $maxprobacount ++;
                                        }
				}
				else
				{# start a new region
					$chr[$totalcount] = $j;
					$startpos[$totalcount] = ($start-1)*$windowSize+1;
					$endpos[$totalcount] = $end*$windowSize;
					$regionwid[$totalcount] = $length *$windowSize;
					$supheight[$totalcount] = $maxheight;	
					$totallogproba[$totalcount] = $sumlogproba;				
					$totalmaxproba[$totalcount] = $maxprobacount;
					$start = $order;
					$end = $order;
					$maxheight = $height;
        	                        $length = 1;
	                                $sumlogproba = $logproba;
                                        if($proba >=1)
                                        {
                                                $maxprobacount =1;
                                        }
                                        else
                                        {
                                                $maxprobacount =0;
                                        }                	                
					$count ++;
					$totalcount ++;
				}
			}
		}
	}
	if($reached ==1)
	{
		$chr[$totalcount] = $j;
        	$startpos[$totalcount] = ($start-1)*$windowSize+1;
        	$endpos[$totalcount] = $end*$windowSize;          
        	$regionwid[$totalcount] = $length *$windowSize;
        	$supheight[$totalcount] = $maxheight;
		$totallogproba[$totalcount] = $sumlogproba;
                $totalmaxproba[$totalcount] = $maxprobacount;
		$totalcount ++;
		print "totalcount = $totalcount\n";
	}
}
close(INPUTFILE);

open(PARASFILE, $parasfilename);
for($j=0;$j<4;$j++)
{
	$datline = <PARASFILE>;           
}
chomp($datline);
@datas = split(" ", $datline); 
$fgread = $datas[0];
$datline = <PARASFILE>;  
$datline = <PARASFILE>;
chomp($datline);             
@datas = split(" ", $datline);
$fragmentWidth = $datas[0];

print "foreground read = $fgread\n";
close(PARASFILE);
$lambda = $fgread*$fragmentWidth/$gcover;
$probcutoff = $siglevel/$totalcount;
$heightcutoff = &decide($lambda,$probcutoff);
print "There are $totalcount regions.\n";
print "lambda= $lambda.\n";
print "max height = $heightcutoff\n";
open(OUTPUTFILE, ">".$outputfilename.".allregions.txt");
for($j=0;$j<$totalcount;$j++)
{
        if(($totalmaxproba[$j] >= $minregionlength)&&($regionwid[$j] > $minregionlength *$windowSize))
        {        
               print OUTPUTFILE "$chr[$j]\t$startpos[$j]\t$endpos[$j]\t$regionwid[$j]\t$supheight[$j]\t$totallogproba[$j]\n";
        }			
}
close(OUTPUTFILE);

sub decide()
{
	local($mylambda,$myprobcutoff)= @_;
	my @proba;
	my $j;
	my $k;
	my $m;
	my $middle = 30;
        my $logsum;
	$proba[0] = exp(-$mylambda);
	for($j = 1;$j<$middle;$j++) 	               
        {
		$logsum =0;
		for($m=1;$m<=$j;$m++)                
        	{                
        		$logsum = $logsum + log($lambda)-log($m);                
        	}
		$proba[$j] = exp(-$lambda+$logsum);
	}  
	for($j = $middle;$j<100;$j++)
        {
		$proba[$j] = 0;
	}
	$sum =0;
	for($j=0;$j<$middle;$j++)
	{
		$sum = $sum + $proba[$j];
	}
	my $sum = 0;   
	for($m=0;$m<=$middle;$m++)
        {
       		$sum = $sum + $proba[$m]*$proba[$m];                
        }
	my $trace = $sum;
	my $result = -1;
	if($trace>(1 - $myprobcutoff))
	{
		$result = 1;
	}
	else
	{
		$sum=0;
        	for($k=1;$k<=$middle;$k++)# just a large number to garantee achieve the significance level.
        	{
                	for($m=0;$m<=$middle;$m++)
                	{
                        	$sum = $sum + $proba[$m]*$proba[$m+$k];
                	}
			if((2*$sum+$trace) > (1-$myprobcutoff))
			{	
				$result = $k;
				last;
			}
		} 
        }
	print "result= $result\n";
	$result;
}
