#!/usr/bin/perl 
#use strict;
# read in HMM results and determines the start and end of enriched regions.

##### modified 04/08/08
##### 1. change allregions from .total.txt information to hmmout total consecutive regions.
##### 2. add minimum region length threshold, 50 bp  by default.
$minregionlength = 2; #2 * 25 = 50 bps
#$windowSize = 25;

#changed by JC at 08/11/2008

$argc = @ARGV;
$argc == 7 || die "Provide 
	1) probability threshold, 
	2) species, 
	3) significance level, 
	4) window size, 
	5) parameter filename, 
	6) input filename,
	7) result file name.\n";

$sp2 = $ARGV[1];
if ($sp2 eq "hg18") 
{
        $nchr2 = 24;
        $fold2 = "hg18";
	$gcover = 0.9*3100000000;
}
elsif ($sp2 eq "human" || "hg19") 
{
        $nchr2 = 24;
        $fold2 = "hg19";
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
$parasfilename = $ARGV[4];
$inputfilename = $ARGV[5];
$outputfilename = $ARGV[6];

open(PARASFILE, $parasfilename);          
$datline = <PARASFILE>;
chomp($datline);
@datas = split(" ", $datline);
$fgread = $datas[0];
for($j=0;$j<4;$j++)
{
        $datline = <PARASFILE>;
}
chomp($datline);
@datas = split(" ", $datline);
$allregions = $datas[0];
$datline = <PARASFILE>;
chomp($datline);
@datas = split(" ", $datline);
$fragmentWidth = $datas[0];

close(PARASFILE);
$lambda = $fgread*$fragmentWidth*1000000/$gcover;          

$total = 0;
for($j=1;$j<($nchr2 + 1);$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	while($datline = <INPUTFILE>)
	{
		chomp($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
                $logproba = $datas[4];
                if($logproba < log(0.00001)) 
		{
 			if($first == 1)
			{
				$start = $order;
				$end = $order;
				$first = 0;
				$length = 1;
				$count = 1;
			}
			else
			{
				if($order == ($end +1))
				{# continue
					$end ++;
					$length ++;
				}
				else
				{# start a new region
					$start = $order;
					$end = $order;
        	                        $length = 1;
                	                $count ++;
				}
			}
		}
	}
	$total = $total + $count;
}
close(INPUTFILE);
print "there are $total regions.\n";
$allregions = $total;

open(OUTPUTFILE, ">$outputfilename.allregions.txt");
$total = 0;
for($j=1;$j<($nchr2 + 1);$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	my $maxprobacount = 0;
	while($datline = <INPUTFILE>)
	{
		chomp($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
		$height = $datas[3];
                $logproba = $datas[4];
		if($logproba < log(0.00001))
		{
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
					if($maxprobacount >= $minregionlength)
					{
                                                print OUTPUTFILE "$j\t",$start*$windowSize+1,"\t",($end+1)*$windowSize,"\t",$length*$windowSize,"\t$maxheight\t$sumlogproba\n";
					}
					else
                                        {         
                                        }
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
				}
			}
		}
	}
	if($maxprobacount >= $minregionlength)
        {	
                print OUTPUTFILE "$j\t",$start*$windowSize+1,"\t",($end+1)*$windowSize,"\t",$length*$windowSize,"\t$maxheight\t$sumlogproba\n";
	}
	$total = $total + $count;
}
close(INPUTFILE);
close(OUTPUTFILE);
