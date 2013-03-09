#!/usr/bin/perl -w
use POSIX;
# read in multiple realign results from ELAND, and combine their genomic locations.
# modified from seqread.pl.
# input files information in a file.
# Steve Qin 11/02/07
##### revised 04/13/08
# changed by JC at 08/11/2008

# how long a ChIP product fragment is suppose to be. 300 - 500.
#$fragmentWidth = 300;
#$windowSize = 25;

$argc = @ARGV;
$argc == 6 || die "Provide species, input all region, all read file name, minimum, maximum DNA fragment length and result file name.\n";

$sp1 = $ARGV[0];
if ($sp1 eq "human" || "hg19") 
{
        $nchr = 24;
        $fold1 = "hg19";
}
elsif ($sp1 eq "hg18") 
{
        $nchr = 24;
        $fold1 = "hg18";
}
elsif ($sp1 eq "mouse" || "mm9") 
{
        $nchr = 21;
        $fold1 = "mm9";
}       
else
{
        print "species error!\n";
        exit(0); 
}       

$inputname = $ARGV[1]; 
open(INPUTFILE, $inputname);  
$inputreadname = $ARGV[2];
open(INPUTREADFILE, $inputreadname);  
$fragmin = $ARGV[3];
$fragmax = $ARGV[4];
$fragmentWidth = $fragmax;
$outputname = ">".$ARGV[5].".basecover.out"; 
open(OUTPUTFILE, $outputname);
$outputpeakname = ">".$ARGV[5].".hpeak.out";
open(OUTPUTPEAKFILE, $outputpeakname);
my @added;

##### first step, count number of regions in each chromosome.
my @regioncountbychr;
my @countersizebychr;
for($j=0;$j<$nchr;$j++)
{
        $regioncountbychr[$j] = 0;
	$countersizebychr[$j] = 0;
}
$total = 0;
while($datline = <INPUTFILE>)
{
        chomp($datline);
	if(length($datline) == 0)
	{
		next;
	}
        @datas = split(" ",$datline);
        $chr = $datas[0];
	$countersizebychr[$chr-1] = $countersizebychr[$chr-1] + $datas[3];
	$regioncountbychr[$chr -1] ++;
	$total = $total + $datas[3];
}
close(INPUTFILE);
#print "hema chr 0, $regioncountbychr[0]\n";
#exit(0);

##### second step, chromosome by chromosome.
my @start;
my @end;
my @offset;
my @counter;
my @logproba;
$plusstrand = 1;
open(INPUTFILE, $inputname);   
open(INPUTREADFILE, $inputreadname);

if ($fragmax == $fragmin)
{
	print "Note: fmin and fmax are set to be equal, so no extension will be performed.\n";
	$decayrate = 0;
}
else {$decayrate = 1/($fragmax-$fragmin);}

for($j=0;$j<($nchr*2);$j++)#48
{
	if($plusstrand)
	{#positive strand
		$chr = $j/2;            
		print "Chr ",$chr+1,"\n";
		for($k=0;$k<$fragmin;$k++)
		{
			$added[$k] = 1;
		}
		for($k=$fragmin;$k<$fragmax;$k++)
		{
			$added[$k] = 1-$decayrate*($k-$fragmin);
		}
		for($k=0;$k<$countersizebychr[$chr]+1;$k++)
		{# go through all counters
			$counter[$k] = 0;
		}
		for($k=0;$k<$regioncountbychr[$chr]+1;$k++)
		{# go through all regions          
			$offset[$k] = 0;
		}
		$current = 0;#to help track offset.
		for($k=0;$k<$regioncountbychr[$chr];$k++)
		{# go through all regions        
			$datline = <INPUTFILE>;
			chomp($datline);
			@datas = split(" ",$datline);
			$sizes[$k] = $datas[3];
			$start[$k] = $datas[1];
			$end[$k] = $datas[2];
			$logproba[$k] = $datas[5];
			$current = $current + $sizes[$k];
			$offset[$k+1] = $current;
		}
		print "all sizes = $current\n";
	}
	else
	{#odd $j, negative strand
		for($k=0;$k<$fragmin;$k++)
		{
			$added[$fragmax - $k - 1] = 1;
		}
		for($k=$fragmin;$k<$fragmax;$k++)
		{
			$added[$fragmax - $k - 1] = 1-$decayrate*($k-$fragmin);
		}      
	}
	$skipstart = 0;
##### changed 04/13/08 START
	$datline = <INPUTREADFILE>;
	if ($datline)
	{
		chomp($datline);
		@datas = split(" ",$datline);
		$readscount = scalar(@datas);
		if ($readscount > 1) ## make sure this line containn at least one read
		{
			for($k=1;$k<$readscount;$k++)
			{
				if($plusstrand)
				{#positive strand		
					$readstart = $datas[$k];
					$readend = $readstart + $fragmentWidth -1;
				}
				else
				{#negative strand
					$readend = - $datas[$k];
					$readstart = $readend - $fragmentWidth +1;
				}	    
				$newstartflag = 0;
#				print "hema $skipstart, $regioncountbychr[$chr]\n";
#				exit;
				for($m=$skipstart;$m<$regioncountbychr[$chr];$m++) 
				{
					if($end[$m] < $readstart)
					{# upstream of $readstart
						$newstart = $m;
						$newstartflag = 1;
						next;
					} 
					elsif($start[$m] >= $readend)
					{# downstream of $readstart
						last;
					}
					elsif(($start[$m] < $readend)&&($start[$m]>=$readstart))
					{#start in the middle of the read
						$startplace = $start[$m] - $readstart;	
						if($end[$m]<=$readend)
						{#region contained in read 
							for($p=0;$p<$sizes[$m];$p++)                           
							{
								$counter[$offset[$m] + $p] = $counter[$offset[$m] + $p] + $added[$startplace + $p];
							}
						}
						else
						{#region goes out of the read
							$stopplace = $readend - $start[$m] + 1; 
							for($p=0;$p<$stopplace;$p++)
							{
								$counter[$offset[$m] + $p] = $counter[$offset[$m] + $p] + $added[$startplace + $p];        
							}    
						}
					}
					elsif(($readstart >= $start[$m])&&($readend <= $end[$m]))
          				{#read contained in the region
						$shift = $readstart - $start[$m];
						for($p=0;$p<$fragmentWidth;$p++)
						{
							$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] + $added[$p];;
						}
					}
					elsif(($end[$m] >= $readstart)&&($end[$m] < $readend))
					{#end in the middle of the read
						$shift = $readstart - $start[$m];
						$len = $end[$m] - $readstart +1;
						for($p=0;$p<$len;$p++)  
						{  
							$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] + $added[$p];
						}
					}
					else
					{
						print "something is wrong, $j\t$k\t$m.\n";
						exit(0);
					}
				}
				if($newstartflag == 1)
				{
					$skipstart = $newstart; 		
				}
			}
		}
		else
		{
			print "no regions in chromosome", $chr+1,".\n";
		}
	}
	else
	{
		print "no regions in chromosome", $chr+1,".\n";
	}
##### begin finding max.
	if(!$plusstrand)
	{
		$count = 0;
		for($k=0;$k<$regioncountbychr[$chr];$k++)
		{
			my @locations;
			my $mymax = 0;
                	$stretch = 0;
                	$tie = 0;
                	$first = 1;
                	for($m=0;$m<$sizes[$k];$m++)
                	{
                        	if($counter[$offset[$k] + $m] > $mymax)
                        	{
                                	$mymax = $counter[$offset[$k] + $m];
                                	$count = 1;
                                	$locations[0] = $m;
                                	$previous = $m; 
                                	$tie=0;
                                	$first =1;
                        	}
                        	elsif(abs($counter[$offset[$k] + $m] - $mymax)<0.001)
                        	{
##### changed 06/12/08
                                	$tie=1;
                                	if(($m-$previous) ==1)#consecutive
                                	{
                                        	$locations[$count] = $m;
                                        	$count ++;
                                        	$previous = $m;
                                	}
                                	else
                                	{
                                        	if(($first==1)||($count > $stretch))
                                        	{
                                                	$stretch = $count;
                                                	for($n=0;$n<$count;$n++)
                                                	{
                                                        	$saved[$n] = $locations[$n];
                                                	}
                                        	}
                                        	$count = 1;
                                        	$first = 0;
                                        	$locations[0] = $m;
                                        	$previous = $m;
                                	}
##### changed 06/12/08
                        	}
                	}
##### changed 06/12/08    
			if($tie == 1)
			{
				if(($first==1)||($count > $stretch))
				{
					$stretch = $count;
					for($n=0;$n<$count;$n++)
					{
						$saved[$n] = $locations[$n];
					}
				}
				$maxpos = $saved[0] + ($stretch-1)/2;
			}
			else
			{
				$maxpos = $locations[0];
				$stretch = 1;
				$saved[0] = $maxpos;
			}
##### changed 06/12/08    
			print OUTPUTFILE $chr+1,"\t$start[$k]\t$end[$k]\t$sizes[$k]\t$mymax\t$maxpos\t$stretch\t";    
			for($m=0;$m<$stretch;$m++)
			{
				print OUTPUTFILE "$saved[$m] ";
			}
			print OUTPUTFILE "\n";
			for($m=0;$m<$sizes[$k];$m++)                      
			{
				print OUTPUTFILE "$counter[$offset[$k] + $m] ";
			}
			print OUTPUTFILE "\n";
#12/28/08		print OUTPUTPEAKFILE $chr+1,"\t$start[$k]\t$end[$k]\t$sizes[$k]\t$maxpos\t$mymax\t$logproba[$k]\t($logproba[$k]/$sizes[$k])\n";#$stretch\n";
                        print OUTPUTPEAKFILE $chr+1,"\t$start[$k]\t$end[$k]\t$sizes[$k]\t$maxpos\t$mymax\t",$logproba[$k]/$sizes[$k],"\n";
		}
	}
	$plusstrand = !$plusstrand;
}
close(INPUTFILE);
close(INPUTREADFILE);
close(OUTPUTFILE);
close(OUTPUTPEAKFILE);
