#!/usr/bin/perl -w
use POSIX;
# read in multiple realign results from ELAND, and combine their genomic locations.
# modified from seqread.pl.
# input files information in a file.
# Steve Qin 11/02/07
##### revised 04/13/08
# changed by JC at 08/11/2008
##### revised 09/24/2008 for control reads uninitialization

# how long a ChIP product fragment is suppose to be. 300 - 500.
#$fragmentWidth = 300;
#$windowSize = 25;

$argc = @ARGV;
$argc == 7 || 
die "Provide 
	1) species, 
	2) input all region, 
	3) & 4) all read file names (case and control), 
	5) & 6)  minimum and maximum DNA fragment length,
	7) result file name.\n";

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
$inputcontrolreadname = $ARGV[3];
open(INPUTCONTROLREADFILE, $inputcontrolreadname);
#$fragmentWidth = $ARGV[4];
$fragmin = $ARGV[4];
$fragmax = $ARGV[5];
$fragmentWidth = $fragmax;
$outputname = ">".$ARGV[6].".basecover.out"; 
open(OUTPUTFILE, $outputname);
$outputpeakname = ">".$ARGV[6].".hpeak.out";
open(OUTPUTPEAKFILE, $outputpeakname);

##### added 08/14/08
open(PARASFILE,$ARGV[6].".paras.txt");
$datline = <PARASFILE>;
$datline = <PARASFILE>;
$datline = <PARASFILE>;
$datline = <PARASFILE>;
chomp($datline);
@datas = split(" ", $datline);
$fgrate = $datas[0];
$datline = <PARASFILE>;
chomp($datline);
@datas = split(" ", $datline);
$bgrate = $datas[0];
$folds = $fgrate/$bgrate;
print "$bgrate\t$fgrate\t$folds\n";
close(PARASFILE);
my @added;

##### add weight for control reads 01/25/2009
open(RATIOFILE,"ratio.minus.txt")||die "can not read from ratio file\n";
$ratioln = <RATIOFILE>;
chomp($ratioln);
$ratio = $ratioln;
close(RATIOFILE);

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
###print "total = $total\n";
close(INPUTFILE);


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
	        for($k=0;$k<$countersizebychr[$chr];$k++)
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
	$skipstart = 0;#?????
##### changed 04/13/08 START
	$datline = <INPUTREADFILE>;
	if ($datline)
	{
		chomp($datline);
		@datas = split(" ",$datline);
		$readscount = scalar(@datas);
		if ($readscount > 1)
		{
			for($k=1;$k<$readscount;$k++)##### add first number as chromosome indicator.
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
#								if ($m==8 && $chr==0) {print "counter: $counter[$offset[$m] + $p]\n";}        
							}    
						}
					}
					elsif(($readstart >= $start[$m])&&($readend <= $end[$m]))
          				{#read contained in the region
          					$shift = $readstart - $start[$m];
						for($p=0;$p<$fragmentWidth;$p++)
            					{
       	    						$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] + $added[$p];
#							if ($m==8 && $chr==0) {print "counter: $counter[$offset[$m] + $shift + $p]\n";}
            					}
		      			}
					elsif(($end[$m] >= $readstart)&&($end[$m] < $readend))
					{#end in the middle of the read
						$shift = $readstart - $start[$m];
						$len = $end[$m] - $readstart +1;
            					for($p=0;$p<$len;$p++)  
            					{  
            						$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] + $added[$p];
#							if ($m==8 && $chr==0) {print "counter: $counter[$offset[$m] + $shift + $p]\n";}
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
			print "regions in chromosome", $chr+1," is $regioncountbychr[$chr] when go through input reads.\n";
		}	
		else
		{
			print "No reads found in chromosome ", $chr+1,".\n";
		}
	}
	else
	{
		print "No reads found in chromosome ", $chr+1,".\n";
	}
#	for($m=0;$m<$sizes[5];$m++)
#	{
#		print "$m is $counter[$offset[5] + $m]\n";
#	}
#	print "bobe done\n";


#####
#####
##### begin minus control 
#####
#####
##### add weight for control reads 01/25/2009
	$skipstart = 0;#reset for control.
##### changed 04/13/08 START
	$datline = <INPUTCONTROLREADFILE>;
	if ($datline)
	{
		chomp($datline);
		@datas = split(" ",$datline);
		$readscount = scalar(@datas);
		if ($readscount > 1)
		{
			for($k=1;$k<$readscount;$k++)##### skip the first number since it is chromosome indicator.
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
            							$counter[$offset[$m] + $p] = $counter[$offset[$m] + $p] - ($added[$startplace + $p] * $ratio);
#								if ($m==8 && $chr==0) {print "CTRL counter: $counter[$offset[$m] + $p]\n";}
            						}
						}
						else
						{#region goes out of the read
          						$stopplace = $readend - $start[$m] + 1; 
				    			for($p=0;$p<$stopplace;$p++)
            						{
            							$counter[$offset[$m] + $p] = $counter[$offset[$m] + $p] - ($added[$startplace + $p] * $ratio);
#								if ($m==8 && $chr==0) {print "CTRL counter: $counter[$offset[$m] + $p]\n";}        
            						}	    
						}
					}
					elsif(($readstart >= $start[$m])&&($readend <= $end[$m]))
        				{#read contained in the region
        					$shift = $readstart - $start[$m];
						for($p=0;$p<$fragmentWidth;$p++)
          					{
       	  						$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] - ($added[$p] * $ratio);
#							if ($m==8 && $chr==0) {print "CTRL counter: $counter[$offset[$m] + $shift + $p]\n";}
          					}
		    			}
					elsif(($end[$m] >= $readstart)&&($end[$m] < $readend))
	      				{#end in the middle of the read
						$shift = $readstart - $start[$m];
						$len = $end[$m] - $readstart +1;
        					for($p=0;$p<$len;$p++)  
          					{  
          						$counter[$offset[$m] + $shift + $p] = $counter[$offset[$m] + $shift + $p] - ($added[$p] * $ratio);
#							if ($m==8 && $chr==0) {print "CTRL counter: $counter[$offset[$m] + $shift + $p]\n";}
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
			print "regions in chromosome", $chr+1," is $regioncountbychr[$chr] when we go through control reads.\n";
		}
		else
		{
			print "NO reads in control chromosome", $chr+1,".\n";
		}
	}
	else
	{
		print "NO reads in control chromosome", $chr+1,".\n";
	}	

##### begin finding max.
	if(!$plusstrand)
	{
		print ("##### begin finding max.\n");
		$count = 0;
		for($k=0;$k<$regioncountbychr[$chr];$k++)
		{
			my @locations;
			my @saved;
			my $mymax = 0;
			$stretch = 0;
			$tie = 0;
			$first = 1;
			$count = 0;##### added 08/13/08
                	for($m=0;$m<$sizes[$k];$m++)
                	{
#				if ($k == 8) {print "$m \t $counter[$offset[$k] + $m] \t $mymax\n";}
				if($counter[$offset[$k] + $m] > $mymax)
				{
					$mymax = $counter[$offset[$k] + $m];
					$count = 1;
					$locations[0] = $m;
					$previous = $m; 
					$tie=0;
					$first =1;
				}
				elsif (($mymax>0)&&(abs($counter[$offset[$k] + $m] - $mymax)<0.001))
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
				else
				{
				### nothing happened :)
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
#				print "hema: region $k position: $m length $sizes[$k] mymax: $mymax\n";
#				if ($mymax == 0)
#				{
#					for($m=0;$m<$sizes[$k];$m++)
#					{
#						print "$m is $counter[$offset[$k] + $m]\n";
#					}
#					exit(0);
#				}
			}
#			print $maxpos;
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
##### changed 06/12/08           
			$maxpos++;
#12/28/08              	print OUTPUTPEAKFILE $chr+1,"\t$start[$k]\t$end[$k]\t$sizes[$k]\t$maxpos\t$mymax\t$logproba[$k]\t($logproba[$k]/$sizes[$k])\n";###\t$stretch\n";
			print OUTPUTPEAKFILE $chr+1,"\t$start[$k]\t$end[$k]\t$sizes[$k]\t$maxpos\t$mymax\t",$logproba[$k]/$sizes[$k],"\n";
##### changed 06/12/08    
		}
	}
	$plusstrand = !$plusstrand;
}
close(INPUTFILE);
close(INPUTREADFILE);
close(INPUTCONTROLREADFILE);
close(OUTPUTFILE);
close(OUTPUTPEAKFILE);

