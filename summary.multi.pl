#!/usr/bin/perl
use POSIX;
use Cwd qw(realpath);
use File::Basename; 
# Steve Qin 11/02/07
# Changed by JC 08/11/2008
#            JC 08/30/2008


#$fragmentWidth = 300;
#$readlength = 25;

my $script_path=realpath($0);
my($ef, $dir) = fileparse($script_path);

$argc = @ARGV;
$argc == 8 || die "Provide 
	1) species info,
	2) input files name that include filename(s) of uniquely-mapped reads (*_eland_result.txt, *_export.txt), 
	3) input files name that include filename(s) of multiple-mapped reads (*_eland_multi.txt), 
	4) format,
	5) read length;
	6) fragment size, 
	7) significant level
	8) result file name.\n";

$sp1 = $ARGV[0];
$innamename = $ARGV[1];
open(INNAMEFILE, $innamename) || die "Can not open $innamename\n";

$inmultiname = $ARGV[2];

$format = $ARGV[3];
$format=lc($format);
$readlength=$ARGV[4];
$fragmentWidth = $ARGV[5];
$siglevel = $ARGV[6];
$sname = $ARGV[7];
open(OUTFILE1, ">".$sname.".location.txt") || die "$!\n";
open(OUTFILE2, ">".$sname.".dis.txt") || die "$!\n";
open(OUTFILE3, ">".$sname.".reads.txt") || die "$!\n";
open(OUTSUMFILE,">".$sname.".total.txt") || die "$!\n";
open(OUTPARFILE,">".$sname.".paras.txt") || die "$!\n";

##### 1. read in unquely-mapped reads:
#####
my @filenames;
while($datline = <INNAMEFILE>)
{
	chomp($datline);
	push(@filenames,$datline) if($datline);
}
print "There are ", scalar(@filenames), " unique-mapped read files.\n";
close(INNAMEFILE);

##### 2. read in multiple-mapped reads.
#####
open(INMULTINAMEFILE, $inmultiname) || die "$!\n";
@multifilenames = ();
while($datline = <INMULTINAMEFILE>)
{         
	chomp($datline);
	if ($datline =~/(.*)\<format\>(.*)\<\/format\>/){$datline=$1;}
	push(@multifilenames,$datline) if($datline);
}
close(INMULTINAMEFILE);
print "There are ",scalar(@multifilenames)," multiple-mapped read files.\n";

##### 3. innitialize species specific parameters.
#####
if ($sp1 eq "human" || "hg19") 
{
	$nchr = 24;
	$gcover = 0.9*3100000000;
}
elsif ($sp1 eq "hg18") 
{
	$nchr = 24;
	$gcover = 0.9*3100000000;
}
elsif ($sp1 eq "mouse" || "mm9") 
{
	$nchr = 21;
	$gcover = 0.7*2655000000;
}
else 
{
	print "species error!\n";
	exit(0);
}

##### 4. initialize hits counts
#####
for($i = 0;$i < $nchr; $i++)
{
	$hitscount[$i] = 0;
	$posreadscount[$i] = 0;
	$negreadscount[$i] = 0;
}

##### 5. read in uniquely-mapped read data.      
#####
$totalcount = 0;
$count = 0;
foreach my $currentline (@filenames)
{
	if ($currentline =~/(.*)\<format\>(.*)\<\/format\>/)
	{
		$tmpformat = lc($2);
		$tmpname = $1;
	}
	else
	{
		$tmpformat = $format;
		$tmpname = $currentline;
	}
	&realignfileread($tmpname,$tmpformat,*totalcount, *count,*hits, *hitscount,*directedreads,*posreadscount, *negreadscount);##### newly added 04/13/08
	print "total read count = $totalcount.\n";
	print "unique read count according to the current datbase = $count.\n";
}
$uniquereadtotal = $count;
$totalread = $count;

##### 6. sort reads for each chromosome
#####
foreach (0..($nchr-1)) {
	my @pos=@{$hits[$_]};
	@{$hits[$_]}=sort {$a<=>$b} @pos if(scalar(@pos)>1);
}

##### 7. read in sorted multi-mapped read data, del 11/18/2008
##### JC later decide not to use the filter

##### 8. read in multi-mapped read data.
#####
########added by jianjun: Assign multi-mapping reads######
########changed by JC 11/22/20 take the whole multifilename array once #################
#print "$totalcount\t$count\n";
$multitotalcount = 0;
&multifileread($sname,*multifilenames, *count, *multitotalcount, *hits, *hitscount, *directedreads, *posreadscount, *negreadscount, *sorthits);
####sort reads for each chromosome again after adding new positions from multi-mapping reads3##
print "multi=$multitotalcount\ttotal count (unique+multi) according to our current database = $count\n";
$totalread = $count;
foreach (0..($nchr-1)) {
	my @pos=@{$hits[$_]};
	@{$hits[$_]}=sort {$a<=>$b} @pos if(scalar(@pos)>1);
}

$count = 0;
for($i = 0;$i < $nchr; $i++)
{
        print OUTFILE3 "chr",$i+1,"\t";
	if($hitscount[$i]==0)
        {
		print OUTFILE3 "\n";
                next;
        }
	elsif($hitscount[$i]==1)
        {
                $count ++;
		next;
        }
        my @sortpos=@{$hits[$i]};
##### newly added 07/13/08 START
        my @directedstart;
        for($j=0;$j<$hitscount[$i];$j++)
        {
                $directedstart[$j] = $directedreads[$i][$j];
        }
        @sortdirectedstart = sort {$a <=> $b} @directedstart;
	for($j=0;$j<$posreadscount[$i];$j++)
        {
		print OUTFILE3 $sortdirectedstart[$negreadscount[$i] + $j],"\t";
        }
        print OUTFILE3 "\nchr",$i+1,"\t";
        for($j=0;$j<$negreadscount[$i];$j++)
        {
		print OUTFILE3 $sortdirectedstart[$negreadscount[$i] - $j - 1],"\t";
        }
        print OUTFILE3 "\n";
##### newly added 07/13/08 END
	for($j=1;$j<$hitscount[$i];$j++)
        {
                if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
                {#extension continues
                        next;
                }
                else
                {#extension stops
                        $count ++;
                }
        }       
		$count ++;
}
close(OUTFILE3);
$allregion = $count;
print "total region = $count.\n";

#$order = 1;
$totalpeaks = 0;
$peakcover = 0;
$totalcover = 0;
$totalinpeakreads = 0;
my @peakwidth; 
my @inpeakreadscount;
my $keyratio = $totalread /$gcover;
for($i = 0;$i < $nchr; $i++)
{
        print OUTPUTFILE "chr",$i+1,"\t";
	if($hitscount[$i]==0)
        {
		next;
	}
	my @sortpos=@{$hits[$i]};
	$count = 1 ;
	$start = $sortpos[0];
	print OUTFILE1 $i+1, "\t$start\t";
	if($hitscount[$i]==1)
	{
		$end = $sortpos[0] + $fragmentWidth -1;                 
	        $regionwid = $end - $start +1;
        	$totalcover = $totalcover + $regionwid; 
		print OUTFILE1 "\n";
		next;
	}

	for($j=1;$j<$hitscount[$i];$j++)
	{
		$dis = $sortpos[$j] - $sortpos[$j-1];
		print OUTFILE2 "$dis\n" if($dis <1000);
		if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
		{#extension continues
	                print OUTFILE1 "$sortpos[$j]\t";
			$count ++;
			next;
		}
		else 
		{#extension stops
			print OUTFILE1 "\n";
			$end = $sortpos[$j-1] + $fragmentWidth -1;
			$regionwid = $end - $start +1;
			$totalcover = $totalcover + $regionwid;
			$lambda = $keyratio * $regionwid;
			$sum = 1;
			for($k=1;$k<$count;$k++)
			{
				$logsum = 0;
				for($m=1;$m<=$k;$m++)
                        	{
					$logsum = $logsum + log($lambda)-log($m);
				}
				$sum = $sum + exp($logsum);
			}
			$prob = 1 - exp(-$lambda)* $sum;
			if($prob < 0.0000000001)
			{
				$prob = 0;
			}
                        if($prob < $siglevel/$allregion)
			{
                	        #print OUTFILE2 "Chr", $i + 1,"\t",$start,"\t", $end, "\t",$regionwid,"\t", $count,"\t", "\n";	
				$totalinpeakreads = $totalinpeakreads + $count;
				$peakwidth[$totalpeaks] = $regionwid;
                                $inpeakreadscount[$totalpeaks] = $count;
				$totalpeaks ++;
				$peakcover = $peakcover + $regionwid;
			}
        	        $start = $sortpos[$j];
		        print OUTFILE1 $i + 1,"\t$start\t";
                        $count = 1;
		}
	} 
##### start of added 04/20/08
        $end = $sortpos[$hitscount[$i]-1] + $fragmentWidth -1;
        $regionwid = $end - $start +1;
	$totalcover = $totalcover + $regionwid;
##### 07/04/08        $lambda = $keyratio * $regionwid;
##### 07/13/08 don't know why the above deletion.
        $lambda = $keyratio * $regionwid;
	$sum = 1;
        for($k=1;$k<$count;$k++)
        {
                $logsum = 0;
                for($m=1;$m<=$k;$m++)
                {
                        $logsum = $logsum + log($lambda)-log($m);
                }
                $sum = $sum + exp($logsum);
        }
        $prob = 1 - exp(-$lambda)* $sum;
        if($prob < 0.0000000001)
        {
                $prob = 0;
        }
        #if($count > $threshold)
        if($prob < $siglevel/$allregion)
        {
                $totalinpeakreads = $totalinpeakreads + $count;
                $peakwidth[$totalpeaks] = $regionwid;
                $inpeakreadscount[$totalpeaks] = $count;
                $totalpeaks ++;
                $peakcover = $peakcover + $regionwid;
         }
##### end of added 04/20/08 
   print OUTFILE1 "\n";
}
$avecover = $totalcover /$totalread;
$totalcover = $totalcover /1000000;
print "total coverage is $totalcover MB.\n";
print "average coverage per read is $avecover KB.\n";
close(OUTFILE1);
close(OUTFILE2);

##### determine the hit count threshold per chromosome.
$lam = $totalread*$fragmentWidth /$gcover;
$probacut = $siglevel/$allregion;
print "lambda=$lam\tprobacut=$probacut\n";
$hitlevel = 0;
$proba =1;
$sum=1;
while($proba > $probacut)
{
	$hitlevel ++;
	$logsum = 0;
	for($m=1;$m<=$hitlevel;$m++)
	{            
		$logsum = $logsum + log($lam)-log($m);
	}
	$sum = $sum + exp($logsum);
	$proba = 1 - exp(-$lam)* $sum;
	print "proba=$proba\n";
}                 
print "totalread: $totalread\tallregion: $allregion\thitlevel: $hitlevel\n";
@sortpeakwidth = sort {$a <=> $b} @peakwidth;
$medianwidth = $sortpeakwidth[floor(0.5*$totalpeaks)];

print "total number of reads = $totalread.\n";
print "total number of uniquely-aligned reads = $uniquereadtotal. (",100*$uniquereadtotal/$totalread," %)\n";
print "total number of multiple-aligned reads = $multitotalcount. (",100*$multitotalcount/$totalread," %)\n";
print OUTSUMFILE "Total sequenced reads (millions)\t",$totalread/1000000,"\n";
close(OUTSUMFILE); 

print OUTPARFILE $totalread/1000000,"\n";
print OUTPARFILE $totalinpeakreads/1000000,"\n";
print OUTPARFILE $peakcover/1000000,"\n";
print OUTPARFILE "$medianwidth\n";
print OUTPARFILE "$totalpeaks\n";
print OUTPARFILE $fragmentWidth,"\n";#to accout for linear descreasing at the tail.
print OUTPARFILE $hitlevel,"\n";
close(OUTPARFILE);

sub parse_match 
{
        local (*sub_multireads,*sub_locs,$fw,$readName,$str,$mode)=@_;
        my @matches=split /\,/,$str;
        foreach my $match (@matches) 
	{
                if($match=~/^hs_ref_chr(\d+|X|Y)\.fa\:(\d+)([FR])$mode/) 
		{
                        my $index=0;
                        my ($chr,$pos,$direct)=($1,$2,$3);
                        $index=scalar(@{$sub_multireads{$readName}}) if(defined $sub_multireads{$readName});
                        $chr=$nchr-1 if($chr eq "X");
                        $chr=$nchr if($chr eq "Y");
                        $sub_multireads{$readName}[$index]=[($chr,$pos,$direct)];#$1:chr,$2:pos,$3:direction/strand
                        if($direct eq "R") 
			{
                                $pos=$pos-$fw+1;
                                $sub_locs{$chr}{$pos}++; #chr,pos
#				print "pos: $pos\n";
#				exit;
                        } 
			else 
			{
                                $sub_locs{$chr}{$pos}++; #chr,pos
                        }
                }
        }
}

sub multifileread 
{
	local ($sname,*files, *subtotalcount, *subcount, *subhits, *subhitscount, *subdirectedreads, *subposreadscount,*subnegreadscount, *subsorthits) =@_;
	my %multireads=%locs=%uniqreads_per_pos=();
	print "Total ",scalar(@files)," multi-mapping read files found\n";
	foreach my $file (@files)
	{
        	open(INPUT,"<$file") || die "Can not open file $file\n";
        	while(my $line=<INPUT>) 
		{
			chomp($line);
			my @ele=split /\s+/,$line;
			next if($ele[3] eq "" || @ele<4);
			#assign chromosome # to reads where multiple reads are located at same chromosome
			if($ele[3]=~/\,\d+[FR]\d/) 
			{
				my @locs=split /\,/,$ele[3];
				my $prefix;
				for(my $i=0;$i<scalar(@locs);$i++) 
				{
					if($locs[$i]=~/^(.+\:)\d+[FR]/) 
					{
						$prefix=$1;
					} 
					elsif($locs[$i]=~/^\d+[FR]/) 
					{
						$locs[$i]=$prefix.$locs[$i];
					}
				}
				$ele[3]=join(",",@locs);
			}
			if($ele[2]=~/^(10|[2-9])\:/) 
			{
				&parse_match(\%multireads,\%locs,$fragmentWidth,$ele[0],$ele[3],0);
			}
			elsif($ele[2]=~/^0\:(10|[2-9])\:/) 
			{
				&parse_match(\%multireads,\%locs,$fragmentWidth,$ele[0],$ele[3],1);
			}
			elsif($ele[2]=~/^0\:0\:(10|[2-9])/) 
			{
				&parse_match(\%multireads,\%locs,$fragmentWidth,$ele[0],$ele[3],2);
			}
		}
	}	
	foreach my $chr (sort {$a<=>$b} keys %locs) 
	{
		my @pos_array=sort {$a<=>$b} keys %{$locs{$chr}};
		my $j=0;
		my @arr=@{$subhits[$chr-1]};
		foreach my $pos (@pos_array)
		{
			my $start=$pos-$fragmentWidth+1;
			my $end=$pos+$fragmentWidth;
			while($j<@arr && $arr[$j]<$start) 
			{
				$j++;
			}
			for(my $i=$j;$i<@arr;$i++) 
			{
				last if($arr[$i]>$end);
				$uniqreads_per_pos{$chr}{$pos}++;
			}
		}	
	}
	%locs=();
	foreach my $read (sort keys %multireads) 
	{
		my @matches=@{$multireads{$read}};
		my $sum=0;
		my @scores=();
		for(my $i=0;$i<@matches;$i++) 
		{
			my ($chr, $pos,$direct)=@{$matches[$i]};
#			print "$chr\t$pos\t$direct\n";
			if($direct eq "F" || $direct eq "+") 
			{
				if($uniqreads_per_pos{$chr}{$pos}>0) 
				{
					$sum+=$uniqreads_per_pos{$chr}{$pos};
					$scores[$i]=$uniqreads_per_pos{$chr}{$pos};
				}
			} 
			else 
			{
				$pos=$pos - $fragmentWidth + $readlength;
				if($uniqreads_per_pos{$chr}{$pos}>0) 
				{
					$sum+=$uniqreads_per_pos{$chr}{$pos};
					$scores[$i]=$uniqreads_per_pos{$chr}{$pos};
				}
			}
		}
		next if($sum<1);
		for(my $i=0;$i<@matches;$i++) 
		{
			if($scores[$i]>0) 
			{
				my ($chr, $pos,$direct)=@{$matches[$i]};
				$locs{$chr}{$pos}{$direct}+=$scores[$i]/$sum;
			}
		}
	}
	undef %multireads;
#	print "bobe see subtotalcount = $subtotalcount\n";
	open(OUT,">".$sname."_multiread.bed") || die "Can not open output bed file:$!\n";
	foreach my $chr (sort {$a<=>$b} keys %locs) 
	{
		foreach my $pos (sort {$a<=>$b} keys %{$locs{$chr}}) 
		{
			foreach my $direct (sort keys %{$locs{$chr}{$pos}}) 
			{
				my $val=$locs{$chr}{$pos}{$direct};
				print OUT "$chr\t$pos\t$direct\t$val\n";
				$val=int($val+0.5);
				$subtotalcount += $val;
				if($direct eq "F" || $direct eq "+")
				{
					$subposreadscount[$chr-1] += $val;
					for(my $j=0;$j<$val;$j++)
					{
						$subhits[$chr-1][$subhitscount[$chr-1]] = $pos;
						$subdirectedreads[$chr-1][$subhitscount[$chr-1]] = $pos;
						$subhitscount[$chr-1] ++;
					}
				} 
				elsif($direct eq "R" || $direct eq "-") 
				{
					$subnegreadscount[$chr-1] += $val;
					for(my $j=0;$j<$val;$j++)#$loading
					{
						$subhits[$chr-1][$subhitscount[$chr-1]] = $pos + $readlength - $fragmentWidth;
						$subdirectedreads[$chr-1][$subhitscount[$chr-1]] = -($pos+$readlength-1);
						$subhitscount[$chr-1] ++;
					}
				}
				$subcount += $val;
			}
		}
	}
	close(OUT);
	print "subcount for multiple hits= $subcount.\n";
}

sub realignfileread {
         local ($subcurrentfilename, $subformat, *subtotalcount, *subcount, *subhits, *subhitscount, *subdirectedreads, *subposreadscount, *subnegreadscount) =@_;
         if($subcurrentfilename =~ /gz$/)
         {
                 system("gunzip ".$subcurrentfilename);
                 $namelength = length($subcurrentfilename);
                 $newname = substr($subcurrentfilename,0,$namelength-3);
                 open(INREADFILE, $newname) || die "$!\n";
         }
         elsif($subcurrentfilename =~ /zip$/)
         {
                 system("unzip ".$subcurrentfilename);
                 $newname = substr($subcurrentfilename,0,$namelength-4);
                 open(INREADFILE, $newname.".txt") || die "$!\n";
         }
         else
         {
                 open(INREADFILE, $subcurrentfilename) || die "$!\n";
         }
	while($datline = <INREADFILE>)
	{
		my $chr=$posit=$direct=$loc=0;
		chomp($datline);
		my @datas=split /[\s\t]+/, $datline;
		$subtotalcount ++;
		if($subformat eq "eland")
		{
			if(($datline=~ /chr/)&&($datline =~/\.fa/))
   			{
				$posit = $datas[7];
        			$direct = $datas[8];
        			$loc=$datas[6];
			}
		}	 
		elsif ($subformat=~/custom\[(.*)\]/)
		{
			$tempo = $1;
			if(($datline=~ /chr/)&&($datline =~/\.fa/))
			{
				my @infos=split /\,\s*/,$tempo;
				$loc=$datas[$infos[0]-1];
				$posit=$datas[$infos[1]-1]+1;
				$direct=$datas[$infos[2]-1];
			}
		} 
		elsif ($subformat eq "bed") 
		{
			$chr =$2 if($datas[0]=~/^(chr)?(\d+|X|Y)$/);
			$posit = $datas[1]+1;
			$direct = $datas[5];
		}
		$chr=$1 if($loc=~/chr(\d+|X|Y)\.fa/ || $loc=~/chromosome\.(\d+|X|Y)\.fa/);
		next if(!$chr);
		$chr=($nchr-1) if($chr eq "X");
		$chr=($nchr) if($chr eq "Y");
		if($direct eq "+" || $direct eq "F")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit;
			$subdirectedreads[$chr-1][$subhitscount[$chr-1]] = $posit;
			$subposreadscount[$chr-1] ++;
		}
		elsif($direct eq "-" || $direct eq "R")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit + $readlength - $fragmentWidth;
			$subdirectedreads[$chr-1][$subhitscount[$chr-1]] = -($posit + $readlength - 1);
			$subnegreadscount[$chr-1] ++;
		}
		$subhitscount[$chr-1] ++;
		$subcount++;
	}
close(INREADFILE);
print "subcount = $subcount.\n";
}

