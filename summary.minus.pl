#!/usr/bin/perl
use POSIX;
# Steve Qin 11/02/07
# changed by JC at 08/11/2008

$argc = @ARGV;

#$readlength = 25;
$argc == 10 || die "Provide 
	1) species info,
	2) multi-tag for treated file,
	3) multi-tag for control file, 
	4) treated files name that include filename(s) of uniquely-mapped reads, 
	5) control files name that include filename(s) of uniqulye-mapped reads, 
	6) format,
	7) read length,
	8) fragment size, 
	9) significant level
	10) result file name.\n";

$sp1 = $ARGV[0];
if ($sp1 eq "human" || "hg19") 
{
	$nchr = 24;
	$fold1 = "hg19";
	$gcover = 0.9*3100000000;
}
elsif ($sp1 eq "hg18") 
{
	$nchr = 24;
	$fold1 = "hg18";
	$gcover = 0.9*3100000000;
}
elsif ($sp1 eq "mouse" || "mm9") 
{
	$nchr = 21;
	$fold1 = "mm9";
	$gcover = 0.7*2655000000;
}       
else
{
	print "species error!\n";
	exit(0); 
}

$tmulti = $ARGV[1];
$cmulti = $ARGV[2];
$sname = $ARGV[9];       
$innamename = $ARGV[3];
if ($tmulti == 1) 
{$innamename = $sname."_treated.unique";}
open(INNAMEFILE, $innamename)|| die "Can not open $innamename\n"; 
$incontrolnamename = $ARGV[4];
if ($cmulti == 1) 
{$incontrolnamename = $sname."_control.unique";}
open(INCONTROLNAMEFILE, $incontrolnamename)|| die "Can not open $incontrolnamename\n";
$format = $ARGV[5];
$readlength = $ARGV[6];
$fragmentWidth = $ARGV[7];
$siglevel = $ARGV[8];
open(OUTPARFILE,">".$sname.".paras.txt");

$count = 0;
my @filenames;
while($datline = <INNAMEFILE>)
{
	chomp($datline);
	if(length($datline) ==0){next;}
	$filenames[$count] = $datline;
	$count ++;
}
$filecount = $count;
print "There are $filecount case files.\n";
for($i = 0;$i < $nchr; $i++)
{
	$hitscount[$i] = 0;
	$multihitscount[$i] = 0;
}
$totalcount = 0;
$count = 0;
for($i=0;$i<$filecount;$i++)
{
	$currentline = $filenames[$i];
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
	&realignfileread($tmpname, $tmpformat,*totalcount, *count,*hits, *hitscount);
	print "total count = $totalcount.\n";
	print "unique read count = $count.\n";
}
$uniquereadtotal = $count;
#### A. sort reads for each chromosome
####
foreach (0..($nchr-1))
{
	my @pos=@{$hits[$_]};
	@{$hits[$_]}=sort {$a<=>$b} @pos if(scalar(@pos)>1);	
}
#### B. read in multi_treated.bed.
####
if (tmulti == 1)
{
	$multitotalread = 0;
	$inmultinamename = $sname."_treated_multiread.bed";
	&multifileread($inmultinamename,*multitotalread, *count, *multihits, *multihitscount, *loading);
	print "total multi read count for treatedfile = $multitotalread.\n";
	print "total read count for treatedfile (unique + multi) = $count.\n";
}
$totalread = $count;

$count = 0;
my @controlfilenames;
while($datline = <INCONTROLNAMEFILE>)
{
	chomp($datline);
	if(length($datline) ==0)
	{
		next;
	}
	$controlfilenames[$count] = $datline;
	$count ++;
}
$filecount = $count;     
print "There are $filecount control files.\n";
for($i = 0;$i < $nchr; $i++)
{
	$controlhitscount[$i] = 0;
	$multictrlhitscount[$i] = 0;
}
$totalcontrolcount = 0;
$controlcount = 0;
for($i=0;$i<$filecount;$i++)
{
	$currentline = $controlfilenames[$i];
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
	&realignfileread($tmpname,$tmpformat,*totalcontrolcount, *controlcount,*controlhits, *controlhitscount);
	print "total control count = $totalcontrolcount.\n";
	print "control count = $controlcount.\n";
}
$uniquectrlreadtotal = $controlcount;

#$ratio = $totalread/$totalcontrolread;
#print "$totalread\t$totalcontrolread\t$ratio\n";

#### B. read in multi_control.bed.
####
$multicontrolcount = 0;
if (cmulti == 1)
{
	$inmultictrlname = $sname."_control_multiread.bed";
	&multifileread($inmultictrlname,*multicontrolcount, *controlcount,*multictrlhits, *multictrlhitscount, *ctrlloading);
	print "total multi-mapped read count for control file = $multicontrolcount.\n";
	print "total reads count in control file (unique + multi) = $controlcount.\n";
}
$totalcontrolread = $controlcount;

$ratio = $totalread/$totalcontrolread;
print "$totalread\t$totalcontrolread\t$ratio\n";

#########JC 01/25/2009
open(RATIOFILE,">ratio.minus.txt")||die "can not write into ratio file\n";
print RATIOFILE "$ratio\n";
close(RATIONFILE);

$count = 0;
for($i = 0;$i < $nchr; $i++)
{
	my @pos;
	for($j=0;$j<$hitscount[$i];$j++)
	{
		$pos[$j] = $hits[$i][$j];
	}
	if (tmulti == 1)
	{
		for ($j0=0;$j0<$multihitscount[$i];$j0++)
		{
			$pos[$j+$j0+1] = $multihits[$i][$j+$j0+1];
		}	
	}
	@sortpos = sort {$a <=> $b} @pos;
	
	for($j=1;$j<scalar(@sortpos);$j++)
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
$allregions = $count;
print "total region = $count.\n";

##### determine the hit count threshold per chromosome.
$lam = $totalread*$fragmentWidth /$gcover;
$probacut = $siglevel/$allregions;
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
print "total reads: $totalread\ttotal regions: $allregions\thit level: $hitlevel\n";

$probcutoff = $siglevel/$allregions;
print "$allregions\n";
$order = 1;
$totalpeaks = 0;
$peakcover = 0;
$totalinpeakreads = 0;
$totalinpeakcontrolreads = 0;
$totalinpeakdif = 0;
my @peakwidth; 
my @inpeakdifcount;
$badcount = 0;
for($i = 0;$i < $nchr; $i++)
{
	print "chr",$i+1,"\n";
	my @pos;
	for($j=0;$j<$hitscount[$i];$j++)
	{
		$pos[$j] = $hits[$i][$j];
	}
	if (tmulti == 1) 
	{
		for ($j0=0;$j0<$multihitscount[$i];$j0++)
		{
			$pos[$j+$j0+1] = $multihits[$i][$j0];		
		}
	}
	@sortpos = sort {$a <=> $b} @pos;

	my @controlpos;
	for($j=0;$j<$controlhitscount[$i];$j++)
	{
		$controlpos[$j] = $controlhits[$i][$j];
	}
	if (cmulti == 1)
	{
		for ($j0=0;$j0<$multictrlhitscount[$i];$j0++)
		{
			$controlpos[$j+$j0+1] = $multictrlhits[$i][$j0];
		}	
	}
	@sortcontrolpos = sort {$a <=> $b} @controlpos;

	$count = 1;
	$start = $sortpos[0];
	my $j;
	$skipstart = 0;
	for($j=1;$j<scalar(@sortpos);$j++) 
	{ 
		if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
		{#extension continues
			$count ++;
			next;
		}
		else 
		{#extension stops
### need to change? 05/11/2009
			$end = $sortpos[$j-1] + $fragmentWidth - 1;
			$regionwid = $end - $start + 1;
#			print "reg: $regionwid $end $start;\n";
			$matchcontrolcount = &checkcontrol($start,$end,scalar(@sortcontrolpos),*skipstart, *sortcontrolpos);
			$baseline = $totalread /$gcover;
			$dif = $count - $ratio * $matchcontrolcount;
			if($dif > 0)
			{
				if($dif > 100)
				{
					$outcome = 1;
				}
				else
				{
					$outcome = &decide($baseline, $regionwid, $dif, $probcutoff);
				}
				if($outcome ==1)
				{
					$totalinpeakreads = $totalinpeakreads + $count;
					$totalinpeakcontrolreads = $totalinpeakcontrolreads + $matchcontrolcount;
					$totalinpeakdif = $totalinpeakdif + $dif;
					$peakwidth[$totalpeaks] = $regionwid;
					$inpeakdifcount[$totalpeaks] = $dif;
					$totalpeaks ++;
					$peakcover = $peakcover + $regionwid;
				}
			}
			$start = $sortpos[$j];
			$order ++;
			$count = 1;
		}
	} 
}
print "badcount = $badcount.\n";

@sortpeakwidth = sort {$a <=> $b} @peakwidth;
$medianPeakwidth = $sortpeakwidth[floor(0.5*$totalpeaks)];

print "total number of uniquely-mapped reads = $totalcount.\n";
print "total number of successfully aligned uniquely-mapped reads = $uniquereadtotal. (",100*$uniquereadtotal/$totalcount," %)\n";
print "total number of successfully aligned reads = $totalread, uniquely mapped = $uniquereadtotal (",100*$uniquereadtotal/$totalread," %), multi-mapped = $multitotalread (",100*$multitotalread/$totalread," %)\n";
print "total number of uniquely-mapped control reads = $totalcontrolcount.\n";
print "total number of successfully aligned uniquely control reads = $uniquectrlreadtotal. (",100*$uniquectrlreadtotal/$totalcontrolcount," %)\n";
print "In control file, total number of successfully aligned control reads = $totalcontrolread, uniquely mapped = $uniquectrlreadtotal (",100*$uniquectrlreadtotal/$totalcontrolread," %), multi-mapped = $multictrltotalread (",100*$multictrltotalread/$totalcontrolread," %)\n";

print OUTPARFILE $peakcover/1000000,"\n";    
print OUTPARFILE $totalinpeakreads,"\n";
print OUTPARFILE $totalinpeakcontrolreads,"\n";    
print OUTPARFILE "$totalread\n";
print OUTPARFILE "$totalcontrolread\n";
print OUTPARFILE $fragmentWidth,"\n";# to account for linear descreasing at the tail.
print OUTPARFILE "$medianPeakwidth\n";
print OUTPARFILE $hitlevel,"\n";   
close(OUTPARFILE);

sub checkcontrol
{
	local($localstart, $localend,$localchrsize, *localskipstart,*localsortcontrolpos) = @_;
	my $j;
	my $count = 0;
	my $newstart = $localskipstart;
	for($j=$localskipstart;$j<$localchrsize;$j++)
	{
		if(($localstart - $localsortcontrolpos[$j]) > $fragmentWidth)
		{
			$newstart = $j;
			next;
		}
		if($localend < $localsortcontrolpos[$j])
		{
			last;
		}
		if((($localstart - $localsortcontrolpos[$j]) <= $fragmentWidth)&&($localend >= $localsortcontrolpos[$j]))
		{ 
			$count ++;
		}
	}
	$localskipstart = $newstart;
	$count;
}

sub realignfileread
{
	local ($subcurrentfilename, $subformat, *subtotalcount, *subcount, *subhits, *subhitscount) = @_; 
	if($subcurrentfilename =~ /gz$/)
	{
		system("gunzip ".$subcurrentfilename);
		$namelength = length($subcurrentfilename);		
		$newname = substr($subcurrentfilename,0,$namelength-3); 
		open(INREADFILE, $newname);
	}
	elsif($subcurrentfilename =~ /zip$/)
	{
		system("unzip ".$subcurrentfilename);
		$newname = substr($subcurrentfilename,0,$namelength-4);
		open(INREADFILE, $newname.".txt");
	}
	else
	{
		print "$subcurrentfilename.\n";
		open(INREADFILE, $subcurrentfilename);		
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
			$direct = $datas[-1];
		}
		$chr=$1 if($loc=~/chr(\d+|X|Y)\.fa/ || $loc=~/chromosome\.(\d+|X|Y)\.fa/);
		next if(!$chr);
		$chr=($nchr-1) if($chr eq "X");
		$chr=$nchr if($chr eq "Y");
		if($direct eq "+" || $direct eq "F")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit;
		}
		elsif($direct eq "-" || $direct eq "R")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit + $readlength - $fragmentWidth;
##### newly added 04/13/08 START
##### newly added 04/13/08 END
		}
		$subhitscount[$chr-1] ++;
		$subcount++;
	}
	close(INREADFILE);
	print "subcount = $subcount.\n";
}

sub decide()
{
	local($mybaseline, $myregionwid, $mydif, $myprobcutoff)= @_;
	my $lambda = $mybaseline * $myregionwid;
	my @proba;
	my $j;
	my $k;
	my $m;
	my $middle = 30;
	my $upper = 150; # need $upper - $middle > 100.
	my $logsum;
	$proba[0] = exp(-$lambda);
	for($j = 1;$j<$middle;$j++) 	               
	{
		$logsum =0;
		for($m=1;$m<=$j;$m++)                
		{                
			$logsum = $logsum + log($lambda)-log($m);                
		}
		$proba[$j] = exp(-$lambda+$logsum);
	}  
	for($j = $middle;$j<$upper;$j++)
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
	$sum=0;
	my $finish = 0;
	my $result = -1;
	for($k=1;$k<=$mydif;$k++)
	{
		for($m=0;$m<=$middle;$m++)
		{ 
			$sum = $sum + $proba[$m]*$proba[$m+$k];
		}
		if((2*$sum+$trace) > (1-$myprobcutoff))
		{
			$result = 1;
			$finish = 1;
			last;
		} 
	}
	if(($finish == 0)&&((2*$sum+$trace) <= (1-$myprobcutoff)))
	{
		$result = 0;
	}
	$result;
}

sub multifileread
{
	local ($subcurrentfilename, *subtotalcount, *subcount, *subhits, *subhitscount, *subloading) = @_; 
	open(INREADFILE, $subcurrentfilename);
	while($datline = <INREADFILE>)
	{
		my $chr=$posit=$direct=$loading=0;
		chomp($datline);
		my @datas=split /[\s\t]+/, $datline;
		$chr = $datas[0];
		$posit = $datas[1];
		$direct = $datas[2];
		$loading = $datas[3];
		$subtotalcount += $loading;
		if($direct eq "+" || $direct eq "F")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit;
			$subloading[$chr-1][$subhitscount[$chr-1]] = $loading;
			$subhitscount[$chr-1] ++;
			$subcount += $loading;
		}
		elsif($direct eq "-" || $direct eq "R")
		{
			$subhits[$chr-1][$subhitscount[$chr-1]] = $posit + $readlength - $fragmentWidth;
			$subloading[$chr-1][$subhitscount[$chr-1]] = $loading;
			$subhitscount[$chr-1] ++;
			$subcount += $loading;
		}
		else
		{
			print "There is strand error. strand is $direct.\n";
			exit(0);
		}
	}
	close(INREADFILE);
}
