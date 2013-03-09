#!/usr/bin/perl 
use POSIX;
use Cwd qw(realpath cwd);
use File::Basename;
my($fn, $dir) = fileparse(realpath($0));
# Steve Qin 05/15/08
# Changed by JC at 08/11/2008

$argc = @ARGV;
$argc == 8 || die "Provide input files names that include species info, multiple realign.txt filenames, file format, minimum and maximum DNA fragment length, read length, window size, and result filename.\n";

$threshold = 1;
#$readlength = 25;

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

$innamename = $ARGV[1];
open(INNAMEFILE, $innamename); 
$format = $ARGV[2];
$fragmin = $ARGV[3];
$fragmax = $ARGV[4];
$fragmentWidth = $fragmax; #####($fragmin + $fragmax)/2;
$readlength = $ARGV[5];
$resolution = $ARGV[6];
$frontname = ">".$ARGV[7].".windowhitscount.";
$selectname = ">".$ARGV[7].".select.";
open(CHROMOFILE,"$dir/data/".$fold1."/chromoends.txt");

$count = 0;
my @filenames;
while($datline = <INNAMEFILE>)
{
        chomp($datline);
        if(length($datline) ==0)
        {
                next;
        }
	$filenames[$count] = $datline;
	$count ++;
}
$filecount = $count;
print "There are $filecount files.\n";
for($i = 0;$i < $nchr; $i++)
{
        $plushitscount[$i] = 0;
	$minushitscount[$i] = 0;
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
	&realignfileread($tmpname,$tmpformat,*totalcount, *count,*plushits, *minushits, *plushitscount, *minushitscount);
	print "total count = $totalcount.\n";
        print "$i\tcount = $count.\n";
}
$totalread = $count;

$count = 0;
my @chromosize;  
while($datline = <CHROMOFILE>)
{
        chomp($datline);
        @datas = split(" ",$datline);
	$chromosize[$count] = $datas[0];
	print "$chromosize[$count]\n";
	$count++;
	if($count == $nchr)
	{
		last;
	}
}
close(CHROMOFILE);

my @added;
$firstpart = $fragmin / $resolution;
$secondpart = ($fragmax - $fragmin) / $resolution;
if ($secondpart == 0)
{
	$totalparts = $fragmax / $resolution;
	$dimension = $totalparts + 1;
	for($j=0;$j<$resolution;$j++)
	{
		$startfrac = 1-$j/$resolution;
		$endfrac = 1-$startfrac;
		$added[$j * $dimension] = $startfrac ;
#		$added[$j * $dimension] = 1;
		for($k = 1;$k < $totalparts;$k++)
		{
			$added[$j * $dimension + $k] = 1;
		}
		$added[$j * $dimension + $totalparts] = $endfrac;
	}
#print "bobe\n@added\n";
#print "dimension $dimension, resolution $resolution $ARGV[5]\n";
#exit(0);
}
else
{
	$inv = 1/$secondpart;
	$totalparts = $fragmax / $resolution;
	$dimension = $totalparts + 1;
	for($j=0;$j<$resolution;$j++) 
	{
		$startfrac = 1-$j/$resolution;
        	$endfrac = 1-$startfrac;         
		$added[$j * $dimension] = $startfrac ;
#		$added[$j * $dimension] = 1;
        	for($k = 1;$k < $firstpart;$k++)
        	{
                	$added[$j * $dimension + $k] = 1;
        	}
		$added[$j * $dimension + $firstpart] = $endfrac + $startfrac * (1 - 0.5 * $startfrac * $inv);        
		for($k = $firstpart + 1;$k < $totalparts;$k++)
        	{
                	$added[$j * $dimension + $k] = 1 - ($k - $firstpart - 0.5 + $startfrac) * $inv;
        	}
		$added[$j * $dimension + $totalparts] = 0.5 * $endfrac * $inv *$endfrac;
	}
#print "bobe\n@added\n";
#exit(0);
}

for($i = 0;$i < $nchr; $i++)
{
	$chromosome = $i +1;
	open(OUTFILE, $frontname."chr".$chromosome.".txt");	
        open(OUTSELECTFILE, $selectname."chr".$chromosome.".txt");
	my @pluspos;
	for($j=0;$j<$plushitscount[$i];$j++)
	{
		$pluspos[$j] = $plushits[$i][$j];
	}
        @sortpos = sort {$a <=> $b} @pluspos;
	$plustotal = $plushitscount[$i];

	$upperlimit = floor($chromosize[$i]/$resolution) + 1;
	for($j = 0;$j < $upperlimit; $j++)
        {
		$heights[$j] = 0;		
	}
	print "chr ", $i+1, " upperlimit = $upperlimit.\n";
	for($j = 0;$j < $plustotal; $j++)
        {
		$lowend = floor(($sortpos[$j]-1) /$resolution);##### 03/26/08 ( -1)
                $residual = $sortpos[$j] - 1 - $lowend *$resolution;
        	if(($lowend >= 0)&&($lowend < ($upperlimit - $totalparts)))
		{
			for($k=0;$k<($totalparts+1);$k++)
        	        {
				$heights[$lowend + $k] = $heights[$lowend + $k] + $added[$residual * $dimension + $k];
#				print "heights ($lowend + $k), ($residual * $dimension + $k),$added[$residual * $dimension + $k]\n";
			}
		}
#		exit(0);
	}
#	print "pos: $sortpos[1]\n";
#	exit(0);
	my @minuspos;
        for($j=0;$j<$minushitscount[$i];$j++)
        {
                $minuspos[$j] = $minushits[$i][$j];
        }
        @sortpos = sort {$a <=> $b} @minuspos;
        $minustotal = $minushitscount[$i];
	for($j = 0;$j < $minustotal; $j++)
        {
                $lowend = floor(($sortpos[$j]-1) /$resolution);##### 03/26/08 ( -1)
                $residual = $sortpos[$j] - 1 - $lowend *$resolution;
		$residual = $resolution -1 - $residual;
                if(($lowend >= 0)&&($lowend < ($upperlimit - $totalparts)))
                {
                        for($k=$totalparts;$k>=0;$k--)
                        {
                                $heights[$lowend + $k] = $heights[$lowend + $k] + $added[$residual * $dimension + ($totalparts - $k)];
			}
                }
        }
#	print "minus: $sortpos[1]\n";
#	exit(0);
	$selected = 0;
        for($j=0;$j<$upperlimit;$j++) 
        { 
		$start = $resolution*$j + 1;
		$end = $start + $resolution -1;
		if($i == ($nchr - 2))
		{
			print OUTFILE "chrX\t",$start,"\t", $end, "\t",$heights[$j],"\n";
		}
		elsif($i == ($nchr-1))
		{
			print OUTFILE "chrY\t",$start,"\t", $end, "\t",$heights[$j],"\n";
		}
		else
		{
			print OUTFILE "chr", $i + 1,"\t",$start,"\t", $end, "\t",$heights[$j],"\n";	
		}
                if($heights[$j]>$threshold)
                {
                        print OUTSELECTFILE "$j\t$heights[$j]\n";
                	$selected ++;
		}
	} 
	print "chr ", $i+1," # nonempty windows = ",$selected,"\n";
	close(OUTFILE);
	close(OUTSELECTFILE);
}

print "total number of reads = $totalcount.\n";
print "total number of successfully aligned reads = $totalread. (",100*$totalread/$totalcount," %)\n";

sub realignfileread()
{
	local ($subcurrentfilename, $subformat, *subtotalcount, *subcount, *subplushits, *subminushits, *subplushitscount, *subminushitscount) = @_; 
	if($subcurrentfilename =~ /gz$/)
	{
		system("gunzip ".$subcurrentfilename);
		$namelength = length($subcurrentfilename);		
		$newname = substr($subcurrentfilename,0,$namelength-3); 
		print "name = $newname\n";
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
				$posit=$datas[$infos[1]-1];
				$direct=$datas[$infos[2]-1];
			}
		}
		elsif ($subformat eq "bed") 
		{
			$chr =$2 if($datas[0]=~/^(chr|c)?(\d+|X|Y)$/);
			$posit = $datas[1]+1;
			$direct = $datas[-1];
		}
		$chr=$1 if($loc=~/chr(\d+|X|Y)\.fa/ || $loc=~/chromosome\.(\d+|X|Y)\.fa/);
		next if(!$chr);
		$chr=($nchr-1) if($chr eq "X");
		$chr=$nchr if($chr eq "Y");
		if($direct eq "+" || $direct eq "F")
		{
			$subplushits[$chr-1][$subplushitscount[$chr-1]] = $posit;
#			print "$posit \n";
#			exit(0);
			$subplushitscount[$chr-1] ++;
		}
		elsif($direct eq "-" || $direct eq "R")
		{
			$subminushits[$chr-1][$subminushitscount[$chr-1]] = $posit + $readlength - $fragmentWidth;
#			print "$posit + $readlength - 1\n";
#			exit(0);
			$subminushitscount[$chr-1] ++;
		}
		else
		{
			print "There is strand error. strand is $direct, position is $posit...............\n\n\n";
			exit(0);
		}
		$subcount++;
	}
close(INREADFILE);
}

