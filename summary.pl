#!/usr/bin/perl
use POSIX;
use Cwd qw(realpath);
use File::Basename; 
# Steve Qin 11/02/07
# Changed by JC 08/11/2008
# add seperated format 01/09/2009 JC

#$fragmentWidth = 300;
#$readlength = 25;

my $script_path=realpath($0);
my($ef, $dir) = fileparse($script_path);

$argc = @ARGV;
$argc == 7 || die "Provide input files names that include species info, multiple realign.txt filenames, format, read length, fragment size, significant level and result file name.\n";

$sp1 = $ARGV[0];
$innamename = $ARGV[1];
open(INNAMEFILE, $innamename) || die "Can not open $innamename\n";
$format = $ARGV[2];
$format=lc($format);
$readlength = $ARGV[3];
$fragmentWidth = $ARGV[4];
$siglevel = $ARGV[5];
open(OUTFILE1, ">".$ARGV[6].".location.txt") || die "$!\n";
open(OUTFILE2, ">".$ARGV[6].".dis.txt") || die "$!\n";
open(OUTFILE3, ">".$ARGV[6].".reads.txt") || die "$!\n";
open(OUTSUMFILE,">".$ARGV[6].".total.txt") || die "$!\n";
open(OUTPARFILE,">".$ARGV[6].".paras.txt") || die "$!\n";

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

if ($sp1 eq "human" || "hg19") 
{
	$nchr = 24;
	$gcover = 0.9*3100000000;
}
if ($sp1 eq "hg18") 
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

for($i = 0;$i < $nchr; $i++)
{
        $hitscount[$i] = 0;
	$posreadscount[$i] = 0;
	$negreadscount[$i] = 0;
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
	&realignfileread($tmpname,$tmpformat,*totalcount, *count,*hits, *hitscount,*directedreads,*posreadscount, *negreadscount);##### newly added 04/13/08
	print "STAT for $tmpname.\n";
	print "total count = $totalcount.\n";
        print "count = $count.\n";
}
$totalread = $count;

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
        my @pos;
        for($j=0;$j<$hitscount[$i];$j++)
        {
                $pos[$j] = $hits[$i][$j];
        }
        @sortpos = sort {$a <=> $b} @pos;
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
	my @pos;
	for($j=0;$j<$hitscount[$i];$j++)
	{
		$pos[$j] = $hits[$i][$j];
	}
        @sortpos = sort {$a <=> $b} @pos;
	$count = 1;
	$start = $sortpos[0];
	print OUTFILE1 $i + 1,"\t$start\t";
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
#                	        print "Ppang Chr", $i + 1,"\t",$start,"\t", $end, "\t",$regionwid,"\t", $count,"\t", "\n";	
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
#	print "bobe $totalcover, $regionwid\n";
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
print "$totalread\t$allregion\t$hitlevel\n";

@sortpeakwidth = sort {$a <=> $b} @peakwidth;
$medianwidth = $sortpeakwidth[floor(0.5*$totalpeaks)];

print "total number of reads = $totalcount.\n";
print "total number of successfully aligned reads = $totalread. (",100*$totalread/$totalcount," %)\n";
print OUTSUMFILE "Total sequenced reads (millions)\t",$totalcount/1000000,"\n";
print OUTSUMFILE "Total number of reads can be mapped in our database (millions)\t",$totalread/1000000," (",100*$totalread/$totalcount," %)\n";
close(OUTSUMFILE); 

print OUTPARFILE $totalread/1000000,"\n";
print OUTPARFILE $totalinpeakreads/1000000,"\n";
print OUTPARFILE $peakcover/1000000,"\n";
print OUTPARFILE "$medianwidth\n";
print OUTPARFILE "$totalpeaks\n";
print OUTPARFILE $fragmentWidth,"\n";#to accout for linear descreasing at the tail.
print OUTPARFILE $hitlevel,"\n";
close(OUTPARFILE);

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
			$chr =$2 if($datas[0]=~/^(chr|c)?(\d+|X|Y)$/);
			$posit = $datas[1]+1;
			$direct = $datas[-1];
		}
### changed by JC 09/30/2008
		$chr=$1 if($loc=~/chr(\d+|X|Y)(\.fa)?/ || $loc=~/chromosome\.(\d+|X|Y)(\.fa)?/);
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
