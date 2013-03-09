#!/usr/bin/perl 
#use strict;
# Steve Qin 11/30/07
# changed by JC at 08/11/2008

use Cwd qw(realpath cwd);
use File::Basename;
my($fn, $dir) = fileparse(realpath($0));

$argc = @ARGV;
$argc == 4 || die "Provide species, input file name, window size and significance level.\n";

$sp1 = $ARGV[0];
if ($sp1 eq "human") 
{
        $nchr = 24;
        $fold1 = "hg19";
}
elsif ($sp1 eq "mouse") 
{
        $nchr = 21;
        $fold1 = "mm9";
}       
else
{
        print "species error!\n";
        exit(0); 
}       

$name = $ARGV[1];
$windowSize = $ARGV[2];
$siglevel = $ARGV[3];

open(SIZEINFOFILE, $dir."/data/".$fold1."/chromoends.txt");
$count = 0;
my @sizes;
while($dataline = <SIZEINFOFILE>)
{
	chomp($dataline);
	@datas = split(" ",$dataline);
	$sizes[$count] = int($datas[0]/$windowSize) + 1;
	$count ++;
}
for($k=1; $k<($nchr + 1); $k++)
{
	if ($sp1 eq "human")
	{
		$com = "$dir/hmmminus ".$name.".dif.chr".$k.".txt ".$name.".paras.txt ".$name.".hmmout.chr".$k.".txt $windowSize ".$sizes[$k-1];
	}
	elsif ($sp1 eq "mouse")
	{
		$com = "$dir/hmmminus ".$name.".dif.chr".$k.".txt ".$name.".paras.txt ".$name.".hmmout.chr".$k.".txt $windowSize ".$sizes[$k-1];
	}
##	$com = "$dir/hmmminus ".$name.".dif.chr".$k.".txt ".$name.".paras.txt ".$name.".hmmout.chr".$k.".txt $windowSize ".$sizes[$k-1];
	print "$com\n";
	system($com);
}
$com = "perl $dir/minus.region.pl 1 $sp1 $siglevel $windowSize ".$name.".hmmout ".$name.".paras.txt ".$name;
print "$com\n";
system($com);

