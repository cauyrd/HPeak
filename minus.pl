#!/usr/bin/perl
use POSIX;
# Steve Qin 12/03/07
# changed by JC at 08/11/2008

$argc = @ARGV;
$argc == 5 || die "Provide species, unsimulated, simulated data files, threshold and result file name.\n";

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

$unsimfilename = $ARGV[1];
$simulfilename = $ARGV[2];
$threshold = $ARGV[3];
$outputfilename = $ARGV[4];

open(PARASFILE,$outputfilename.".paras.txt");
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

for($j=1;$j<=$nchr;$j++)
{
	print "Chromosome $j:\n";
	&subtract($j,$unsimfilename,$simulfilename,$outputfilename,$folds,$threshold);
}

sub subtract()
{
	local($mychr, $myunsimfilename,$mysimulfilename,$myoutputfilename, $myfolds, $mythreshold) = @_;
	$localunsimfilename = $myunsimfilename.".windowhitscount.chr".$mychr.".txt";
	$localsimulfilename = $mysimulfilename.".windowhitscount.chr".$mychr.".txt";
	$localoutputfilename = $myoutputfilename.".dif.chr".$mychr.".txt";
	open(UNSIMFILE, $localunsimfilename); 
	open(SIMFILE, $localsimulfilename); 
        open(OUTPUTFILE, ">".$localoutputfilename);
	my $count = 0;
	while($datline = <UNSIMFILE>)
	{
        	chomp($datline);
		@datas = split(" ", $datline);
		$bghits = $datas[3];
		$datline = <SIMFILE>;
		chomp($datline);
        	@datas = split(" ", $datline);
        	$hits = $datas[3];
		$dif = $hits - $bghits * $myfolds;	
		if($dif >= $mythreshold)
		{
			print OUTPUTFILE "$count\t$dif\t$hits\t",$bghits * $myfolds,"\n";
		}
		$count ++;
	}
close(UNSIMFILE);
close(SIMFILE);
close(OUTPUTFILE);
}
