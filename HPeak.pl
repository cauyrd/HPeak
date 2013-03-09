#!/usr/bin/env perl -w
#use strict;
#Written by: Steve Qin and Jianjun Yu 
# 05/20/08
#Changed by JC at mouse 08/11/2008
#              at multi 08/30/2008

use Cwd qw(realpath cwd);
use File::Basename;
use Getopt::Long;
Getopt::Long::Configure ('bundling_override');
$|=1;
my($fn, $dir) = fileparse(realpath($0));
my $pwd=cwd();
my $cmd=$0." ".join(" ",@ARGV); ####command line copy
my ($hformat, $format, $sp, $help, $inputfilename, $inputcontrolfilename, $PET, $fmin, $fmax, $readlength, $windowSize, $siglevel, $wigswitch,$seqswitch, $annotation, $name, $interfiles);

#GetOptions('h'=>\$help,'hformat:s'=>\$hformat,'format:s'=>\$format,'sp:s'=>\$sp,'t:s'=>\$inputfilename,'c:s'=>\$inputcontrolfilename,'pet:s'=>\$PET, 'fmin:i'=>\$fmin,'fmax:i'=>\$fmax,'lread|r:i'=>\$readlength,'w|window:i'=>\$windowSize,'s|sig:f'=>\$siglevel,'wig:s'=>\$wigswitch,'seq'=>\$seqswitch,'ann:s'=>\$annotation,'n|name:s'=>\$name, 'interfiles:s'=>\$interfiles);
GetOptions('h'=>\$help,'hformat:s'=>\$hformat,'format:s'=>\$format,'sp:s'=>\$sp,'t:s'=>\$inputfilename,'c:s'=>\$inputcontrolfilename,'pe'=>\$PET, 'fmin|isize:i'=>\$fmin,'fmax:i'=>\$fmax,'lread|r:i'=>\$readlength,'w|window:i'=>\$windowSize,'s|sig:f'=>\$siglevel,'wig'=>\$wigswitch,'seq'=>\$seqswitch,'ann'=>\$annotation,'n|name:s'=>\$name, 'interfiles'=>\$interfiles, 'sec'=>\$SEC, 'nodup'=>\$NODUP);

if(defined $help || (!$hformat && (!$name || !$inputfilename))) {
	print <<DES;
Usage: $fn <-format FORMAT -t TFILE -n NAME> [Options]

Options:
  -h               Show this help message
  -hformat         Help identify the custom format parameter you need, file name file need to be provided
  -format          Format of tag file, "ELAND", "BED" or "Custom". REQUIRED
                   When choosing "Custom", the following info needs to be 
                   specified in the order:
                   1. The column to store genome file where match was found
                   2. The column to store base position of match
                   3. The column to store direction of match
                   The column numbers start at 1 and seperate with comma.
                   Example: the BED format can be presented as -format custom[1,2,6]
                   DEFAULT: BED
  -sp              Speices Option. Currently support "HUMAN" or "MOUSE". DEFAULT: HUMAN 
  -t TFILE         Treatment file name. REQUIRED
                   The name for unique hit file its corresponding multiple hit file must be in one line, and seperated by ';' (unique map file first)
                   and you can also specify the individual format for each line at the end as <format>INDIVIDUAL FORMAT</format>
                   Example:     unique reads file(<format>INDIVIDUAL FORMAT</format>)
                            OR  unique reads file;multi reads file(<format>INDIVIDUAL FORMAT</format>)
  -n, -name        Experiment name used to generate output file. REQUIRED
  -c CFILE         Control file name.
  -pet             Whether ChIP data is pair-ended.
  -pec             Wheter control data is pair-ended.
  -fmin            Minimal DNA fragment size. DEFAULT: 100. If pair-ended, DEFAULT: 200.
  -fmax            Maximal DNA fragment size. DEFAULT: 300. If pair-ended, DEFAULT: 200.
  -r, -lread       Read Length. DEFAULT: 36. if pair-ended, DEFAULT: 100.
  -w, -window      Window size (bp). DEFAULT: 25
  -s, -sig         P-value threshold for peak detection. DEFAULT: 1e-3
  -wig             Whether to generate WIG file for UCSC genome brower
  -seq         	   Whether to extract peak sequences
  -ann         	   Whether to extract nearest gene information for peaks
  -interfiles      Whether to keep the intermediate files
  -nodup           Whether to remove PCR (or other) duplicates

Example: $fn -format BED -sp HUMAN -t example.txt -n Example -fmin 100 -fmax 300 -w 25 -s 1e-3 -wig -seq -ann

> To Idetify the proper custom format you need : $fn -hformat example.txt
         

DES
exit;
}

#### added 10/21/2008
if (defined $hformat)
{
	print "Welcome to the format guide, hope you can find the proper custom format here. \n";
	print "Your Input Name File is $hformat.\n";
	&testformat($hformat);
	exit;
}
else
{
$sp = "HUMAN" if (!$sp);
$sp = lc($sp);
#$PET = lc($PET);
$format="BED" if (!$format);
$format=lc($format);
#if ($PET eq "yes") {
if (defined $PET){
	$PEC = $PET if (!$SEC);
	$fmin = 200 if (!$fmin); 
	$fmax = $fmin;
}
else {
	$fmin = 100 if(!$fmin);
    $fmax = 300 if(!$fmax);
}

if ($fmin > $fmax)
{        
	print "Error: your fmax need to be no less than the fmin value. \n";
        exit;
}
#if ($PET  eq "yes"){
if (defined $PET){
	$readlength = 100 if (!$readlength);
}
else {$readlength = 36 if(!$readlength);}

$windowSize = 25 if(!$windowSize);
$siglevel = 0.001 if(!$siglevel);
$fmin=int($fmin/$windowSize+0.5)*$windowSize;
$fmax=int($fmax/$windowSize+0.5)*$windowSize;
my $fragmentWidth = ($fmin+$fmax)/2;

open(SUMFILE, ">$name.sum");
my $logname = $name.".log";
open(LOGFILE, ">$logname");
print SUMFILE "program started at: ";    
$timestamp = localtime();
print "$timestamp\n"; 
print SUMFILE "$timestamp\n\n";

print SUMFILE "The command used is:\n";
print SUMFILE "$cmd\n\n";
print SUMFILE "Parameters specified:\n";
print SUMFILE "We are dealing with $sp genome data.\n";
print SUMFILE "input treated name: $inputfilename\n";

#######08/30/2008
$t_multi = 0; ###multi tag for treated file
#print SUMFILE "treated sample data is in:\n";
open(INNAMEFILE, $inputfilename);
while($datline = <INNAMEFILE>)
{
	chomp($datline);    
	if(length($datline) ==0)       
	{            
		next;           
	}
	if($datline=~/(.*)\<format\>(.*)\<\/format\>/)
	{
		print SUMFILE "format specified for individual file is detected : $2\n";
		$tmpnames=$1;
		$tmpformat=lc($2);
	}
	else 
	{
		print SUMFILE "Using the same format for all treated data: $format\n";
		$tmpnames=$datline;
		$tmpformat=$inputfilename;
	}
	@names = split(";",$tmpnames);
	$n_name = @names;
	if ($n_name == 2)
	{
		$t_multi = 1;
		print SUMFILE "Multi-mapping file is detected for treatment file, multi-file analysis will be lauched.\n";
		print SUMFILE "unique mapping file: $names[0]\n";
		print SUMFILE "corresponded multi-mapping file: $names[1]\n";
	}
	elsif ($n_name == 1)
	{ 
			print SUMFILE "The treated sample is:\n";
        	print SUMFILE "$datline\n";
        	next;
	}
	else
	{
		die "Error in number of file detected in one line, it should be either one or two files there, please check the manual to correct the input file format.\n";
		exit;
	}   
}
close(INNAMEFILE);
print "finish reading input files names.\n";

if($inputcontrolfilename)
{
	print SUMFILE "input control name: $inputcontrolfilename\n";
#	print SUMFILE "control sample data is in:\n";

	$c_multi = 0 ;###multi tag for control file
	open(INNAMEFILE, $inputcontrolfilename); 
	while($datline = <INNAMEFILE>)
	{
		chomp($datline);          
		if(length($datline) ==0)
		{          
			next;   
		}
		if($datline=~/(.*)\<format\>(.*)\<\/format\>/)
		{
			$tmpnames=$1;
			$tmpformat=lc($2);
			print SUMFILE "format specified for individual file is detected : $2\n";
		}
		else
		{
			$tmpnames=$datline;
			$tmpformat=$format;
			print SUMFILE "Using the same format for all control data: $format\n";
		}
		@c_names = split(";",$tmpnames);
		$nc_name = @c_names;
		if ($nc_name == 2)
		{
			$c_multi = 1;
			print SUMFILE "Multi-mapping file is detected for control files, multi-file analysis will be lauched for control reads.\n";
			print SUMFILE "unique mapping file: $c_names[0]\n";
			print SUMFILE "corresponded multi-mapping file: $c_names[1]\n";
		}
		elsif ($nc_name == 1)
		{
			print SUMFILE "The Control sample is:\n";
			print SUMFILE "$datline\n";
		}
		else
		{
			die "Error in number of file detected in one line, it should be either one or two files there, please check the manual to correct the input control file format.\n";
			exit;
		}
              
	}
	close(INNAMEFILE);
}

print SUMFILE "output file prefix: $name\n";
print SUMFILE "format: $format\n";
#if ($PET eq "yes"){
if (defined $PET){
	print SUMFILE "Data is pair-ended.\n";
	print SUMFILE "Average insert size : $fmin bp.\n";
}
else {
	print SUMFILE "Data is single-ended. \n";
	print SUMFILE "Fragement width : $fmin - $fmax bp.\n";
}
print SUMFILE "Read length is : $readlength bp\n";
print SUMFILE "window size: $windowSize\n";
print SUMFILE "significance level: $siglevel\n\n";
if (defined $NODUP){
	print SUMFILE "Dedupping is requested. \n";
}
my $difcutoff = 1;#####0.01;##### can potentially be specified by the users.

my $tname=$cname=$name;
if($inputcontrolfilename) {
	$tname=$name."_treated";
	$cname=$name."_control";
}

##### first step, process raw data in treated sample and get summary.
print "##### process raw data in treated sample and get summary.\n";
if ($t_multi)
{
	open(INNAMEFILE, $inputfilename);
	open(INUNIQUE, ">".$tname.".unique");
	open(INMULTI, ">".$tname.".multi");
	while ($datline = <INNAMEFILE>)
	{
		chomp($datline);
		if (length($datline) == 0)
		{
			next;
		}
		if($datline=~/(.*)\<format\>(.*)\<\/format\>/)
		{
			$tmpnames=$1;
			$tmpformat=lc($2);
		}
		else
		{
			$tmpnames=$datline;
			$tmpformat=$format;
		}
		@names = split(";",$tmpnames);
		if ($names[0]) 
		{
			print INUNIQUE "$names[0]\<format\>$tmpformat\<\/format\>\n"
		}
		if ($names[1])
		{
			print INMULTI "$names[1]\n"
		}
	}
	close(INMULTI);
	close(INUNIQUE);
	close(INNAMEFILE);
		
	print LOGFILE "perl $dir/summary.multi.pl $sp $tname.unique $tname.multi $format $readlength $fragmentWidth $siglevel $tname \n";
	`perl $dir/summary.multi.pl $sp $tname.unique $tname.multi $format $readlength $fragmentWidth $siglevel $tname`;
	if($?>>8 >0)
	{
		die "Application Error: child process (summary.multi.pl) exited abnormally.\n";
		exit;
	}
}
else
{
	print LOGFILE "perl $dir/summary.pl $sp $inputfilename $format $readlength $fragmentWidth $siglevel $tname \n";
	`perl $dir/summary.pl $sp $inputfilename $format $readlength $fragmentWidth $siglevel $tname`;
	if($?>>8 >0)
	{
		die "Application Error: child process (summary.pl) exited abnormally.\n";
		exit;
	}
}
print SUMFILE "in treated sample:\n";
open(INFOFILE,"$tname.total.txt");
while($datline = <INFOFILE>)
{
	print SUMFILE "$datline";
}
print SUMFILE "\n";
close(INFOFILE);
print "done.\n";
$timestamp = localtime();
print "$timestamp\n";

##### second step, get genome-wide coverage profile.
print "##### get genome-wide coverage profile.\n";
if ($t_multi)
{
	print LOGFILE "perl $dir/chrwisewindow.multi.pl $sp $tname.unique $format $fmin $fmax $readlength $windowSize $tname \n";
	`perl $dir/chrwisewindow.multi.pl $sp $tname.unique $format $fmin $fmax $readlength $windowSize $tname`;
	if($?>>8 >0)
	{
		die "Application Error: child process (chrwisewindow.multi.pl) exited abnormally.\n";
		exit;
	}
}
if (defined $NODUP)
{
	print LOGFILE "java -d64 -Xmx1024m -jar MarkDuplicates.jar I = $inputfilename O = $tname."nodup.sam" AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true";
	`java -d64 -Xmx1024m -jar MarkDuplicates.jar I = $inputfilename O = "$tname.nodup.sam" AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true`;
	print LOGFILE "samtools view -bS "$tname.nodup.sam"  
#if ($PET eq "yes")
if (defined $PET)
{
	print LOGFILE "python $dir/chrwisewindow.py $sp $inputfilename $format $readlength $windowSize $tname \n";
	`python $dir/chrwisewindow.py $sp $inputfilename $format $readlength $windowSize $tname`;
	if($?>>8 >0)
	{
		die "Application Error: child process (chrwisewindow.py) exited abnormally.\n";
		exit;
	}
}
else{
		print LOGFILE "perl $dir/chrwisewindow.pl $sp $inputfilename $format $fmin $fmax $readlength $windowSize $tname \n";
		`perl $dir/chrwisewindow.pl $sp $inputfilename $format $fmin $fmax $readlength $windowSize $tname`;
        if($?>>8 >0) 
	       {   
			 die "Application Error: child process (chrwisewindow.multi.pl) exited abnormally.\n";
             exit;
           }   
}
print "done.\n";
$timestamp = localtime();
print "$timestamp\n";

##### third step, use HMM to get enriched regions.
if(!$inputcontrolfilename)
{
	print "##### use HMM to get enriched regions.\n";
 	print LOGFILE "perl $dir/run.pl $sp $tname $windowSize $siglevel \n";
	`perl $dir/run.pl $sp $tname $windowSize $siglevel`;
	if($?>>8 >0)
	{
		die "Application Error: child process (run.pl) exited abnormally.\n" if($?>>8 >0);
		exit;
	}
##### new 07/13/08
        print LOGFILE "perl $dir/maxsinglebasecover.pl $sp $tname.allregions.txt $tname.reads.txt $fmin $fmax $tname \n";
		`perl $dir/maxsinglebasecover.pl $sp $tname.allregions.txt $tname.reads.txt $fmin $fmax $tname`;
	if($?>>8 >0)
	{
		die "Application Error: child process (maxsinglebasecover.pl) exited abnormally.\n";
		exit;
	}
##### new 07/13/08
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
}
##### fourth step, if there is control, get information.
if($inputcontrolfilename)
{
	print "##### process raw data in control sample and get summary.\n";
	if ($c_multi)
	{
		open(INNAMEFILE, $inputcontrolfilename);
		open(INUNIQUE, ">$cname.unique");
		open(INMULTI, ">$cname.multi");
		while ($datline = <INNAMEFILE>)
		{
			chomp($datline);
			if (length($datline) == 0)
			{
				next;
			}
			if($datline=~/(.*)\<format\>(.*)\<\/format\>/)
			{
				$tmpnames=$1;
				$tmpformat=lc($2);
			}
			else
			{
				$tmpnames=$datline;
				$tmpformat=$format;
			}
			@c_names = split(";",$tmpnames);
			if ($c_names[0]) 
			{
				print INUNIQUE "$c_names[0]\<format\>$tmpformat\<\/format\>\n"
			}
			if ($c_names[1])
			{
				print INMULTI "$c_names[1]\n"
			}
		}
		close(INMULTI);
		close(INUNIQUE);
		close(INNAMEFILE);
		print LOGFILE "perl $dir/summary.multi.pl $sp $cname.unique $cname.multi $format $readlength $fragmentWidth $siglevel $cname \n";
		`perl $dir/summary.multi.pl $sp $cname.unique $cname.multi $format $readlength $fragmentWidth $siglevel $cname`;
		if($?>>8 >0)
		{
			die "Application Error: child process (summary.multi.pl) exited abnormally.\n";
			exit;
		}
	}
	else
	{
		print LOGFILE "perl $dir/summary.pl $sp $inputcontrolfilename $format $readlength $fragmentWidth $siglevel $cname \n";
		`perl $dir/summary.pl $sp $inputcontrolfilename $format $readlength $fragmentWidth $siglevel $cname`;
		if($?>>8 >0)
		{
			die "Application Error: child process (summary.pl) exited abnormally.\n";
			exit;
		}
	}
	print SUMFILE "in control sample:\n";
	open(INFOFILE,"$cname.total.txt");
	while($datline = <INFOFILE>)
    	{
		print SUMFILE "$datline";
	}
	print SUMFILE "\n";
	close(INFOFILE);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print "##### get genome-wide coverage profile.\n";

	if ($c_multi)
	{
		print LOGFILE "perl $dir/chrwisewindow.multi.pl $sp $cname.unique $format $fmin $fmax $readlength $windowSize $cname \n";
		`perl $dir/chrwisewindow.multi.pl $sp $cname.unique $format $fmin $fmax $readlength $windowSize $cname`;
		if($?>>8 >0)
		{
			die "Application Error: child process (chrwisewindow.multi.pl) exited abnormally.\n";
			exit;
		}
	}
	if (defined $PEC)
#	if ($PET eq "yes")
	{
		print LOGFILE "python $dir/chrwisewindow.py $sp $inputcontrolfilename $format $readlength $windowSize $cname \n";
		`python $dir/chrwisewindow.py $sp $inputcontrolfilename $format $readlength $windowSize $cname`;
		if($?>>8 >0)
		{
			die "Application Error: child process (chrwisewindow.py) exited abnormally.\n";
			exit;
		}
	}
	else {
		print LOGFILE "perl $dir/chrwisewindow.pl $sp $inputcontrolfilename $format $fmin $fmax $readlength $windowSize $cname \n";
		`perl $dir/chrwisewindow.pl $sp $inputcontrolfilename $format $fmin $fmax $readlength $windowSize $cname`;
		if($?>>8 >0)
		{
			die "Application Error: child process (chrwisewindow.pl) exited abnormally.\n";
			exit;
		}
	}

	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";

	print "##### get summary of differences between treated and control samples.\n";
	print LOGFILE "perl $dir/summary.minus.pl $sp $t_multi $c_multi $inputfilename $inputcontrolfilename $format $readlength $fragmentWidth $siglevel $name\n"; 
	`perl $dir/summary.minus.pl $sp $t_multi $c_multi $inputfilename $inputcontrolfilename $format $readlength $fragmentWidth $siglevel $name`;
	if($?>>8 >0)
	{
		die "Application Error: child process (summary.minus.pl) exited abnormally.\n";
		exit;
	}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
   	print "##### get genome-wide coverage profile for the differences.\n";

	print LOGFILE "perl $dir/minus.pl $sp $cname $tname $difcutoff $name \n";
	`perl $dir/minus.pl $sp $cname $tname $difcutoff $name`;
	if($?>>8 >0)
	{
		die "Application Error: child process (minus.pl) exited abnormally.\n";
		exit;
	}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print "##### use HMM to get enriched regions.\n";
	print LOGFILE "perl $dir/minus.run.pl $sp $name $windowSize $siglevel \n";
	`perl $dir/minus.run.pl $sp $name $windowSize $siglevel`;
	if($?>>8 >0)
	{
		die "Application Error: child process (minus.run.pl) exited abnormally.\n";
		exit;
	}
##### new 07/13/08              
        print LOGFILE "perl $dir/maxsinglebasecover.minus.pl $sp $name.allregions.txt $tname.reads.txt $cname.reads.txt $fmin $fmax $name \n";
		`perl $dir/maxsinglebasecover.minus.pl $sp $name.allregions.txt $tname.reads.txt $cname.reads.txt $fmin $fmax $name`;
	if($?>>8 >0)
	{        
		die "Application Error: child process (maxsinglebasecover.minus.pl) exited abnormally.\n";
		exit;
	}
##### new 07/13/08
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
}
##### fifth step, find out number of region and total coverage.
open(REGIONFILE, "$name.allregions.txt");
my $count =0;
my $totalsize = $minwidth = $maxwidth = $minheigh = $maxheigh = 0;
while($datline = <REGIONFILE>)      
{
    chomp($datline);
    my @datas = split(" ", $datline);
    my $start = $datas[1];
    my $end = $datas[2];
    my $width = $datas[3];
    my $heigh = $datas[4];
    if(($count==0)||($width > $maxwidth))
    {
        $maxwidth = $width;    
    }
    if(($count==0)||($width < $minwidth))
    {
        $minwidth = $width;
    }
    if(($count==0)||($heigh > $maxheigh))
    {
        $maxheigh = $heigh;
    }
    if(($count==0)||($heigh < $minheigh))
    {
        $minheigh = $heigh;
    }
    $count ++;
    $totalsize = $totalsize + $end-$start+1;
}

close(REGIONFILE);
$totalsize = $totalsize/1000000;
print SUMFILE "Number of enriched regions is $count.\n";
print SUMFILE "Total length of genome covered by the enriched regions: $totalsize Mb.\n";
print SUMFILE "size of the enriched regions range from $minwidth to $maxwidth.\n";
print SUMFILE "coverage for the enriched regions ranges from $minheigh to $maxheigh.\n";  
print SUMFILE "\nResult files:\n";
print SUMFILE "enriched regions are stored in $name.hpeak.out\n";

if(defined $wigswitch)
#if (defined $wigswitch)
{
	if(!$inputcontrolfilename) 
	{
        	print "##### generate wig file for enriched regions.\n";
       		print LOGFILE "perl $dir/wig.region.pl $sp $name.allregions.txt $name $name $windowSize \n";
			`perl $dir/wig.region.pl $sp $name.allregions.txt $name $name $windowSize`;
		if($?>>8 >0)
		{
			die "Application Error: child process (wig.region.pl) exited abnormally.\n";
			exit;
		}
	} 
	else 
	{
        	print "##### generate wig file for enriched regions between treated and control samples.\n";
        	print LOGFILE "perl $dir/wig.minus.region.pl $sp $name.allregions.txt $tname $cname $name $windowSize \n";
			`perl $dir/wig.minus.region.pl $sp $name.allregions.txt $tname $cname $name $windowSize`;
		if($?>>8 >0)
		{
			die "Application Error: child process (wig.minus.region.pl) exited abnormally.\n";
			exit;
		}
	}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "WIG file of enriched regions is in $name.wig\n";
}

if(defined $seqswitch)
{
	print "##### obtain DNA sequences for enriched regions.\n";		
	print LOGFILE "perl $dir/get_seq.pl $sp $name.allregions.txt \n";
	`perl $dir/get_seq.pl $sp $name.allregions.txt`;
	if($?>>8 >0)
	{
		die "Application Error: child process (get_seq.pl) exited abnormally.\n";
		exit;
	}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "FASTA format sequences in enriched regions are stored in $name.seq.fa\n";         
}

if (defined $annotation)
#if($annotation == 'Yes')
{
	print "##### produce detailed information for enriched regions.\n";
	if ($sp eq "human" || "hg19")
	{
		print LOGFILE "perl $dir/peakann_hg19.pl -a -o $name.annotation.txt $name.allregions.txt \n";
		`perl $dir/peakann_hg19.pl -a -o $name.annotation.txt $name.allregions.txt`;
		if($?>>8 >0)
		{
			die "Application Error: child process (peakann.pl) exited abnormally.\n";
			exit;
		}
	}
	elsif ($sp eq "hg18")
	{
		print LOGFILE "perl $dir/peakann_hg18.pl -a -o $name.annotation.txt $name.allregions.txt \n";
		`perl $dir/peakann_hg18.pl -a -o $name.annotation.txt $name.allregions.txt`;
		if($?>>8 >0)
		{
			die "Application Error: child process (peakann.pl) exited abnormally.\n";
			exit;
		}
	}
	elsif ($sp eq "human" || "mm9")
	{
        print LOGFILE "perl $dir/peakann_mm9.pl -a -o $name.annotation.txt $name.allregions.txt\n";
        `perl $dir/peakann_mm9.pl -a -o $name.annotation.txt $name.allregions.txt`;
        if($?>>8 >0)
		{
			die "Application Error: child process (peakann.pl) exited abnormally.\n";
			exit;
		}
    }
	else {print "Error: species info error when generating annotation files.\n";}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "Detailed annotation of enriched regions is in $name.annotation.txt\n";         
}
close(LOGFILE);
#exit(0);

#### sixth step, clean up intermediate files, output wig, seq and annotation file if required.
if (defined $interfiles)
#if (defined $interfiles)
{
	print "choose to keep all the intermediate files.\n";
}
else
{
	print "##### clean up intermediate files.\n";
	unlink <$tname.windowhitscount.chr*.txt>;
	unlink <$tname.select.chr*.txt>;
	unlink <$tname.location.txt>;
	unlink <$tname.paras.txt>;
	unlink <$tname.total.txt>;
	unlink <$tname.reads.txt>;
	unlink <$tname.dis.txt>;  
	unlink <$cname.windowhitscount.chr*.txt> if($inputcontrolfilename);
	unlink <$cname.select.chr*.txt> if($inputcontrolfilename);
	unlink <$cname.location.txt> if($inputcontrolfilename);
	unlink <$cname.paras.txt> if($inputcontrolfilename);
	unlink <$cname.total.txt> if($inputcontrolfilename);   
	unlink <$cname.reads.txt> if($inputcontrolfilename);
	unlink <$cname.dis.txt> if($inputcontrolfilename);
	unlink <$name.paras.txt> if($inputcontrolfilename);
	unlink <$name.hmmout.chr*.txt>; 
	unlink <$name.dif.chr*.txt>;
	unlink <$name.ratio.minus.txt>; 
	#unlink <poisson.out>;
	unlink <$name.basecover.out>;
}
print "done.\n";
$timestamp = localtime();
print "$timestamp\n";
print SUMFILE "\nprogram ended at: ";
$timestamp = localtime(); 
print "$timestamp\n";
print SUMFILE "$timestamp\n";
close(SUMFILE);
}

#### added 10/21/2008 JC
sub testformat
{
	local ($subhformat) = @_;
#	print $subhformat;
	open(INTMP, $subhformat) || die "$!\n";
	while ($datline = <INTMP>)
	{
		chomp($datline);
		if (length($datline) != 0)
		{
			@names = split(";",$datline);
			$tmpname = $names[0];
			last;
		}
		print "error!! there is no file name in your input file.\n";	
	}
	close(INTMP);
	print "the file we take here is $tmpname\n";
	print "here we just show a typical dataline in your file and the corresponding number of each item.\n";
	print "what you want should be 1. chromosome number (like chr1 or chr1.fa ), 2. start site, 3. starnd direction (+ - F R)\n";
	if($tmpname =~ /gz$/)
	{
		system("gunzip ".$tmpname);
		$namelength = length($tmpname);
		$newname = substr($tmpname,0,$namelength-3);
		open(INREADFILE, $newname) || die "$!\n";
	}
	elsif($tmpname =~ /zip$/)
	{
		system("unzip ".$tmpname);
		$newname = substr($tmpname,0,$namelength-4);
		open(INREADFILE, $newname.".txt") || die "$!\n";
	}
	else
	{
		open(INREADFILE, $tmpname) || die "$!\n";
	}
	while ($datline = <INREADFILE>)
	{
		chomp($datline);
		if (index($datline,'.fa')==-1){
		    next if not (($datline =~ /chr(\d+|X|Y)?/)||($datline =~ /chromosome\.(\d+|X|Y)?/));
		}else{
		next if not (($datline =~ /chr(\d+|X|Y)(\.fa)?/)||($datline =~ /chromosome\.(\d+|X|Y)(\.fa)?/));
	    }
		my @datas=split /[\s\t]+/, $datline;
		$tmpn=1;
		print "bobe: $datline\n";
		foreach $a (@datas)
		{
			print "        the $tmpn th item: $a\n";
			$tmpn += 1;
		}
		last;
	}
	close(INREADFILE);
	print "please choose custom[a,b,c] according to the list content of your file. thanks again to choose HPeak :)\n";
}

