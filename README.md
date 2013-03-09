Description
===========

HPeak is a hidden Markov model-based approach that can accurately pinpoint regions to where significantly more sequence reads map. Testing on real data shows that these regions are indeed highly enriched by the right protein binding sites.

Command (single-end):

    perl /compbio/software/HPeak3/HPeak.pl -sp HUMAN/MOUSE -format BED -t TREATMENT.inp -c CONTROL.inp -n OUTPUTPREFIX -fmin 100 -fmax 300 -r 36 -ann -wig -seq -interfiles
Command (pair-end):

    perl /compbio/software/HPeak3/HPeak.pl -sp HUMAN/MOUSE -format BED -pe TREATMENT.inp -c CONTROL.inp -n OUTPUTPREFIX -isize 200 -r 100 -pec (if control is PE) -ann -wig -seq –interfiles
note:

1.  Default species is HUMAN. Also supports MOUSE. Can add any other genome if in need.
2.   Default format is BED. Also support ELAND. Will add SAM and BAM.
3.   –r specifies read length (this is import through my experience).
4.   –pe indicates pair-ended data.
5.   If data is pair-ended, -isize refers to insert size (total length of a pair).


Features
========
* Accurate peak calling with detailed information including summit location and p-value
* Thorough peak annotations including enriched genes symbol, refseq ID, summit distance and even conversation score
* Fully supports pair-end data by treating two ends together as ONE fragment
* Prepares necessary outputs ready for down-steam analysis such as enrichment curve and heatmap

