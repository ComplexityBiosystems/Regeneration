#!/bin/bash

#######################
# Here we report exmples of scripts that can be used for reciprocal alignment of transcriptomes using lastz or
# alignment of sequences on Swwissprot database using Blastx.
# ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
# http://www.bx.psu.edu/~rsharris/lastz/
#######################

lastz \
$path'Dataset1.fa'[multiple] $path'Dataset2.fa' \
--noytrim --coverage=20 --identity=30 \
--format=general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2,identity,coverage \
>> $path'Alignment_Dataset1_Dataset2.fa'

blastx -query $path'Dataset1.fa' \
-db $path'/Swissprot/swissprot' \
-outfmt '6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore qframe' \
-evalue 1e-2 -num_threads 16 >> $path'/Alignment_swissprot_Dataset1'


