#!/bin/bash

#######################
#This program is an example of the pipeline used to annotate to Swissprot a transcriptome. 
#input: the result of blastx alignment of transcripts on Swissprot database as in 0_Transcriptome_alignment.sh
#output: a list of selected Swispprot hits on genes
#######################

rm Intermediate* 

###Uses CD-Hit to cluster transcripts into genes
### Note that CD-Hit produces more than an output file see https://github.com/weizhongli/cdhit

$path_CDHIT/cdhit-master/cd-hit-est \
 -i <fasta_infile> -o <outfile> -c 0.80 -T 8 -M 0

awk '{if($1==">Cluster"){cluster=$2; }else{d[substr($3,2,length($3)-4)]=cluster;}}END{for(i in d){print i, "Cluster_"d[i]}}' <outfile>.clstr > <file_genes_and_transcripts>


###Save the size of "cluster of transcripts" from the output file of CD-Hit <file_genes_and_transcripts> ######

awk '{if($1 == ">Cluster"){if(count>0){print clust_name, count}; clust_name=substr($1,2)"_"$2;  count=0;}\
     else{count +=1;}}END{print clust_name, count}' <file_genes_and_transcripts> \
     > <outfile_CountsClusterSize>

### Sort the start/stop coordinates and substitute tanscripts with genes
###  <file_genes_and_transcripts> contains association genes/transcripts 
###  <alignment_swissprot_database> contains the reuslt of blastx alignment

awk \
'{if(NR==FNR){c[$2]=$1}\
  else{\
     if($1 in c){if($7 < $8){start=$7; stop=$8}else{start=$8; stop=$7};\
      print c[$1],$1,$3,start,stop,$9,$10 > "./Intermediate.txt"; next;};\
      }}' k=0\
      <file_genes_and_transcripts> <alignment_swissprot_database>


### Select eukaryotes entries
###  <speclist> contains Uniprot association between species abbreviation and name as downloaded from Uniprot database
###  e.g. PIG   E    9823: N=Sus scrofa  
###  "E" means eukaryotes, and "S=Hepatitis E virus genotype" is the only exception


awk \
 '{if(NR==FNR){if($2=="E"){if(length($1)==5){c["_"$1]=substr($3, 1,length($3-1));}\
                           else{if($1!="S=Hepatitis"){d["_"$1]=substr($3, 1,length($3-1))}}}}\
   else{if(substr($3, length($3)-5, length($3)) in c){print $0;}\
        else{\
            if(substr($3, length($3)-3, length($3))in d){print $0;}\
            else{\
                if(substr($3, length($3)-4, length($3))in d){print $0;}}}}}'\
   <speclist> ./Intermediate.txt \
    > Intermediate_eukaryotes.txt

### Averages the alignment over all the tanscripts and select best matches

sort -k1,1 -k3,3 -k2,2 -k6,6g  ./Intermediate_eukaryotes.txt > ./Intermediate_eukaryotes_sorted.txt

awk '{if(NR==FNR){c[$1]=$2;next;\
      }else{if(isogr==$1){\
                    if(gene==$3){\
                            if(isotig != $2){\
                                count++; start+=$4; stop+=$5; pval+=$6; score+=$7; isotig=$2;};\
                     }else{\
                            print isogr, gene, int(start/count), int(stop/count), pval/count, score/count,(100*count/c[isogr]);\
                            gene=$3; isotig=$2; start=$4; stop=$5; pval=$6; score=$7; count=1}\
         }else{\
               if(count>0){print isogr,gene, int(start/count), int(stop/count),pval/count, score/count,(100*count/c[isogr]);}\
               isogr=$1; gene=$3; isotig=$2; start=$4; stop=$5; pval=$6; score=$7; count=1;}\
     };}' isogr="isogr0" gene="gene0" isotig="tig0" start=0 stop=0 pval=0 score=0 count=0 \
     <outfile_CountsClusterSize> ./Intermediate_eukaryotes_sorted.txt\
    | sort -k1,1 -k7,7rg -k6,6rg > ./Alignment_swissprot_eukaryotes.txt


rm Intermediate*

awk '{if(isogr!=$1){print $0 > "Best_matches_euk.txt"; isogr=$1; count=0}\
 else{if(count==0){print $0 > "Second_best_matches_euk.txt"; count=1;}}}'\
     isogr="-1" count=0 ./Alignment_swissprot_eukaryotes.txt

awk '{if($5<=1e-15){print $0 > "Best_matches_euk_em15txt"}}' Best_matches_euk.txt

awk '{if($5<=1e-15){d[$1][$2]=$1}}END{for(m in d){printf m; printf "\t";for(k in d[m]){printf k; printf "\t";}printf "\n";}}'  \
    ../Alignment_swissprot_database_eukaryotes_cluster_averaged.txt > Best_matches_list_euk_em15.txt

awk '{if((isogr!=$1)&&(substr($2, length($2)-4, length($2))=="HUMAN")){print $0 > "Best_matches_euk_Human.txt"; isogr=$1; count=0}}'\
    isogr="-1" count=0 ./Alignment_swissprot_eukaryotes.txt

awk '{if((isogr!=$1)&&(substr($2, length($2)-4, length($2))=="HUMAN")&&($5 <=1e-15)){print $0 > "Best_matches_euk_Human_1em15.txt"; isogr=$1; count=0}}'\
    isogr="-1" count=0 ./Alignment_swissprot_eukaryotes.txt






