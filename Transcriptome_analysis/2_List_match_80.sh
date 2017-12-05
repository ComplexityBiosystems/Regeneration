#!/bin/bash

path='./'

#### List match annotation


python 1_List_match_em15_intra.py  $path 'Hydra' './Hydra/Best_matches_list_euk_em15.txt' 

awk '{if(($3>10) || (($5 > 0.8) && ($7 > 0.8))){print $1,$2}}' ./Clustered_Best_matches_list_euk_em15.txt > ./List_selected_couples_em15.txt

awk '{c[$1]+=1; c[$2]+=1;}END{for(i in c){print c[i], i}}' ./Clustered_Best_matches_list_euk_em15.txt > ./Frequencies_Selected_couples_em15.txt


python 2_List_match_em15_inter.py $path 'Hydra' './Hydra/Best_matches_list_euk_em15.txt'  'Smed' './Smed/Best_matches_list_euk_em15.txt' 

