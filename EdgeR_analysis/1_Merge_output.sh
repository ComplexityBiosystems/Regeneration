#!/bin/bash

list_time=("1 2 3 4")
path="./"

for time in $list_time
do                                        
 awk 'NR>1{{print $1,$2,$5}}'  $path'Time0_'$time'.de' | sort -u -k 1,1 > 'buffer_'$time
done

paste buffer_1  buffer_2 buffer_3 buffer_4| \
awk '{if(($1==$4)&&($1==$7)&&($1==$10)){print $1,$2,$5,$8,$11}}' > $path'EdgeR_logFC.dat'

paste buffer_1  buffer_2 buffer_3 buffer_4| \
awk 'function abs(v) {return v < 0 ? -v : v}{for(i=2; i<= NF; i+=3){if((abs($i)>2)&&($(i+1)<0.05)){printf $1; for(j=2; j<=NF; j+=3){printf " ";printf $j}; printf "\n"; next;} } }'> $path'EdgeR_logFC_DE2_pval005.dat'

done
