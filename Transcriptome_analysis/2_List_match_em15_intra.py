
import collections as col
#import cProfile as cP
# -*- coding: utf-8 -*-

##################################
## This program takes the list of hits on Swissprot
## and gives the list of matches between annotated genes G as:
## G_1          G_2  #of_matches #tot_len_list #%_su_len_list #tot_len_list     #%_su_len_list        Name
##   0	       11991	    1	        1	        1.0	            9	          0.1111111 	gi|338..|sp|Q9U..|..._HUMAN
##   0	        245	        1	        1	        1.0	            57	          0.0175438 	gi|338..|sp|Q9U..|..._HUMAN

##################################
import sys
print sys.argv

path=sys.argv[1]
name1=sys.argv[2]
filename1=sys.argv[3]


cluster=col.defaultdict(list)
#list_matches=col.defaultdict(list)
f_in_best_matches=open(path+name1+filename1, 'r')
f_out=open('./Clustered_Best_matches_list_euk_em15.txt', 'w')

while(True):
    l=f_in_best_matches.readline()
    if not l: break
    s=l.strip().split()
    #print s
    for i in range(1,len(s)):
        cluster[str(s[0][8:])].append(str(s[i]))

count =0
for k in sorted(cluster):
    for j in sorted(cluster):
        if(j>k):
            list_matches=list()
            for m in range(0,len(cluster[k])):
                for n in range(0,len(cluster[j])):
                    if(cluster[k][m]==cluster[j][n]):
                        count+=1
                        list_matches.append(cluster[k][m])
            if(count>0):
                f_out.write(str(k)+'\t'+str(j)+'\t'+str(count) +'\t'+ str(len(cluster[k]))+'\t'+ str(float(count)/len(cluster[k]))+'\t'+str(len(cluster[j]))+'\t'+str(float(count)/len(cluster[j])))
                for h in list_matches:
                    f_out.write('\t'+str(h))
                f_out.write('\n')
            count =0
            del list_matches
 
