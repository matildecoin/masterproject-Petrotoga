#!/usr/bin/env python3

import re, sys, os
import pandas as pd
#a file with all the proteomes concatenated (Prokka output)
prot_file=open('/home/mcoin/Petrotoga/count/final/Petrotogaceall-faa-conc.faa')
#Roary Petrotogaceae pangenome output
pang_sum= pd.read_csv('/home/mcoin/Petrotoga/count/final/Petrotogaceall/gene_presence_absence.csv')
output=open('/home/mcoin/Petrotoga/count/final/petrotogaceall.faa', 'w')
index_list=list()
write=True
#number of gene families reported to check later with the pangenome output if it corresponds
g=0
for line in prot_file:
    #it's the name of the gene <organism code>_<gene number>
    if line.startswith('>'):
        s=re.search(r'>([A-Z]+_[0-9]+) (\S+)', line)
        cds=s.group(1)
        idxl=pang_sum.loc[pang_sum.apply(lambda col: col.str.contains(cds, na=False), axis=1).any(axis=1)].index.tolist()
        if len(idxl)==0:
            write=False
            continue
        idx=idxl[0]
        #if this gene belongs to a family that has not been reported yet, it will be copied in the output
        if idx not in index_list:
            index_list.append(idx)
            write=True
            line='>{}\n'.format(pang_sum.loc[idx]['Gene'])
            g+=1
        else:
            write=False
    if write:
        output.write(line)

output.close()
prot_file.close()
