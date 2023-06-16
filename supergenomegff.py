#!/usr/bin/env python3

import re
import pandas as pd

#proteome of one of the species of the pangenome
ref_gff=open('Petrotoga/gene_annotation/gff/GEO-AES.gff')
#roary pangenome output table
pang_sum= pd.read_csv('Petrotoga/pangenome/roary/GOT/gene_presence_absence.csv')

output_gff=open('/Petrotoga/pangenome/geopang.gff', 'w')

#if write is TRUE the line will be copied in the output
for line in ref_gff:
    write=True
    s=re.search(r'ID=([A-Z]+_[0-9]+);', line)
    if line.startswith('N'):
        #all the gene lines start with N in this .gff file, it's the first letter of the names of the contigs
        write=False
    if s:
        write=True
        cds=s.group(1)
        #find the gene in the pangenome table
        idxl=pang_sum[pang_sum['GEO-AES']==cds].index.values
        if len(idxl)==0:
            print('missing gene {}'.format(cds))
            continue
        idx=idxl[0]
        #just copy the line if it's a core gene 
        if pang_sum.at[idx, 'No. sequences']!=3:
            write=False
    if write:
        output_gff.write(line)
output_gff.close()
ref_gff.close()
