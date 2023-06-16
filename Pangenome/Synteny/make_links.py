#!/usr/bin/env python3

import re, sys, os
import pandas as pd

#pangenome output
pang_sum= pd.read_csv('/home/mcoin/Petrotoga/structure/gene_presence_absence.csv')
pet_list=['PET-OLEO', 'PET-SIBI', 'PET-MEX', 'PET-MIO', 'PET-HALO', 'PET-MOBI', 'PET-JAP']
links=open('/home/mcoin/Petrotoga/structure/links-new.csv', 'a')
links.write('link1\tlink2\n')

for i in range(6):
    #open first gff-format proteome
    file1=open('/home/mcoin/Petrotoga/structure/core-gff/{}.gff'.format(pet_list[i]))
    for line1 in file1:
        g=re.search('##sequence-region (P\.[a-z]+)', line1)
        if g:
            genome1=g.group(1)
        else:
            #find gene name, index in the pangenome table, start/end positions in the genome
            s=re.search(r'locus_tag=([A-Z]{8}_[0-9]{5})', line1)
            if s:
                cds1=s.group(1)
                idxl=pang_sum.loc[pang_sum.isin([cds1]).any(axis=1)].index.tolist()
                idx=idxl[0]
                cdss=pang_sum.loc[idx, pet_list[i+1]]
                pos=re.search(r'CDS\t(\d+)\t(\d+)\t', line1)
                #sometimes it's rna instead of CDS
                if pos:
                    start1=pos.group(1)
                    end1=pos.group(2)
                #open second proteome
                file2=open('/home/mcoin/Petrotoga/structure/core-gff/{}.gff'.format(pet_list[i+1]))
                for line2 in file2:
                    h=re.search('##sequence-region (P\.[a-z]+)', line2)
                    if h:
                        genome2=h.group(1)
                    else:
                        c=re.search(r'locus_tag=([A-Z]{8}_[0-9]{5})', line2)
                        if c:
                            cd=c.group(1)
                            if cd in cdss:
                                po=re.search(r'CDS\s+(\d+)\s+(\d+)\s+', line2)
                                start2=po.group(1)
                                end2=po.group(2)
                                links.write('{}, {}, {}\t{}, {}, {}\n'.format(genome1, start1, end1, genome2, start2, end2))
                file2.close()
    file1.close()
links.close()
