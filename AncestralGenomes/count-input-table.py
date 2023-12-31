#!/usr/bin/env python3

import re, sys, os
import pandas as pd

#output table from Roary analysis
pang=pd.read_csv('/home/mcoin/Petrotoga/count/final/gene_presence_absence.csv')

out=open('/home/mcoin/Petrotoga/count/final/input-tab_wagnername.csv', 'w' )
out.write('COGname\tD.tunisiensis\tG.aestuarianus\tG.petraea\tM.sp.1197\tM.sp.38H-ov\tM.lauensis\tM.piezophila\tP.sp.9PWA.NaAc.5.4\tP.halophila\tP.japonica\tP.mexicana\tP.miotherma\tP.mobilis\tP.olearia\tP.sibirica\tT.spiralis\n')
#short species names with which i saved the files
pets=['DEF-TUN', 'GEO-AES', 'GEO-PET', 'PET-HALO', 'PET-JAP', 'PET-MEX', 'PET-MIO', 'PET-MOBI', 'PET-OLEO', 'PET-SIBI', 'TEP-SPI', 'MAR-PIE', 'MAR-LAU', 'MAR-1197', 'PET-9PWANAC54', 'MAR-38HOV']

pdict=dict()

for pang_idx in pang.index.values.tolist() :
    #for each gene family in the pangenome
    name=pang.loc[pang_idx, 'Gene']
    if name not in pdict:
        pdict[name]=dict()

    for p in pets:
        #localize the cell containing the CDS of the specie associated to the gene family
        genes=pang.loc[pang_idx, p]
        try:
            #if there are CDS count how many
            num=genes.count('_')
        except:
            num=0
            
        if p not in pdict[name]:
            pdict[name][p]=0
        pdict[name][p]+=num

#print table in the output
for key1 in sorted(pdict.keys()):
    petsl=''
    print(sorted(pdict[key1].keys()))
    for key2 in sorted(pdict[key1].keys()):
        petsl+='\t{}'.format(pdict[key1][key2])
    out.write('{}{}\n'.format(key1, petsl))

out.close()
