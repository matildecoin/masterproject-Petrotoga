#!/usr/bin/env python3

import re
import pandas as pd
import xml.etree.ElementTree as ET

class Hit():
    def __init__(self):
        self.id='id'
        self.num=0
        self.cover=0
        self.identity=0
        self.evalue=0
        self.name=''

class Query():
    def __init__(self, name):
        self.name=name
        self.hits=[]
        self.len=0

    def add_hit(self, hit):
        self.hits.append(hit)

sp_list=['P.halophila', 'P.japonica', 'P.mexicana', 'P.miotherma', 'P.mobilis', 'P.olearia', 'P.sibirica', 'P.sp.ShatinDStank', 'P.sp.HKApet45', 'P.sp.HWPT5563', 'P.sp.HAHPT5561', 'P.sp.9T1HF07', 'P.sp.8T1HF07', 'P.sp.sl27', 'P.sp.9PA5551', 'P.sp.9PWANaAc54', 'D.tunisiensis', 'G.aestuarianus', 'G.petraea', 'T.spiralis', 'M.lausiensis', 'M.piezophila']
#matrix for the number of shared genes
pa_dict=dict()
#dictionary of the gene number
gene_n=dict()
for sp1 in sp_list:
    gene_n[sp1]=0
    pa_dict[sp1]=dict()
    for sp2 in sp_list:
        pa_dict[sp1][sp2]=0
#percent of conserved protein results
pocp=pa_dict.copy()
#all against all BLASTP
file=open('/home/mcoin/Petrotoga/pocp/pcop1.xml')

tree = ET.parse(file)
root = tree.getroot()

for iteration in root.iter('Iteration'):
    for gene in iteration:
        if gene.tag == 'Iteration_query-def':
            cgene = Query(gene.text)
            #each gene in each proteome-so every query-is called <specie name>_<number of the gene>
            sp1=re.search('([A-Z]\..+)_\d{5}', cgene.name).group(1)
            #increase the count of  the genes of the specie1
            gene_n[sp1]+=1
        elif gene.tag == 'Iteration_query-len':
            cgene.len=int(gene.text)
        elif  gene.tag=='Iteration_hits':
            for ithit in gene.iter('Hit'):
                #list all the hits with their parameters
                chit=Hit()
                for hit in ithit:
                    if hit.tag == 'Hit_num':
                        chit.num=int(hit.text)
                    elif hit.tag == 'Hit_id':
                        chit.id = hit.text
                    for ithsp in hit.iter('Hsp'):
                        for hsp in ithsp:
                            if hsp.tag == 'Hsp_num':
                                Hsp_num=int(hsp.text)
                            elif hsp.tag == 'Hsp_evalue' and Hsp_num==1:
                                chit.evalue = float(hsp.text)
                            elif hsp.tag == 'Hsp_query-from'and Hsp_num==1:
                                qstart = int(hsp.text)
                            elif hsp.tag == 'Hsp_query-to'and Hsp_num==1:
                                qend = int(hsp.text)
                                chit.cover = round(100*((qend-qstart)+1)/cgene.len,2)
                            elif hsp.tag == 'Hsp_identity'and Hsp_num==1:
                                nId = int(hsp.text)
                            elif hsp.tag == 'Hsp_align-len'and Hsp_num==1:
                                length = int(hsp.text)
                                chit.identity = round(100 * nId / length, 2)
                                cgene.add_hit(chit)
            #when it reaches the field 'Iteration stat' it went through all the hits                      
        elif gene.tag=='Iteration_stat':
            #make a list not to count just the best hit of every other specie with the query
            list_el=sp_list.copy()
            for hit in cgene.hits:
                if hit.identity>50 and hit.cover>50 and hit.evalue<0.00001:
                    sp2=re.search('([A-Z]\..+)_\d{5}', hit.id).group(1)
                    if sp2 in list_el:
                        pa_dict[sp1][sp2]+=1
                        list_el.remove(sp2)

df=pd.DataFrame.from_dict(pa_dict)
df.to_csv('/home/mcoin/Petrotoga/conserved_proteins1.csv')

for i in range(len(sp_list)):
    sp_1=sp_list[i]
    for j in range(i, len(sp_list)):
        sp_2=sp_list[j]
        if sp1!=sp2:
            pocp[sp_2][sp_1]=(pa_dict[sp_1][sp_2]+pa_dict[sp_2][sp_1])*100/(gene_n[sp_1]+gene_n[sp_2])
            pocp[sp_1][sp_2]=pocp[sp_2][sp_1]
            print('{}\t{}\t{}'.format(sp_1, sp_2, pocp[sp_1][sp_2], pocp[sp_2][sp_1]))
        else:
            pocp[sp_1][sp_2]=100
df=pd.DataFrame.from_dict(pa_dict)
df.to_csv('/home/mcoin/Petrotoga/pocp/pocp1.csv')
