#!/usr/bin/env python3

import re, sys, os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET
from functools import reduce
from operator import concat, attrgetter

class Hit():
    def __init__(self):
        self.id='id'
        self.num=0
        self.cover=0
        self.identity=0
        self.evalue=0
        self.org=''


class Query():
    def __init__(self, name):
        self.name=name
        self.hits=[]
        self.len=0

    def add_hit(self, hit):
        self.hits.append(hit)


if len(sys.argv)!=3:
    sys.errot.write('{} <input file> <output name>')
    exit(1)

directory=os.getcwd()
fname=sys.argv[1]
oname='{}/{}.csv'.format(directory, sys.argv[2])
o2name='{}/all_genes2.csv'.format(directory, sys.argv[2])


name_dict={'HALO': 'P. halophila',
'JAPO': 'P. japonica',
'MEXI': 'P. mexicana',
'MIOT': 'P. miotherma',
'MOBI': 'P. mobilis',
'OLEO': 'P. olearia',
'SIBI': 'P. sibirica'}


def parserBlastResults(file, output, output2):
    orga_dict=dict()
    tree = ET.parse(file)
    root = tree.getroot()
    n=0
    for iteration in root.iter('Iteration'):
        for gene in iteration:
            if gene.tag == 'Iteration_query-def':
                cgene = Query(gene.text)
                n+=1
            elif gene.tag == 'Iteration_query-len':
                cgene.len=int(gene.text)
            elif  gene.tag=='Iteration_hits':
                for ithit in gene.iter('Hit'):
                    chit=Hit()
                    for hit in ithit:
                        if hit.tag == 'Hit_num':
                            chit.num=int(hit.text)
                        elif hit.tag == 'Hit_id':
                            chit.id = hit.text
                        elif hit.tag == 'Hit_def':
                            anot = hit.text
                            a=re.findall(r'\[([a-z\.A-Z0-9 _\)\(:\/\=\-]+)\]$', anot)
                            print(a)

                            if len(a)>0:
                                org=a[0]
                            else:
                                print('nope', anot)

                            if len(org)==0:
                                print(anot)
                            if re.search(r'Petrotoga ', anot):
                                print('Petrotoga present:', anot)

                            chit.org=org

                        for ithsp in hit.iter('Hsp'):
                            besthsp=False
                            for hsp in ithsp:
                                if hsp.tag=='Hsp_num':
                                    if int(hsp.text)==1:
                                        besthsp=True
                                    else:
                                        besthsp=False
                                if besthsp and hsp.tag == 'Hsp_evalue':
                                    chit.evalue = float(hsp.text)
                                elif besthsp and hsp.tag == 'Hsp_query-from':
                                    qstart = int(hsp.text)
                                elif besthsp and hsp.tag == 'Hsp_query-to':
                                    qend = int(hsp.text)
                                    chit.cover = round(100*((qend-qstart)+1)/cgene.len,2)
                                elif besthsp and hsp.tag == 'Hsp_identity':
                                    nId = int(hsp.text)
                                    # print(cover)
                                elif hsp.tag == 'Hsp_align-len':
                                    length = int(hsp.text)
                                    chit.identity = round(100 * nId / length, 2)
                                    cgene.add_hit(chit)
            elif gene.tag=='Iteration_stat':
                best_hits=sorted(cgene.hits, key=attrgetter('num'))
                try:
                    best_hit=best_hits[0]
                except:
                    results=[cgene.name, str(len(cgene.hits)), 'no hits', 'no hits', 'no hits', 'no hits']
                    output.write('\t'.join(results) + '\n')
                    output2.write('\t'.join(results) + '\n')
                    continue
                if best_hit.identity<30:
                    output.write('\t'.join([cgene.name, str(len(cgene.hits)),  str(best_hit.org), str(best_hit.cover), str(best_hit.identity), str(best_hit.evalue)]) + '\n')
                output2.write('\t'.join([cgene.name, str(len(cgene.hits)),  str(best_hit.org), str(best_hit.cover), str(best_hit.identity), str(best_hit.evalue)]) + '\n')

                o=re.findall(r'[A-Z][a-z]+ ', str(best_hit.org))
                if o:
                    on=o[-1]
                else:
                    on=best_hit.org
                if best_hit.evalue<0.001 and best_hit.cover>50 and best_hit.identity>50:
                    if on not in orga_dict:
                        orga_dict[on]=0
                    orga_dict[on]+=1



    print('genes: ', n)
    return orga_dict


file=open(fname)
output=open(oname, 'w')
output.write('query\tn.hits\torganism\tcoverage\tidenity\te value\n')
output2=open(o2name, 'w')
output2.write('query\tn.hits\torganism\tcoverage\tidenity\te value\n')
organisms=parserBlastResults(file, output, output2)

k=list(sorted(organisms.keys()))
x_vals=[organisms[key] for key in k]
pal=sns.color_palette("husl", len(x_vals))

f, ax1 = plt.subplots()
ax=sns.barplot(x=x_vals, y=k, palette=[pal[i] for i in range(len(k))], ax=ax1)
ax1.tick_params(axis='both', which='major', labelsize=4)
ax1.spines['right'].set_visible(False)
ax1.bar_label(ax.containers[0], fontsize=4)

plt.savefig('/home/mcoin/Petrotoga/blastp/hist_blastp_{}'.format(sys.argv[2]), dpi=300, bbox_inches='tight')
