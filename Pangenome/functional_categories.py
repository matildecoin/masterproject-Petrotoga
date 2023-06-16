#!/usr/bin/env python3

import re, sys, os
import pandas as pd
import xml.etree.ElementTree as ET
from operator import attrgetter

class Hit():
    def __init__(self):
        self.id='id'
        self.COGname=''
        self.COGcat=[]
        self.funct=''
        self.num=0
        self.cover=0
        self.identity=0
        self.evalue=0

    def set_categories(self, categories):
        self.COGcat=categories

class Query():
    def __init__(self, name):
        self.name=name
        self.hits=[]
        self.len=0

    def add_hit(self, hit):
        self.hits.append(hit)

def add_to_dict(categories, n, f):
    if categories!=['NS']:
        f+=1
    for cat in categories:
        n=n+1
        if cat not in ann_dict:
            ann_dict[cat]=0
        ann_dict[cat]+=1
    return n, f

def parserBlastResults(file, g):
    a=0
    f=0
    n=0
    ns=0
    save=True
    tree = ET.parse(file)
    root = tree.getroot()

    for iteration in root.iter('Iteration'):
        for gene in iteration:
            if gene.tag == 'Iteration_query-def':
                g+=1
                cgene = Query(gene.text)
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
                            bioprocs = re.findall(r'\[([a-zA-Z0-9 \-,;/]+)\]', anot)
                            chit.funct=' / '.join(bioprocs)

                            chit.COGname=str(re.search(r'(COG[0-9]+)', anot).group(1))

                            if chit.COGname in cog_categories:
                                chit.COGcat=[*cog_categories[chit.COGname]]
                            elif chit.num==1:
                                chit.COGcat=['NF']

                        for ithsp in hit.iter('Hsp'):
                            for hsp in ithsp:

                                if hsp.tag == 'Hsp_evalue':
                                    chit.evalue = float(hsp.text)
                                elif hsp.tag == 'Hsp_query-from':
                                    qstart = int(hsp.text)
                                elif hsp.tag == 'Hsp_query-to':
                                    qend = int(hsp.text)
                                    chit.cover = round(100*((qend-qstart)+1)/cgene.len,2)
                                elif hsp.tag == 'Hsp_identity':
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
                    output.write('\t'.join([cgene.name, 'na', 'na', 'na', 'na', 'na', 'na']) + '\n')
                    ns+=1
                    n, f=add_to_dict(['NS'], n, f)
                    continue
                
                a+=1
                output.write('\t'.join([cgene.name, best_hit.COGname, str(best_hit.cover), str(best_hit.identity), str(best_hit.evalue), best_hit.funct, cat_line(best_hit.COGcat)]) + '\n')
                if best_hit.evalue<0.001:
                    n, f=add_to_dict(best_hit.COGcat, n, f)
                    if best_hit.COGcat==['NF']:
                        not_found[chit.COGname]=''
                else:
                    ns+=1
                    n, f=add_to_dict(['NS'], n, f)

    for key in ann_dict.keys():
        value=ann_dict[key]
        perc=value*100/g
        output2.write("{}\t{}\t{}\n".format(key,value,perc))
    
    print(key)
    print('specific genes: ', g, ' genes with no hits: ',g-a, ' functions found: ', n, ' genes with a function: ', f)
    print('NS: ', ns)
    return g

def main():

    #import xml file with the RPS-BLAST results
    if len(sys.argv)!=3:
        sys.errot.write('{} <input xml file> <output name>')
        exit(1)

    #create COG-function dictionary from the table of all COGs 
    cog_categories=dict()
    for line in open('/home/mcoin/Petrotoga/functional_annotation/cog-20.def.tab'):
        form=line.split('\t')
        cog_categories[form[0]]=form[1]

    oname=sys.argv[2]
    filename=sys.argv[1]

    output=open('/home/mcoin/Petrotoga/functional_annotation/{}{}.csv'.format(oname, node.text), 'a')
    output.write('\t'.join(['Core gene', 'COG best hit', 'coverage %', 'identity %', 'e-value', 'function', 'category']) + '\n')
    output2=open('/home/mcoin/Petrotoga/functional_annotation/{}{}-summary.csv'.format(oname, node.text), 'a')
    output2.write("category\tn.genes\tpercent\n")

    ann_file=open(filename)

    ann_dict=dict()
    g=0
    not_found=dict()

    g=parserBlastResults(ann_file, g)

    output2.close()
    output.close()

    print('Keys not found:', ','.join(list(not_found.keys())))
