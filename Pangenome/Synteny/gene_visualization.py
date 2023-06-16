#!/usr/bin/env python3

import re, os
import pandas as pd

node_dir='OneDrive\Desktop\GhostKoala\Pathwayall'
node_files=os.listdir(node_dir)

path=dict()
node_list=list()
node_list.append('GeneralCat')
node_list.append('SpecificCat')

for node_file in node_files:
    file=open('{}\{}'.format(node_dir, node_file))
    node=re.search('node(.+).txt', node_file).group(1)
    node_list.append(node)

    l=0
    for line in file:
        if not line.startswith('M0'):
            if l>0:
                cat=prev_l.rstrip()
                prev_l=line
            else:
                prev_l=line
                l+=1
            
        if line.startswith('M0'):
            if l>0:
                scat=prev_l.rstrip()
            l=0
            search=re.search('M\d{5} (.+) \(\d+\)  \(.+ (\d+)\/(\d+)\)', line)
            if not search:
                print(line)
            name=search.group(1)
            num=int(search.group(2))
            den=int(search.group(3))
            r=num/den

            if name not in path:
                path[name]=dict()
                
            for n in node_list:
                if n not in path[name]:
                    if n=='GeneralCat':
                        path[name][n]=cat
                        print( name,'g', cat)
                    elif n=='SpecificCat':
                        path[name][n]=scat
                        print('s', scat)
                    else:
                        path[name][n]=0.0
            path[name][node]=r

        

for k in list(path.keys()): 
    if node not in path[k]:
        path[k][node]=0.0

df=pd.DataFrame.from_dict(path)
df.T.to_csv('OneDrive\Desktop\GhostKoala\paths_all_petrotogaceae.csv', sep='\t')

#print(path)
''''
csv=''
for k1 in list(path.keys())[0]:
    csv+=' \t'
    for k2 in list(path[k1].keys()):
        csv+='{}\t'.format(k2)
    csv.rstrip()
    csv+='\n'

for k1 in list(path.keys()):
    csv+=''.format(k1)
    for k2 in list(path[k1].keys()):
        csv+='{}\t'.format(path[k1][k2])
    csv.rstrip()
    csv+='\n'

print(csv)
'''''
