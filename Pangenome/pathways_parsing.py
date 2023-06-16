#!/usr/bin/env python3

import re, os
import pandas as pd

node_dir='OneDrive\Desktop\GhostKoala\Pathwayall'
node_files=os.listdir(node_dir)

path=dict()
node_list=list()
#adds the name of the first 2 columns of the table
#broad metabolic map to which the pathway belongs (e.g. Aminoacid biosynthesis)
node_list.append('GeneralCat')
#more specific metabolic map (e.g. Proline biosynthesis)
node_list.append('SpecificCat')


for node_file in node_files:
    file=open('{}\{}'.format(node_dir, node_file))
    node=re.search('node(.+).txt', node_file).group(1)
    node_list.append(node)
    #adds a specie to the table

    l=0
    #represents the number of previous lines that didn't start with 'M0', a pathway
    for line in file:
        if not line.startswith('M0'):
            
            if l>0:
                #if the previous line was not a pathway saves the content of the previous line as a category
                cat=prev_l.rstrip()
                prev_l=line
            else:
                #save the line if it's not starting with 'M0'
                prev_l=line
                l+=1
            
        if line.startswith('M0'):
            if l>0:
                scat=prev_l.rstrip()
                #specific category is always before the firts pathway line
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

        
#fill with 0 all the cells that have been missed
for k in list(path.keys()): 
    if node not in path[k]:
        path[k][node]=0.0

df=pd.DataFrame.from_dict(path)
df.T.to_csv('OneDrive\Desktop\GhostKoala\paths_all_petrotogaceae.csv', sep='\t')

