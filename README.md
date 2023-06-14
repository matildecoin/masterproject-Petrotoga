# Masterproject-Petrotoga

In this repository are found the scripts used as part of the comparative genomic analysis of the genus _Petrotoga_.
The analysis was mainly based on preexisting software applications such as stand-alone BLAST [1], Count [2], Prokka [3], and Roary [4].

The project aims to conduct a comparative genomics analysis of the seven cultivated and sequenced _Petrotoga_ species, along with their closest phylogenetic neighbors, to identify genomic signatures associated with their extreme living environment. 

The analysis is divided in three parts: 

i)  Preliminary analyses to define genus and species boundaries

ii) Comparative analysis based on pangenome. Functional analysis of core and specific genes, identification of marker genes.

iii) Reconstruction of the evolutionary history of pathways of interest.


## Specie and Genus determination

The proteome of each species was predicted with Prokka.



### Specie determination
To identify species relationships, I utilized the Average Nucleotide Identity (ANI) metric with BLAST algorithm, specifically \ac{ANI}b, from the web-based software JSpeciesWS. [4]
### Genus determination
To identify genus boundaries I used two methods, one developed by Qin _et al._ 2014 [5] defined as Percent of Conserved Proteins (POCP) and one defined by Barco _et al._ 2020 [6].

POCP: First step is a all agains all BLASTP between all the predicted coding sequences

## Pagenome analysis

## Evolutionary history
