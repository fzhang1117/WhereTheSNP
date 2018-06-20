
WhereTheSNP is a python written script set that helps to predict the function of SNP.

Usage:
---
This script should be run in linux envrionment. Python is needed.
```
run.sh snpfile hmpfile gff cds prefix
```
Input:
---
snpfile:

The list of SNPs you are interested

hmpfile:

A hapmap format file of these snps, or a larger dataset

gff:

A gff file the orgnisms

cds:

A fasta format file of the coding sequences, the cds file must corresponded with gff file in version

prefix:

Prefix of the output. 

output:
---
prefix_cds.txt：

A file contains snps and relevant informations in coding regions.

prefix_not_in_cds.txt：

A file contains snps and relevant informations in gene region but in in coding regions.

prefix_promoter_region.txt：

A file contains snps and relevant informations in less 2000bp upstream of a gene.

prefix_others.txt：

A file contains snps which are not in gene regions and promoter regions.
