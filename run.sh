#! /bin/sh
snpfile=$1
hmpfile=$2
gff=$3
cds=$4
prefix=$5

python ./SNPInfoExtract.py $snpfile $hmpfile
python ./WhereTheSnp.py ./4WhereTheSNP.txt $gff $cds $prefix
rm -rf ./4WhereTheSNP.txt
