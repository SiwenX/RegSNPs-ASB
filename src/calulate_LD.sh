#! /bin/bash
# The 1000 Genomes vcf files was downloaded at: hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
# Then convert vcf file to plink file using: plink2 --vcf .vcf --out "plink_chr*", we can get .pgen/.psam/.pvar

n=0
cat gene.bed | while read line
do
line_arr=( $line )
chr=${line_arr[0]}
id1=${line_arr[1]}
id2=${line_arr[2]}
plink2 --pfile "plink2_$chr" --ld 'id1' 'id2'
grep --after-context=6  "alleles:" plink2.log > $id2

done
