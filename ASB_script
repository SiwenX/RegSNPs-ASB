#!/bin/bash

n=0
cat SNP_in_TFBS.bed | while read line
do
	line_arr=( $line )
	pos=${line_arr[0]}
	sta=${line_arr[1]}
	end=${line_arr[2]}	
	samtools view -F 16 path/to/.bam  $pos:$sta-$end | awk -v OFS="\t" '{print $3,$4,$4+length($10)-1,$10,$6,$4+5}' > fragments_for.bed
	grep -v "S" fragments_for.bed > fragments_for_noS.bed
	grep -v "*" fragments_for_noS.bed > fragments_for_noSS.bed
	samtools view -f 16 path/to/.bam  $pos:$sta-$end | awk -v OFS="\t" '{print $3,$4,$4+length($10)-1,$10,$6,$4+length($10)-1-4}' > fragments_rev.bed
	grep -v "S" fragments_rev.bed > fragments_rev_noS.bed
	grep -v "*" fragments_rev_noS.bed > fragments_rev_noSS.bed
	cat fragments_for_noSS.bed fragments_rev_noSS.bed > fragments_noSS.bed	
	((n=n+1))
	echo $n
	Rscript ASB.R $pos $sta $n

done


