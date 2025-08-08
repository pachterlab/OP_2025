#! /bin/bash

while read p; do
	echo $p
	grep 'u234187' count_matrices_${p}/bus_df.csv | awk '{print $2}' | grep -A 1 -f - ${p}_fastq/${p}_2.fastq > ${p}_extracted_u234187.fastq
done <SRR_F0K.txt