#! /bin/bash

ref=/home/oakesc/references/rat
kallisto=/home/oakesc/.local/lib/python3.9/site-packages/kb_python/bins/linux/kallisto/
bustools=/home/oakesc/.local/lib/python3.9/site-packages/kb_python/bins/linux/bustools/

while read p; do
	prefetch "$p"
	fasterq-dump $p -O ${p}_fastq --include-technical
	f1=${p}_fastq/${p}_1.fastq
	f2=${p}_fastq/${p}_2.fastq
	f3=${p}_fastq/${p}_3.fastq
	output_dir=count_matrices_$p
	mkdir -p $output_dir/tmp
	mkdir -p $output_dir
	$kallisto/kallisto bus -x "-1,0,0:0,0,0:1,0,0,2,0,0" -i $ref/index.idx -o $output_dir -t 40 --paired $f1 $f2 $f3 
	$bustools/bustools sort -o $output_dir/tmp/output.s.bus -T $output_dir/tmp -t 40 $output_dir/output.bus
	$bustools/bustools inspect -o $output_dir/inspect.json $output_dir/tmp/output.s.bus

	mkdir -p $output_dir/counts_unfiltered_total

	$bustools/bustools count -o $output_dir/counts_unfiltered_total/cells_x_tcc -g $ref/t2g.txt -e $output_dir/matrix.ec -t $output_dir/transcripts.txt --multimapping $output_dir/tmp/output.s.bus

	mkdir -p $output_dir/quants_unfiltered_total

	$kallisto/kallisto quant-tcc -o $output_dir/quants_unfiltered_total -i $ref/index.idx -e $output_dir/counts_unfiltered_total/cells_x_tcc.ec.txt -g $ref/t2g.txt -t 40 -f $output_dir/flens.txt --matrix-to-directories $output_dir/counts_unfiltered_total/cells_x_tcc.mtx

	rm -rf $output_dir/tmp
	rm -r $p
	rm -r ${p}_fastq
done <SRR_list.txt
