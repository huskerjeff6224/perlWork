#!/bin/sh
#loops though directory to find directories to loop through and run a command

list=`ls | grep -P 'thout_*'`
for dir in $list
do
echo $dir
samtools sort -n $dir/accepted_hits.bam  $dir/accepted_hits_sorted &
done
#samtools view $dir/accepted_hits_sorted.bam > $dir/accepted_hits_sorted.sam
#/safer/investigator_data/sanger/python/bin/htseq-count -q -s no $dir/accepted_hits_sorted.sam /safer/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > $dir/counts.txt 
