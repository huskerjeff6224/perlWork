#!/bin/bash
#This script lets you run a series of commands then wait for those commands to finish
# before running another set of commands. Good for jobs that depend on the prevcous commands output

FAIL=0

echo "starting"

list=`ls | grep -P 'thout_*'`
for dir in $list
 do
echo $dir
 samtools sort -n $dir/accepted_hits.bam  $dir/accepted_hits_sorted &
 done

for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
    done

    echo $FAIL

    if [ "$FAIL" == "0" ];
    then
    echo "YAY!"
    else
    echo "FAIL! ($FAIL)"
    fi

st=`ls | grep -P 'thout_*'`
for dir in $list
 do
 echo $dir
   samtools view $dir/accepted_hits_sorted.bam > $dir/accepted_hits_sorted.sam
   done

   for job in `jobs -p`
   do
    echo $job
    wait $job || let "FAIL+=1"
   done

   echo $FAIL

   if [ "$FAIL" == "0" ];
      then
        echo "YAY!"
      else
        echo "FAIL! ($FAIL)"
   fi


list=`ls | grep -P 'thout_*'`
for dir in $list
   do
     echo $dir
    /safer/investigator_data/sanger/python/bin/htseq-count -q -s no $dir/accepted_hits_sorted.sam /safer/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > $dir/counts.txt &
     done

    for job in `jobs -p`
     do
     echo $job
     wait $job || let "FAIL+=1"
     done

echo $FAIL

  if [ "$FAIL" == "0" ];
    then
      echo "YAY!"
    else
      echo "FAIL! ($FAIL)"
 fi

