#!/bin/bash
#QUANTIFICATION

#Preparing directory and files...
echo Preparing directories and files...

directory=$1
genomeversion=$2
libtypeFeatureCounts=$3
libPE=$4

#validacion del numero de argumentos
if [ $# -eq 0 ]; then
   echo "No se ha utilizado ningún parámetro"
   exit
elif [ $# -ne 4 ]; then
   echo "Not equal to 4 parameters"
   exit
else
   while [ $# -gt 0 ]; do
      echo $1
      shift
   done
fi
#NO MODIFICAR
gtffile=/home/LAB/lab_glb_2/rnaseq-pipeline/ref-genome/mus_musculus/gtf/v${genomeversion}/*.gtf
echo ${gtffile}
############## SubReads -> FeatureCounts
echo .| echo .
echo featureCounts...
#Without filter by MAPQ > 0
 
#featureCounts
featureCounts $libPE -s $libtypeFeatureCounts -T 6 -t exon -g gene_id -a $gtffile \
-o ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt ${directory}/filesSTAR/BAMfiles/*.bam

#-s < intorstring > (isStrandSpecific)  0 (un-stranded), 1 (stranded) and 2 (reversely stranded)


#https://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/x700.html 
#Counting Files in the Current Directory
NofSample=$(ls -1 ${directory}/filesSTAR/BAMfiles/*.bam | grep -v ^l | wc -l)
echo "Number of BAM files = " $NofSample

cut -f 1,$(seq -s , 7 $((7+$(($NofSample - 1))))) ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt > ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.mat 







