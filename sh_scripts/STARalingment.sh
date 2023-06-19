#!/bin/bash

start=$(date +%s)

echo .| echo .
#Preparing directory and files...
echo Preparing directories and files...

##### MODIFICAR 
#argumentos de entrada
directory=$1
genomeversion=$2

#NO MODIFICAR

genomeSTAR=/home/LAB/lab_glb_2/rnaseq-pipeline/indexMus_musculus/v${genomeversion} #ruta genoma STAR 
#genomeigv=/home/lpa/igv/genomes/mm39.genome #ruta genoma igv 
gtffile=/home/LAB/lab_glb_2/rnaseq-pipeline/ref-genome/mus_musculus/gtf/v${genomeversion} #ruta para el archivo .gtf del genoma
######################################
echo $genomeSTAR
echo $gtffile
echo .| echo .
#Alingment with STAR
echo Alingment with STAR? Yes, of course\!\!

for sample in $directory/*_1.fq.gz ; do
  echo ${sample}
   ruta=$(dirname $sample)
   describer=$(basename $sample | cut -d '_' -f 1)
   echo ${describer}
   echo $ruta

   	STAR --runThreadN 16 --genomeDir $genomeSTAR --genomeLoad LoadAndKeep --limitBAMsortRAM 200000000000\
      	--readFilesIn $directory/${describer}_1.fq.gz $directory/${describer}_2.fq.gz  \
      	--outFileNamePrefix ${describer}_ \
      	--outSAMtype BAM SortedByCoordinate \
     	--quantMode TranscriptomeSAM GeneCounts \
    	--readFilesCommand zcat
	#ordenamos todos los archivos

done
STAR --genomeDir $genomeSTAR --genomeLoad Remove



	mv -v ./*.tab $directory/filesSTAR/BAMfiles/logSTAR/
	mv -v ./*.out $directory/filesSTAR/BAMfiles/logSTAR/
	mv -v ./*Aligned.sortedByCoord.out.bam $directory/filesSTAR/BAMfiles/
	mv -v ./*.toTranscriptome.out.bam $directory/filesSTAR/BAMfiles/toTranscriptome/


end=$(date +%s)
duration=$((end - start))

echo "El programa ha tarda en ejecutarse ${duration}" > durationSTAR.txt
