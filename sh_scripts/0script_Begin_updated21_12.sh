#!/bin/bash

#### RUN directly from SCRIPT FOLDER
fecha=$(date --rfc-3339=date ) #Ejecuta date y lo almacena en la variable fecha


logfile=LoganalysisRNA_${fecha}.log

echo it has begun at $(date)

echo .| echo .
#Preparing directory and files...
echo Preparing directories and files...

#######
# Abrir 
#directory=./Files_2019_05_STAR #$ ./SCRIPTS/script_Begin.sh 
directory1=~/rnaseq-pipeline/analisisMar_P0

#$ ./script_Begin.sh #Inside of the folder

#crea 2 directorios, FastQCfiles y filesSTAR. Dentro de fileSTAR crea 3 directorios,BAMfiles,FeatureCounts,TDFfiles.
#Dento de BamQCFiles crea BamQCFiles, logSTAR,toTranscriptome.
mkdir -vp ${directory1}/{FastQCfiles,trimFastQCfiles,TrimFiles,filesSTAR/{BAMfiles/{BamQCFiles,logSTAR,toTranscriptome},FeatureCounts,TDFfiles}}


#${directory}/{FastQCfiles,filesSTAR/{BAMfiles/logSTAR,HTseq,TDFfiles}}

#ruta donde se encuentra el raw data
##### MODIFICAR 
fastqfiles1=~/rnaseq-pipeline/RawData_RNAseq/RawData_Mar/P0/muestrasP0/X204SC21112540-Z01-F001/raw_data/muestras

#version del genoma
genomeversion=104
#tipo de libreria que usará FeatureCounts
libtypeFeatureCounts=0 #-s < intorstring > (isStrandSpecific)  0 (un-stranded), 1 (stranded) and 2 (reversely stranded)
libPE=-p
#https://www.tldp.Corg/HOWTO/Bash-Prompt-HOWTO/x700.html 
#busca los archivos que empiezen por l sitados en el directorio Fastq y luego cuenta el numero de archivos que hay
NofSample1=$(ls -1 $fastqfiles1 | grep -v ^l | wc -l)



########################################
echo .| echo .
#echo FASTQC Quality Control 
#bucle sobre todos los archivos .fq.gz ordenados que no estén repetidos


for sampleraw in $fastqfiles1/*.fq.gz;do
   echo $sampleraw
   fastqc -t 16 $sampleraw -o ${directory1}/FastQCfiles/
done

#Trim-galore
echo "comienza el trim_galore"
output_dir=$directory1/TrimFiles
for file in $fastqfiles1/*_1.fq.gz ; do

#nombre del archivo
sample=$(basename $file|cut -d '_' -f 1)
fastq1=${fastqfiles1}/${sample}_1.fq.gz
fastq2=${fastqfiles1}/${sample}_2.fg.gz
echo ${fastq1}
echo ${fastq2}

#directorio salida


#mkdir -p ${output_dir_sample}

 	trim_galore --paired ${fastq1} ${fastq2} -o ${output_dir}
done

#echo "comienza el FASTQC"
##FastQc a los archivos trimeados
for trimfile in ${directory1}/TrimFiles/*.fq.gz;do

	 fastqc -t 12 $trimfile -o ${directory1}/trimFastQCfiles
	echo "el fastqc ha terminado para $trimfile"
done
#llamamos al programa STARalingment para que ejecute el alineamiento por STAR 
cd
cd rnaseq-pipeline
echo RNA-seq ANALYSIS has begun at $(date)
time ./STARalingment.sh $directory1 $genomeversion |& tee -a $logfile #cambiar si se cambia de alineador
#directory=$1
#fastqfiles=$2tee
#genomeversion=$3
echo STAR ANALYSIS has finished at $(date)

#Perform strand-specific read counting (use '-s 2' if reversely stranded, '-s 1' if forwardly stranded): 
time ./quantSTAR_FeatureCounts.sh $directory1 $genomeversion $libtypeFeatureCounts $libPE |& tee -a $logfile
#directory=$1
#genomeversion=$2
#libtypeFeatureCounts=$3
#NofSample=$4
echo Quantification with Featurecounts has finished at $(date)

echo Kallisto Analysis!!
time ./Kallisto_alignment.sh $directory1/TrimFiles $directory1/kallisto_analisis|& tee -a $logfile
echo \n
echo Salmon Analysis!!
time ./Salmon_alignment.sh $directory1/TrimFiles $directory1/salmon_analisis|& tee -a $logfile

echo Multiqc Report

multiqc ${directory1}/FastQCfiles/ ${directory1}/TrimFiles ${directory1}/trimFastQCfiles ${directory1}/filesSTAR/BAMfiles/logSTAR/ ${directory1}/filesSTAR/FeatureCounts/ ${directory1}/salmon_analisis ${directory1}/kallisto_analisis -o $directory1 

echo "This script is sent to the log file, \"$logfile\"." 

