#!/bin/bash
start=$(date +%s)
# Set the paths to input files and Kallisto executable
fastq_dir=$1
output_dir=$2
mkdir -p $output_dir
#kallisto index
genome_index=/home/LAB/lab_glb_2/rnaseq-pipeline/ref-genome/mus_musculus/indexKALLISTO/Mus_musculus.GRCm39.cdna.all.104.idx

# Iterate through each pair of FASTQ files in the input directory
for file in ${fastq_dir}/*_1.fq.gz ; do
  # Extract sample name from file name
  sample=$(basename $file | cut -d '_' -f 1)

  # Define input file paths
  fastq1=${fastq_dir}/${sample}_1.fq.gz
  fastq2=${fastq_dir}/${sample}_2.fq.gz
echo ${fastq1}
echo ${fastq2}
  # Define output file paths
  output_dir_sample=${output_dir}/${sample}
 
echo $output_dir_sample
  # Create output directory if not exists
  mkdir -p ${output_dir_sample}

  # Run Kallisto pseudoalignment with RAM usage monitoring
  kallisto quant\
    -i ${genome_index} \
    -t 16 \
    -o ${output_dir_sample} \
    ${fastq1} ${fastq2} |& tee -a /${output_dir_sample}/${sample}_logfile.log 

mapping_file=~/rnaseq-pipeline/analisisMar2/kallistoAnalisis/${sample}/transcripts.txt


  # Print completion message
  echo "Pseudoalignment complete for sample ${sample}. Results stored in ${abundance_file}."
done 

echo "All samples processed."
end=$(date +%s)
duration=$((end-start))

echo "Kallisto ha tardado en alinear las muestras ${duration}">${output_dir}/duracion_kallisto.txt
