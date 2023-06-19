#!/bin/bash
start=$(date +%s)
# Set the paths to input files and Salmon executable
fastq_dir=$1
output_dir=$2
mkdir -p $output_dir
#kallisto index
genome_index=/home/LAB/lab_glb_2/rnaseq-pipeline/salmon_GRC39_V104_index

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
  salmon quant\
    -i ${genome_index} \
    -l A \
    -p 16 \
    -o ${output_dir_sample} \
    --validateMappings \
    --gcBias\
    --posBias\
    -1 ${fastq1} -2 ${fastq2} |& tee -a /${output_dir_sample}/${sample}_logfile.log 

done 

echo "All samples processed."
end=$(date +%s)
duration=$((end-start))

echo "Salmon ha tardado en alinear las muestras ${duration}">${output_dir}/duracion_salmon.txt
