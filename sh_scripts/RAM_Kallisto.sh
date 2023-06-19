#!/bin/bash

# Ruta al ejecutable de STAR
#EXECUTABLE="kallisto"

# Ruta al directorio del genoma
GENOME_DIR="/home/LAB/lab_glb_2/rnaseq-pipeline/ref-genome/mus_musculus/indexKALLISTO/Mus_musculus.GRCm39.cdna.all.104.idx"

# Ruta al archivo FASTQ de muestra
SAMPLE_FASTQ_1="/home/LAB/lab_glb_2/rnaseq-pipeline/RawData_RNAseq/RawData_Mar/X204SC22101701-Z01-F003/01.RawData/CA81/CA81_FKRN230005364-1A_HFWVTDSX5_L3_1.fq.gz"
SAMPLE_FASTQ_2="/home/LAB/lab_glb_2/rnaseq-pipeline/RawData_RNAseq/RawData_Mar/X204SC22101701-Z01-F003/01.RawData/CA81/CA81_FKRN230005364-1A_HFWVTDSX5_L3_2.fq.gz"
describer=$(basename $SAMPLE_FASTQ_1 | cut -d '_' -f 1-4)
# Opciones de STAR
#OPTIONS= "quant\
 #	-i $GENOME_DIR\
#	-t 16\
#	 -o $describer\
#	 $SAMPLE_FASTQ_1 $SAMPLE_FASTQ_2"
# Función para obtener el uso de RAM en megabytes
get_ram_usage() {
  local pid=$1
  local ram_usage=$(pmap -x $pid | tail -1 | awk '{print $3}')
  echo $((ram_usage / 1024))
}
#toma un ID de proceso (pid) como argumento. Utiliza el comando pmap para obtener información detallada sobre el uso de memoria del proceso y luego extrae el valor de la memoria en la última línea. Luego, divide ese valor por 1024 para obtener el uso de RAM en megabytes

# Archivo de salida para la información de RAM
OUTPUT_FILE="uso_ram_kallisto.txt"

# Ejecutar STAR y medir el uso de RAM
start_time=$(date +%s)
kallisto quant -i $GENOME_DIR -t 16 -o $describer $SAMPLE_FASTQ_1 $SAMPLE_FASTQ_2 & # & al final hace que se ejecute en segundo plano.
pid=$! #El ID del proceso se guarda en la variable star_pid utilizando el $! que representa el ID del proceso más reciente en segundo plano

while kill -0 $pid 2>/dev/null; do
  ram_usage=$(get_ram_usage $pid)
  echo "Uso de RAM actual: $ram_usage MB"
  echo "Uso de RAM actual: $ram_usage MB" >> $OUTPUT_FILE
  sleep 10
done

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "El proceso de alineación ha finalizado."
echo "Tiempo de ejecución: ${execution_time} segundos."
echo "Tiempo de ejecución: ${execution_time} segundos." >> $OUTPUT_FILE

