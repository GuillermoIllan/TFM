library(readr)
library(tibble)
library(dplyr)
resultados_alineadores <- read.csv("Escritorio/datosrecuperados/scripts/Rscripts/AED/P0/P0_analisis/resultados_alineadores.csv")


resultados_alineadores <- lapply(resultados_alineadores, function(x) gsub("%", "", x))
resultados_alineadores<-as.data.frame(resultados_alineadores)
resultados_alineadores <- lapply(resultados_alineadores[,-1], function(x) as.numeric(x))

star_align_media<-mean(resultados_alineadores$STAR...Aligned)
salmon_align_media<-mean(resultados_alineadores$SALMON...Aligned)
kal_align_media<-mean(resultados_alineadores$KALLISTO...Aligned)

star_reads_media<-mean(resultados_alineadores$STAR.M.Aligned)
salmon_reads_media<-mean(resultados_alineadores$SALMON.M.Aligned)
kal_reads_media<-mean(resultados_alineadores$KALLISTO.M.Aligned)
