#Analisis P0 LFC0.137_padj01#############################################################
###Preparación de los datos#####################################################
library("AnnotationDbi")
library("org.Mm.eg.db")
library("clusterProfiler")
library("DESeq2")
library("gplots")
library("ggplot2")
library("biomaRt")
library("genefilter")
library("ggrepel")
library("pheatmap")
library("ashr")
library("ggVennDiagram")
library("dplyr")
library("ComplexHeatmap")
library("readr")
library("tibble")
library("apeglm")
library("GeneOverlap")
library("plotly")
library("openxlsx")
library("gplots")
library("circlize")
library("eulerr")
setwd("~/Escritorio/datosrecuperados/scripts/Rscripts/AED")

cts <- as.matrix(read.csv("featurecountstotal_STAR_2.tsv",sep="",row.names="Geneid"))
coldata <- read.csv("metaDataMar.csv",sep="," ,row.names=1)
coldata <- coldata[,c("Tratamiento","Region")]
coldata$Tratamiento <- factor(coldata$Tratamiento)
coldata$Region <- factor(coldata$Region)

#library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Region + Tratamiento)
###PRE-FILTERING################################################################
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
###COOKS DISTANCE###############################################################
#par(mar=c(8,5,2,2))
#boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

###log transformation###########################################################
rld <- rlog(dds)
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)

plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
par(mfrow=c(1,1))

###Distance#####################################################################
sampleDists <- dist( t( assay(rld) ) )
sampleDists
###Heatmap######################################################################
#library("gplots")
#library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Tratamiento,rld$Region, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE, main= "Distancias")

##PCA##########################################################################
data<-plotPCA(rld, intgroup = c("Tratamiento","Region"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, shape=Region, data=data, label=colnames(rld)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "figuras/PCA_CAUTvsCTRL_general_L2FC0137_01padj.png",
    width=500, height = 500)
p + geom_text(aes(label=colnames(rld)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) + theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(legend.title = element_text(size=15))+
  theme(legend.text = element_text(size=12))+
  guides(color = guide_legend(override.aes = list(size = 10)),shape = guide_legend(override.aes = list(size = 10)))+
  theme(plot.background = element_rect(fill = "white"))+
  theme(panel.background = element_rect(fill = 'white'))+
  geom_point(size = 6)
dev.off()



rv <- rowVars(assay(rld))
# select the ntop genes by variance
# select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld)[select,]))
summary(pca)
var_exp <- summary(pca)$importance[2,]
cum_var_exp <- cumsum(var_exp)
#elbowplot varianza
plot(var_exp, type="b", xlab="Número de componentes principales", ylab="Varianza explicada")

#elbowplot varianza acumulada
plot(cum_var_exp, type="b", xlab="Número de componentes principales", ylab="Varianza explicada acumulada")


#Differential expression analysis#############################################
dds$Tratamiento <- relevel(dds$Tratamiento, "Control")
dds <- DESeq(dds)
res<- results(dds, alpha = 0.1, lfcThreshold= 0)
resShr<- lfcShrink(dds,res = res, coef="Tratamiento_Cauterizado_vs_Control")#Tabla resultado
#Tabla resultado
dim(resShr)
head(resShr)
summary(resShr)
sigs <- na.omit(resShr)
sigs<- sigs[sigs$padj < 0.1,]
df <- as.data.frame(sigs)
dim(df)

###Anotación con biomaRt#######################################################

ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                         host =  'https://may2021.archive.ensembl.org')

mouse <- useDataset("mmusculus_gene_ensembl", ensembl_genes)
listEnsemblArchives()
# ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "104")
rownames_ordenados <- sort(rownames(df))
datos_ordenados <- df[rownames_ordenados, ]
values <- rownames(datos_ordenados)
# lista_atributos<-listAttributes(ensembl)
tabla_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(df),
                       mart = mouse)

dim(tabla_anotacion)
dim(df)
# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combined <- merge(df, tabla_anotacion %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_final<- unique(df_combined)
rownames(df_final)<-df_final$Row.names
df_final$Row.names<-NULL
head(df_combined)

###Heatmap#####################################################################
#filtramos los genes para quedarnos con los significativos
df.top <- df_final[(abs(df_final$log2FoldChange) > 0.137),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
#numero de upregulados=2275 y downregulados=2230
df.top$diffexpressed <- "NO"
df.top$diffexpressed[df.top$log2FoldChange > 0.137] <- "UP"
df.top$diffexpressed[df.top$log2FoldChange < -0.137 ] <- "DOWN"
table(df.top$diffexpressed)
#creamos la matriz
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata)] #sig genes x samples
rownames(mat)<- df.top$mgi_symbol
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- paste0(rld$Tratamiento,"-",rld$Region)

#heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Tratamiento ]

heatmap.2(mat.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row",main = "embWPC_vs_Ctrl_2682up_2662down")

###Enriquecimiento#############################################################

#Filtramos los significativos

genes_to_test <- rownames(df.top)
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 15))
png(filename = "dotplot_CAUTvsCTRL_general_L2FC0137_01padj.png",
    width=1000, height = 1000)
fit2<- plot(dotplot(GO_results, showCategory = 15))
dev.off()

###Volcano plot con ggplot#######################################################
head(df)

#seleccionamos los genes expresados diferencialmente
df_final$diffexpressed <- "NO"
df_final$diffexpressed[df_final$log2FoldChange > 0.137 & df_final$padj < 0.1] <- "UP"
df_final$diffexpressed[df_final$log2FoldChange < -0.137 & df_final$padj < 0.1] <- "DOWN"

#etiquetamos los genes diferenciados
df_final$delabel <- NA
df_final$delabel[df_final$diffexpressed != "NO"] <- df_final$mgi_symbol[df_final$diffexpressed != "NO"]
table(df_final$diffexpressed)

#representamos con ggplot
png(filename = "volcano_CAUTvsCTRL_general_L2FC0137_01padj.png.png",
    width=1000, height = 1000)
ggplot(data=df_final, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.137, 0.137), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off
#Controles######################################################################
#Subset de las muestras en la matriz
controles<-cts[,1:8]
controles_coldata<-coldata[1:8,]

#library("DESeq2")
ddsC <- DESeqDataSetFromMatrix(countData = controles,
                               colData = controles_coldata,
                               design = ~ Region)

ddsC
ddsC<-estimateSizeFactors(ddsC)
keepC <- rowSums(counts(ddsC)) >= 10
ddsC <- ddsC[keepC,]
ddsC <- DESeq(ddsC)
resC<-results(ddsC, alpha = 0.1, lfcThreshold= 0)
summary(resC)
head(resC)
resC_apglm<- lfcShrink(ddsC,res = resC, coef = "Region_PMBSF_vs_ALBSF")#Tabla resultado
summary(resC_apglm)
head(resC_apglm)
rldC<-rlog(ddsC)

##PCA#################
data<-plotPCA(rldC, intgroup = "Region",returnData=TRUE,
              ntop = 500)


percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Region, data=data, label=colnames(rldC)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p + geom_text(aes(label=colnames(rldC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Region),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75) 

ggsave("controlesPCA.png", device="png", path = "./AED-Mar_files/graficas-reunionMar-20-03")

##Filtrar los datos#################
sigsC <- na.omit(resC_apglm)
sigsC<- sigsC[sigsC$padj < 0.1,]
sigsC
dfC <- as.data.frame(sigsC)
dim(dfC)

##plotMA#################
plotMA(resC_apglm, ylim=c(-5,5))

##Anotación con biomaRt############
#library(biomaRt)
#ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "104")
rownames_ordenadosC <- sort(rownames(dfC))
datos_ordenados_C <- dfC[rownames_ordenadosC, ]
valuesC <- rownames(datos_ordenados_C)
#lista_atributos<-listAttributes(ensembl)
tabla_anotacionC<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = valuesC,
                        mart = mouse)

# Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combinedC <- merge(dfC, tabla_anotacionC %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_finalC<- unique(df_combinedC)
rownames(df_finalC)<-df_finalC$Row.names
df_finalC$Row.names<-NULL
df_final_Tissue<-rownames_to_column(df_finalC,var = "EnsemblID")

#etiquetado de genes
tabla_anotacionC_Tissue <- tabla_anotacionC %>%
  mutate(region = ifelse(ensembl_gene_id %in% genes_comunes_PMBSF_P0_P8_lor , "PMBSFspc:ALBSFembWPC P8&P0",
                         ifelse(ensembl_gene_id %in% genes_comunes_ALBSF_P0_P8_lor, "ALBSFspc:ALBSFembWPC P8&P0",
                                ifelse(ensembl_gene_id %in% genesPMBSF_ALBSF_emb_WPC, "PMBSFspc:ALBSFembWPC P8",
                                  ifelse(ensembl_gene_id %in% genesALBSF_ALBSF_emb_WPC, "ALBSFspc:ALBSFembWPC P8",
                                    ifelse(ensembl_gene_id %in% genes_comunes_PMBSFspc_P0_P8_lor , "PMBSFspc P8&P0",
                                       ifelse(ensembl_gene_id %in% genes_comunes_ALBSFspc_P0_P8_lor, "ALBSFspc P8&P0",
                                              ifelse(ensembl_gene_id %in% genesPMBSFspc_P8, "PMBSFspc P8",
                                                     ifelse(ensembl_gene_id %in% genesALBSFspc_P8, "ALBSFspc P8","Neg")))))))))
write.xlsx(tabla_anotacionC_Tissue, file = "Tissue_shApeglm.xlsx", rowNames = FALSE)
df_final_Tissue<-df_final_Tissue[,-c(8,9)]
table(tabla_anotacionC_Tissue$region)

##Heatmap####################
#library("genefilter")
#filtramos los genes para quedarnos con los significativos
df.topC <- df_finalC[(abs(df_finalC$log2FoldChange) > 0.137),]
df.topC <- df.topC[order(df.topC$log2FoldChange, decreasing = TRUE),]
df.topC$diffexpressed <- "NO"
df.topC$diffexpressed[df.topC$log2FoldChange > 0.137 ] <- "UP"
df.topC$diffexpressed[df.topC$log2FoldChange < -0.137 ] <- "DOWN"
table(df.topC$diffexpressed)

genesPMBSFspc_P8<-rownames(df.topC)[which(df.topC$diffexpressed == "UP")]
genesALBSFspc_P8<-rownames(df.topC)[which(df.topC$diffexpressed == "DOWN")]
#creamos la matriz
rlog_outC <- rlog(ddsC, blind=FALSE) #get normalized count data from dds object
matC<-assay(rlog_outC)[rownames(df.topC), rownames(controles_coldata)] #sig genes x samples
rownames(matC)<- df.topC$mgi_symbol 
base_mean <- rowMeans(matC)
head(matC)
matC.scaled <- t(apply(matC, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(matC.scaled) <- paste0(rlog_outC$Tratamiento,"-",rlog_outC$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rlog_outC$Region ]
heatmap.2(matC.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row", main= "ALBSFctrl_vs_PMBSFctrl_431up_vs_284down")

##Enriquecimiento###################
genes_to_test <- rownames(df.topC)
GO_resultsC <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_resultsC)
fitC <- plot(barplot(GO_resultsC, showCategory = 15))
png(filename = "dotplot_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0137_01padj.png",
    width=1000, height = 1000)
fitC2<- plot(dotplot(GO_resultsC, showCategory = 15))
dev.off()
##Volcano plot con ggplot#######################

#seleccionamos los genes expresados diferencialmente
df_finalC$diffexpressed <- "NO"
df_finalC$diffexpressed[df_finalC$log2FoldChange > 0.137 & df_finalC$padj < 0.1] <- "UP"
df_finalC$diffexpressed[df_finalC$log2FoldChange < -0.137 & df_finalC$padj < 0.1] <- "DOWN"
table(df_finalC$diffexpressed)

#etiquetamos los genes diferenciados
df_finalC$delabel <- NA
df_finalC$delabel[df_finalC$diffexpressed != "NO"] <- df_finalC$mgi_symbol[df_finalC$diffexpressed != "NO"]

#representamos con ggplot
png(filename = "./figuras/volcano_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0137_01padj.png",
    width=700, height = 700)
ggplot(data=df_finalC, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(min.segment.length = 0, size=7) +
  scale_color_manual(values=c("#00a4a2", "grey","#b05d95")) +
  geom_vline(xintercept=c(-0.137, 0.137), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(panel.background = element_rect(fill="white"))
dev.off()

##Complex Heatmap##############################################################################
dfC_anot<-getBM(attributes =c("ensembl_gene_id","mgi_symbol","gene_biotype") ,
                filters = "ensembl_gene_id",
                values = rownames(dfC),
                mart = mouse)
dfc2<-rownames_to_column(dfC,var = "ensembl_gene_id")
dfC_biotype <- inner_join(dfc2,dfC_anot,"ensembl_gene_id")
dfC_filtrado <- dfC_biotype %>% filter(gene_biotype == "protein_coding")

top_25_genesC <- arrange(dfC_filtrado, desc(log2FoldChange)) %>% head(25)
low_25_genesC <- dfC_filtrado[order(dfC_filtrado$log2FoldChange), ]
low_25_genesC <- head(low_25_genesC, 25)
top_50_genesC <- rbind(top_25_genesC,low_25_genesC)

rlog_outC <- rlog(ddsC, blind=FALSE) #get normalized count data from dds object
matC<-assay(rlog_outC)[rownames(df.topC), rownames(controles_coldata)] #sig genes x samples
rownames(matC)<- df.topC$mgi_symbol 
base_mean <- rowMeans(matC)
matC.scaled <- t(apply(matC, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(matC.scaled) <- paste0(rlog_outC$Tratamiento,"-",rlog_outC$Region)

genes_anotacion_C<-match(rownames(matC.scaled), top_50_genesC$mgi_symbol)

# Crear el vector de posiciones
posiciones_matrizC <- which(!is.na(genes_anotacion_C))


label_genesC <- vector("character", length(genes_anotacion_C))
for (i in 1:length(genes_anotacion_C)) {
  if (!is.na(genes_anotacion_C[i])) {
    label_genesC[i] <- top_50_genesC$mgi_symbol[genes_anotacion_C[i]]
  }
}
label_genesC<-subset(label_genesC, nzchar(label_genesC))
print(label_genesC)
anotaciones_c = rowAnnotation(foo = anno_mark(at = posiciones_matrizC,
                                              labels =label_genesC))

col_fun= colorRamp2(seq(-1.5,1.5, length.out=100),
                    (colorRampPalette(c("darkorchid1","black","yellow1"))(100)),
                    space="RGB")
pdf("figuras/Heatmap_ALBSF_ctrl_VS_PMBSF_ctrl_01padj.pdf",
    width=10, height = 10)
Heatmap(matC.scaled,
        col=col_fun,
        na_col = "black",
        column_title = "PMBSFspc(431)_vs_ALBSFspc(284)",
        column_title_gp = gpar(fontsize = 18),
        right_annotation = anotaciones_c,
        row_names_gp = gpar(fontsize = 14),
        show_row_names = FALSE,
        cluster_columns = TRUE,)
dev.off()

#ALBSF_EMB_VS_ALBSF_ctrl########################################################
ALBSF<-cts[,1:4]
cauterizados<-cts[,9:12]
ALBSF_cts<-cbind(ALBSF,cauterizados)

ALBSF_coldat<-coldata[1:4,]
coldata_cauterizados<-coldata[9:12,]
ALBSF_coldata<-rbind(ALBSF_coldat,coldata_cauterizados)
ALBSF_coldata
#library("DESeq2")
ddsA <- DESeqDataSetFromMatrix(countData = ALBSF_cts,
                               colData = ALBSF_coldata,
                               design = ~ Tratamiento)

ddsA
ddsA<-estimateSizeFactors(ddsA)
keepA <- rowSums(counts(ddsA)) >= 10
ddsA <- ddsA[keepA,]
rldA <- rlog(ddsA)
##PCA###########################################################################
data<-plotPCA(rldA, intgroup = "Tratamiento",returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, data=data, label=colnames(rldA)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p + geom_text(aes(label=colnames(rldA)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75
) 

##Expresión diferencial#########################################################
ddsA$Tratamiento <- relevel(ddsA$Tratamiento, "Control")
ddsA <- DESeq(ddsA)
resA<-results(ddsA,alpha=0.1, lfcThreshold= 0)
resA_apeglm<- lfcShrink(ddsA,res=resA,coef = "Tratamiento_Cauterizado_vs_Control")#Tabla resultado
summary(resA_apeglm)

#Filtrar los datos
sigsA <- na.omit(resA_apeglm)
sigsA<- sigsA[sigsA$padj < 0.1,]
dfA <- as.data.frame(sigsA)
dim(dfA)
summary(sigsA)
#plotMA
plotMA(resA_apeglm, ylim=c(-5,5))

#Anotación con biomaRt
rownames_ordenadosA <- sort(rownames(dfA))
datos_ordenados_A <- dfA[rownames_ordenados, ]
valuesA <- rownames(datos_ordenados_A)
tabla_anotacionA<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = valuesA,
                        mart = mouse)

# Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combinedA <- merge(dfA, tabla_anotacionA %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)
#eliminamos duplicados
df_finalA<- unique(df_combinedA)
rownames(df_finalA)<-df_finalA$Row.names
df_finalA$Row.names<-NULL

##Heatmap#######################################################################

#filtramos los genes para quedarnos con los significativos
df.topA <- df_finalA[(abs(df_finalA$log2FoldChange) > 0.137),]
df.topA <- df.topA[order(df.topA$log2FoldChange, decreasing = TRUE),]
df.topA$diffexpressed <- "NO"
df.topA$diffexpressed[df.topA$log2FoldChange > 0.137 ] <- "UP"
df.topA$diffexpressed[df.topA$log2FoldChange < -0.137] <- "DOWN"
table(df.topA$diffexpressed)

#creamos la matriz
rlog_outA <- rlog(ddsA, blind=FALSE) #get normalized count data from dds object
matA<-assay(rlog_outA)[rownames(df.topA), rownames(ALBSF_coldata)] #sig genes x samples
rownames(matA)<- df.topA$mgi_symbol 
base_mean <- rowMeans(matA)
head(matA)
matA.scaled <- t(apply(matA, 1, scale)) #center and scale each column (Z-score) then transpose
matA.scaled
colnames(matA.scaled) <- paste0(rlog_outA$Tratamiento,"-",rlog_outA$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rlog_outA$Tratamiento ]
heatmap.2(matA.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row",main= "ALBSFctrl_vs_ALBSFemb_WPC_1827up_1330down")

##Enriquecimiento##########################################################
genes_to_test <- rownames(df.topA)
GO_resultsA <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results)
fitA <- plot(barplot(GO_resultsA, showCategory = 15))
png(filename = "dotplot_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0137_01padj.png",
    width=1000, height = 1000)
fitA2<- plot(dotplot(GO_resultsA, showCategory = 15))
dev.off()

##Volcano plot con ggplot################################################

#seleccionamos los genes expresados diferencialmente
df_finalA$diffexpressed <- "NO"
df_finalA$diffexpressed[df_finalA$log2FoldChange > 0.137 ] <- "UP"
df_finalA$diffexpressed[df_finalA$log2FoldChange < -0.137] <- "DOWN"
table(df_finalA$diffexpressed)

#etiquetamos los genes diferenciados
df_finalA$delabel <- NA
df_finalA$delabel[df_finalA$diffexpressed != "NO"] <- df_finalA$mgi_symbol[df_finalA$diffexpressed != "NO"]

#representamos con ggplot
png(filename = "volcano_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0137_01padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalA, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.137, 0.137), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

dev.off()

#Genes específicos modificados por la cauterización
##Diagrama de Venn##############################################################
PMBSF_SpC <- resC_apglm[which(resC_apglm$padj <= 0.1 & 
                                resC_apglm$log2FoldChange > 0.137 ),]
nrow(PMBSF_SpC)

ALBSF_SpC <- resC_apglm[which(resC_apglm$padj <= 0.1 & 
                                resC_apglm$log2FoldChange < -0.137 ),]
nrow(ALBSF_SpC)

ALBSF_mbWPC <- resA_apeglm[which(resA_apeglm$padj <= 0.1 & 
                                   abs(resA_apeglm$log2FoldChange) > 0.137 ),]
nrow(ALBSF_mbWPC)
saveRDS(rownames(ALBSF_mbWPC),"genes_ALBSF_embWPC_P8_0137lfc_01padj.rds")

#library(ggVennDiagram)
ven1 <-list(ALBSFspc=rownames(ALBSF_SpC), 
            PMBSFspc=rownames(PMBSF_SpC), 
            ALBSF_emb_WPC = rownames(ALBSF_mbWPC))

#Fisher's test
go.objALBSF_SpC <- newGeneOverlap(rownames(ALBSF_SpC),
                                  rownames(ALBSF_mbWPC),
                                  nrow(cts) )
go.objALBSF_SpC <- testGeneOverlap(go.objALBSF_SpC)
print(go.objALBSF_SpC)
go.objPMBSF_SpC<- newGeneOverlap(rownames(PMBSF_SpC),
                                 rownames(ALBSF_mbWPC),
                                 nrow(cts) )
go.objPMBSF_SpC <- testGeneOverlap(go.objPMBSF_SpC)
print(go.objPMBSF_SpC)

#DiagVenn
diag_venn1=ggVennDiagram(ven1,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P8 ALBSF_embWPC") + 
  theme(text = element_text(size = 16)) + scale_fill_gradient(low="blue",high = "red")

png(filename = "./figuras/diagVenn_Overlapping_genes_ALBSF_embWPC.png",
    width=1000, height = 1000)
diag_venn1
dev.off()
pdf("./figuras/diagVenn_Overlapping_genes_ALBSF_embWPC.pdf",
    width=15, height = 15)
diag_venn1
dev.off()

pdf(file ='./figuras/venn_genesoverALBSFembWPClfcP8.pdf', width = 8, height = 7)
plot(euler(ven1, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 16,
     legend = F)
dev.off()

png(filename ='./figuras/venn_genesoverALBSFembWPClfcP8.png', width = 600, height =600)
plot(euler(ven1, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 16,
     legend = F)
dev.off()

###Lista Overlapping genes#######################################################
library(gplots)
ven_intersect <- venn(ven1)
genes_especificos<-attr(ven_intersect,"intersections")
genesPMBSF_ALBSF_emb_WPC<-genes_especificos$`PMBSF:ALBSF_emb_WPC`
genesALBSF_ALBSF_emb_WPC<-genes_especificos$`ALBSF:ALBSF_emb_WPC`


saveRDS(genesPMBSF_ALBSF_emb_WPC, "genesPMBSF_ALBSF_emb_WPC_01padj_0137lfc.rds")
saveRDS(genesALBSF_ALBSF_emb_WPC, "genesALBSF_ALBSF_emb_WPC_01padj_0137lfc.rds")

#PMBSFspc_ALBSFembWPC_ALBSFspc######################################################
#seleccion de datos
#cts
cts_PMBSF_ALBSF_emb_WPC<-cts[genesPMBSF_ALBSF_emb_WPC,]
cts_ALBSF_ALBSF_emb_WPC<-cts[genesALBSF_ALBSF_emb_WPC,]
dim(cts_ALBSF_ALBSF_emb_WPC)
dim(cts_PMBSF_ALBSF_emb_WPC)
cts_genes_especif<-rbind(cts_PMBSF_ALBSF_emb_WPC,cts_ALBSF_ALBSF_emb_WPC)
dim(cts_genes_especif)
#coldata
controlesALBSF<-cts_genes_especif[,1:4]
controlesPMBSF<-cts_genes_especif[,5:8]
ALBSF_EMB1<-cts_genes_especif[,9:12]
genes_especificos_data <- cbind(controlesALBSF,ALBSF_EMB1,controlesPMBSF)
controles_coldata_vs_ALBSF_EMB_coldata<-coldata[1:12,]
controlesALBSFcoldata<-coldata[1:4,]
controlesPMBSFcoldata<-coldata[5:8,]
ALBSF_EMB1coldata<-coldata[9:12,]
controles_coldata_vs_ALBSF_EMB_coldata <- rbind(controlesALBSFcoldata,ALBSF_EMB1coldata,controlesPMBSFcoldata)
dim(genes_especificos_data)
###dds###################################################################
dds_Spc_EmbWPC <- DESeqDataSetFromMatrix(countData = genes_especificos_data,
                                         colData = controles_coldata_vs_ALBSF_EMB_coldata,
                                         design = ~ Tratamiento + Region)
###prefiltering############################################################
keep_Spc_EmbWPC <- rowSums(counts(dds_Spc_EmbWPC)) >= 10
dds_Spc_EmbWPC <- dds_Spc_EmbWPC[keep_Spc_EmbWPC,]
##PCA######################################################################
rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC)
data_Spc_EmbWPC<-plotPCA(rld_Spc_EmbWPC, intgroup = c("Tratamiento","Region"),returnData=TRUE)
percentVar <- round(100 * attr(data_Spc_EmbWPC, "percentVar"))
p<-qplot(PC1, PC2, color=Tratamiento, shape=Region, data=data_Spc_EmbWPC, label=colnames(rld_Spc_EmbWPC)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "./figuras/PCA_SPCvsALBSF_EMB_WPC.png",
    width=500, height = 500)
p + geom_text(aes(label=colnames(rld_Spc_EmbWPC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                              alpha = 0.2, 
                                                                                              
                                                                                              show.legend = FALSE, 
                                                                                              
                                                                                              level = 0.75
) + theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(legend.title = element_text(size=15))+
  theme(legend.text = element_text(size=12))+
  guides(color = guide_legend(override.aes = list(size = 10)),shape = guide_legend(override.aes = list(size = 10)))+
  theme(plot.background = element_rect(fill = "white"))+
  theme(panel.background = element_rect(fill = 'white'))+
  geom_point(size = 6)
dev.off()

##PCA PC1 y PC3#####################################
rv2 <- rowVars(assay(rld_Spc_EmbWPC))

# select the ntop genes by variance
# select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
select2 <- order(rv2, decreasing=TRUE)[seq_len(min(500, length(rv2)))]

# perform a PCA on the data in assay(x) for the selected genes
pca_SPC_ALBSFembWPC <- prcomp(t(assay(rld_Spc_EmbWPC)[select2,]))
summary(pca_SPC_ALBSFembWPC)
# svg('PCA_elbow.svg', width = 10, height = 8)
plot(pca_SPC_ALBSFembWPC, type="lines")
# dev.off()
screeplot(pca_SPC_ALBSFembWPC)
# biplot(pca)

# the contribution to the total variance for each component
percentVar2 <- round(pca_SPC_ALBSFembWPC$sdev^2 / sum( pca_SPC_ALBSFembWPC$sdev^2 )*100)

# if (!all("Nucleus" %in% names(colData(rld)))) {
#   stop("the argument 'intgroup' should specify columns of colData(dds)")
# }

intgroup.df <- as.data.frame(colData(rld_Spc_EmbWPC)[, c("Tratamiento","Region"), drop=FALSE])

# add the intgroup factors together to create a new grouping factor
group <- if (length(c("Tratamiento","Region")) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(rld_Spc_EmbWPC)[["Region"]]
}
# pca$x
df_pca <- data.frame(PC1=pca_SPC_ALBSFembWPC$x[,1], PC2=pca_SPC_ALBSFembWPC$x[,2], PC3=pca_SPC_ALBSFembWPC$x[,3], 
                 group=group, intgroup.df, name=colnames(rld_Spc_EmbWPC))
df_pca
p<-qplot(PC1, PC3, color=Tratamiento, shape=Region, data=df_pca, label=colnames(rld_Spc_EmbWPC)) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar2[3],"% variance"))
png(filename = "figuras/PCA_PC1yPC3_SPC_ALBSembWPC_L2FC0137_01padj.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rld_Spc_EmbWPC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                              alpha = 0.2, 
                                                                                              
                                                                                              show.legend = FALSE, 
                                                                                              
                                                                                              level = 0.75
) 
dev.off()
##PCA PC2 y PC3##################################### 
p<-qplot(PC2, PC3, color=Tratamiento, shape=Region, data=df_pca, label=colnames(rld_Spc_EmbWPC)) +
  xlab(paste0("PC2: ",percentVar2[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar2[3],"% variance"))
png(filename = "figuras/PCA_PC2yPC3_SPC_ALBSembWPC_L2FC0137_01padj.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rld_Spc_EmbWPC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                              alpha = 0.2, 
                                                                                              
                                                                                              show.legend = FALSE, 
                                                                                              
                                                                                              level = 0.75
) 
dev.off()

##3Dplot#######################################################################
# https://github.com/mikelove/DESeq2/blob/master/R/plots.R
# plotPCA.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=FALSE)
# {
# calculate the variance for each gene
library("plot3D")
library("plot3Drgl")
rv <- rowVars(assay(rld_Spc_EmbWPC))

# select the ntop genes by variance
# select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld_Spc_EmbWPC)[select,]))
summary(pca)
# svg('PCA_elbow.svg', width = 10, height = 8)
plot(pca, type="lines")
# dev.off()
screeplot(pca)
# biplot(pca)

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# if (!all("Nucleus" %in% names(colData(rld)))) {
#   stop("the argument 'intgroup' should specify columns of colData(dds)")
# }
as.data.frame(colData(rld_Spc_EmbWPC))
intgroup.df <- as.data.frame(colData(rld_Spc_EmbWPC)[, c("Tratamiento","Region"), drop=FALSE])

# add the intgroup factors together to create a new grouping factor
group <- if (length(c("Tratamiento","Region")) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(rld_Spc_EmbWPC)[["Region"]]
}
# pca$x
d2 <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], 
                 group=group, intgroup.df, name=colnames(rld_Spc_EmbWPC))
d2
pdf('./figuras/PCA3D_top500rldBT_Tissue3label2_sin_texto_vista2.pdf', width = 6, height = 6)
scatter3D( d2$PC1, d2$PC2, d2$PC3,pch = 19, cex = 2,
           bty = "f", #"b2" "f"
           phi = 210, 
           theta =-27,
           main = NULL, 
           xlab = paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
           ylab = paste0("PC2: ",round(percentVar[2] * 100),"% variance"), 
           zlab = paste0("PC3: ",round(percentVar[3] * 100),"% variance"),
           # col = c("grey65", "#198ff7"),
           # colkey = list(at = c(2, 3, 4), side = 1),
           colkey = F,
           colvar = as.integer(d2$group),
           col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
           ticktype = "detailed")
text3D(d2$PC1, d2$PC2, d2$PC3,  labels = rownames(d2),
       add = TRUE, colkey = FALSE, cex = 1 )
dev.off()

options(rgl.printRglwidget = TRUE)
plotrgl()
rglwidget()
png(filename = "./figuras/PCA3D_top500rldBT_Tissue3label_2.png", width = 480, height = 480)
scatter3Drgl( d2$PC1, d2$PC2, d2$PC3,pch = 19, cex = 3,
              bty = "f", #"b2" "f"
              phi = -210, theta = -27,
              main = NULL, 
              xlab = paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
              ylab = paste0("PC2: ",round(percentVar[2] * 100),"% variance"), 
              zlab = paste0("PC3: ",round(percentVar[3] * 100),"% variance"),
              # col = c("grey65", "#198ff7"),
              # colkey = list(at = c(2, 3, 4), side = 1),
              colkey = FALSE,
              colvar = as.integer(d2$group),
              col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
              ticktype = "detailed")
dev.off()
text3Drgl(d2$PC1, d2$PC2, d2$PC3,  labels = rownames(d2),
          add = TRUE, colkey = FALSE, cex = 1 )

play3d( spin3d( axis = c(0, 0, 1), rpm = 2), duration = 5 )

movie3d(spin3d(axis = c(0, 0, 1), rpm = 2), duration = 7,
        # dir = getwd(), 
        movie = "pca3D1", type = "gif", dir = "./figuras/")

##Heatmap##########################################################
rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC, blind=FALSE)
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[rownames(genes_especificos_data), rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
mat_Spc_EmbWPC.scaled <- t(apply(mat_Spc_EmbWPC, 1, scale)) #center and scale each column (Z-score) then transpose
mat_Spc_EmbWPC.scaled<-(mat_Spc_EmbWPC - rowMeans(mat_Spc_EmbWPC))/rowSds(mat_Spc_EmbWPC)
colnames(mat_Spc_EmbWPC.scaled) <- paste0(rld_Spc_EmbWPC$Tratamiento,"-",rld_Spc_EmbWPC$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld_Spc_EmbWPC$Tratamiento ]
heatmap.2(mat_Spc_EmbWPC.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="none", main= "Spcf_vs_ALBSFemb_WPC_197up_vs_139down")

##Anotacion######
valuescomplex <- rownames(genes_especificos_data)
#lista_atributos<-listAttributes(ensembl)
tabla_anotacionComplex<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                              filters = "ensembl_gene_id",
                              values = valuescomplex,
                              mart = mouse)
data_anotadaComplex <- merge(genes_especificos_data, tabla_anotacionComplex %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE, sort=FALSE)
data_anotadaComplex

##Complex Heatmap##############################################################################
ALBSFembWPC_SPC<-rownames_to_column(as.data.frame(cts_genes_especif),"ensembl_gene_id")
dfC_anot<-getBM(attributes =c("ensembl_gene_id","mgi_symbol","gene_biotype") ,
                filters = "ensembl_gene_id",
                values = rownames(dfC),
                mart = mouse)
dfc2<-rownames_to_column(dfC,var = "ensembl_gene_id")
dfC_biotype <- inner_join(dfc2,dfC_anot,"ensembl_gene_id")
dfC_filtrado <- dfC_biotype %>% filter(gene_biotype == "protein_coding")
dfc_ALBSFembWPC_SPC<-inner_join(dfC_filtrado,ALBSFembWPC_SPC[1],"ensembl_gene_id")

top_25_genes <- arrange(dfc_ALBSFembWPC_SPC, desc(log2FoldChange)) %>% head(25)
low_25_genes <- dfc_ALBSFembWPC_SPC[order(dfc_ALBSFembWPC_SPC$log2FoldChange), ]
low_25_genes <- head(low_25_genes, 25)
top_50_genes <- rbind(top_25_genes,low_25_genes)


rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC, blind=FALSE)
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[rownames(genes_especificos_data), rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
mat_Spc_EmbWPC.scaled <- t(apply(mat_Spc_EmbWPC, 1, scale)) #center and scale each column (Z-score) then transpose
mat_Spc_EmbWPC.scaled<-(mat_Spc_EmbWPC - rowMeans(mat_Spc_EmbWPC))/rowSds(mat_Spc_EmbWPC)


rownames(mat_Spc_EmbWPC.scaled)<- data_anotadaComplex$mgi_symbol 
colnames(mat_Spc_EmbWPC.scaled) <- paste0(rld_Spc_EmbWPC$Tratamiento,"-",rld_Spc_EmbWPC$Region)
genes_anotacion_m<-match(rownames(mat_Spc_EmbWPC.scaled), top_50_genes$mgi_symbol)

# Crear el vector de posiciones
posiciones_matriz <- which(!is.na(genes_anotacion_m))


label_genes <- vector("character", length(genes_anotacion_m))
for (i in 1:length(genes_anotacion_m)) {
  if (!is.na(genes_anotacion_m[i])) {
    label_genes[i] <- top_50_genes$mgi_symbol[genes_anotacion_m[i]]
  }
}
label_genes<-subset(label_genes, nzchar(label_genes))
print(label_genes)
anotaciones_h = rowAnnotation(foo = anno_mark(at = posiciones_matriz,
                                              labels =label_genes))


col_fun= colorRamp2(seq(-1.5,1.5, length.out=100),
                    (colorRampPalette(c("darkorchid1","black","yellow1"))(100)),
                    space="RGB")
pdf("figuras/Heatmap_ALBSF_EMB_VS_ALBSF_ctrl_VS_PMBSF_ctrl_01padj.pdf",
    width=10, height = 10)
Heatmap(mat_Spc_EmbWPC.scaled,
        col=col_fun,
        na_col = "black",
        column_title = "Genes específicos modificados por la cauterización 197_PMBSF 139_ALBSF ",
        column_title_gp = gpar(fontsize = 18),
        right_annotation = anotaciones_h,
        row_names_gp = gpar(fontsize = 12),
        show_row_names = FALSE,
        cluster_columns = FALSE)
dev.off()

png(filename="figuras/Heatmap_ALBSF_EMB_VS_ALBSF_ctrl_VS_PMBSF_ctrl_01padj.png",
    width=750, height = 750)
Heatmap(mat_Spc_EmbWPC.scaled,
        col=col_fun,
        na_col = "black",
        column_title = "Genes específicos modificados por la cauterización 197_PMBSF 139_ALBSF ",
        column_title_gp = gpar(fontsize = 18),
        right_annotation = anotaciones_h,
        row_names_gp = gpar(fontsize = 12),
        show_row_names = FALSE,
        cluster_columns = FALSE)
dev.off()

#PMBSF_EMB_VS_PMBSF_ctrl########################################################
PMBSF_ctrl<-cts[,5:8]
PMBSF_emb_WPC<-cts[,13:16]
PMBSF_cts<-cbind(PMBSF_ctrl,PMBSF_emb_WPC)

PMBSF_coldat<-coldata[5:8,]
PMBSF_ctrl_coldata<-coldata[13:16,]
PMBSF_coldata<-rbind(PMBSF_coldat,PMBSF_ctrl_coldata)
##dds#######################################################
ddsP <- DESeqDataSetFromMatrix(countData = PMBSF_cts,
                               colData = PMBSF_coldata,
                               design = ~ Tratamiento)

ddsP
ddsP<-estimateSizeFactors(ddsP)

keepP <- rowSums(counts(ddsP)) >= 10
ddsP <- ddsP[keepP,]
rldP <- rlog(ddsP)

##PCA###########################################################################
dataP<-plotPCA(rldP, intgroup = "Tratamiento",returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, data=dataP, label=colnames(rldP)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p + geom_text(aes(label=colnames(rldP)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75
) 
##Expresión diferencial#########################################################
ddsP$Tratamiento <- relevel(ddsP$Tratamiento, "Control")
ddsP <- DESeq(ddsP)
resP<-results(ddsP,alpha=0.1, lfcThreshold= 0)
resP_apeglm<- lfcShrink(ddsP,res=resP,coef = "Tratamiento_Cauterizado_vs_Control")#Tabla resultado
summary(resP_apeglm)

#Filtrar los datos
sigsP <- na.omit(resP_apeglm)
sigsP<- sigsP[sigsP$padj < 0.1,]
dfP <- as.data.frame(sigsP)
dim(dfP)
summary(sigsP)
#plotMA
plotMA(resP_apeglm, ylim=c(-5,5))

##Anotación con biomaRt########################################################
#library(biomaRt)
ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "104")
rownames_ordenadosP <- sort(rownames(dfP))
datos_ordenados_P <- dfP[rownames_ordenados, ]
valuesP <- rownames(datos_ordenados_P)
#lista_atributos<-listAttributes(ensembl)
tabla_anotacionP<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = valuesP,
                        mart = mouse)

# Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combinedP <- merge(dfP, tabla_anotacionP %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)
#eliminamos duplicados
df_finalP<- unique(df_combinedP)
rownames(df_finalP)<-df_finalP$Row.names
df_finalP$Row.names<-NULL

##Heatmap#######################################################################
#library("genefilter")
#filtramos los genes para quedarnos con los significativos
df.topP <- df_finalP[(abs(df_finalP$log2FoldChange) > 0.137),]
df.topP <- df.topP[order(df.topP$log2FoldChange, decreasing = TRUE),]
df.topP$diffexpressed <- "NO"
df.topP$diffexpressed[df.topP$log2FoldChange > 0.137 ] <- "UP"
df.topP$diffexpressed[df.topP$log2FoldChange < -0.137] <- "DOWN"
table(df.topP$diffexpressed)
#creamos la matriz
rlog_outP <- rlog(ddsP, blind=FALSE) #get normalized count data from dds object
matP<-assay(rlog_outP)[rownames(df.topP), rownames(PMBSF_coldata)] #sig genes x samples
rownames(matP)<- df.topP$mgi_symbol 
base_mean <- rowMeans(matP)
head(matP)
matP.scaled <- t(apply(matP, 1, scale)) #center and scale each column (Z-score) then transpose
matP.scaled
colnames(matP.scaled) <- paste0(rlog_outP$Tratamiento,"-",rlog_outP$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rlog_outP$Tratamiento ]
png(filename = "Heatmap-PMBSF_EMB_VS_PMBSF_ctrlL2FC0137_01padj.png",
    width=1000, height = 1000)
heatmap.2(matP.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row", main= "PMBSF_EMB_VS_PMBSF_ctrl_2476up_vs2746down")
dev.off()
##Enriquecimiento##########################################################
genes_to_testP <- rownames(df.topP)
GO_resultsP <- enrichGO(gene = genes_to_testP, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_resultsP)
fitP <- plot(barplot(GO_resultsP, showCategory = 15))
png(filename = "dotplot-PMBSF_EMB_VS_PMBSF_ctrlL2FC0137_01padj.png",
    width=1000, height = 1000)
fitP2<- plot(dotplot(GO_resultsP, showCategory = 15))
dev.off()
##Volcano############################################
df_finalP$diffexpressed <- "NO"
df_finalP$diffexpressed[df_finalP$log2FoldChange > 0.137 ] <- "UP"
df_finalP$diffexpressed[df_finalP$log2FoldChange < -0.137] <- "DOWN"

#etiquetamos los genes diferenciados
df_finalP$delabel <- NA
df_finalP$delabel[df_finalP$diffexpressed != "NO"] <- df_finalP$mgi_symbol[df_finalP$diffexpressed != "NO"]
#representamos con ggplot
#library(ggrepel)
png(filename = "VOLCANO-PMBSF_EMB_VS_PMBSF_ctrlL2FC0137_01padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalP, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.137, 0.137), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()
#Diagrama de Venn#####################################################
PMBSF_mbWPC <- resP_apeglm[which(resP_apeglm$padj <= 0.1 & 
                                   abs(resP_apeglm$log2FoldChange) > 0.137 ),]
nrow(PMBSF_mbWPC)

venPMBSF_emb_WPC <-list(ALBSF=rownames(ALBSF_SpC), 
                        PMBSF=rownames(PMBSF_SpC), 
                        PMBSF_embWPC = rownames(PMBSF_mbWPC))
#library(ggplot2)
png(filename = "diag_vennPMBSF_emb_WPC_L2FC0137_01padj.png")
diag_vennPMBSF_emb_WPC=ggVennDiagram(venPMBSF_emb_WPC) + scale_fill_gradient(low="blue",high = "red")
dev.off()
png(filename = "diag_vennPMBSF_emb_WPC_L2FC0137_01padj.png")
diag_vennPMBSF_emb_WPC
dev.off()
##Fisher's Test##################################################
go.objALBSF_SpC2 <- newGeneOverlap(rownames(ALBSF_SpC),
                                  rownames(PMBSF_mbWPC),
                                  nrow(cts) )
go.objALBSF_SpC2 <- testGeneOverlap(go.objALBSF_SpC2)
print(go.objALBSF_SpC2)
go.objPMBSF_SpC2<- newGeneOverlap(rownames(PMBSF_SpC),
                                 rownames(PMBSF_mbWPC),
                                 nrow(cts) )
go.objPMBSF_SpC2 <- testGeneOverlap(go.objPMBSF_SpC2)
print(go.objPMBSF_SpC2)
##Listas de genes#########################################################
ven_intersect2 <- venn(venPMBSF_emb_WPC)
genes_especificosPMBSF_emb_WPC<-attr(ven_intersect2,"intersections")
genesPMBSF_PMBSF_emb_WPC<-genes_especificosPMBSF_emb_WPC$`PMBSF:PMBSF_embWPC`
genesALBSF_PMBSF_emb_WPC<-genes_especificosPMBSF_emb_WPC$`ALBSF:PMBSF_embWPC`

#GENES ESPECIFICOS PMBSF_embWPC################################################
cts_PMBSF_PMBSF_emb_WPC<-cts[genesPMBSF_PMBSF_emb_WPC,]
cts_ALBSF_PMBSF_emb_WPC<-cts[genesALBSF_PMBSF_emb_WPC,]
dim(cts_ALBSF_PMBSF_emb_WPC)
dim(cts_PMBSF_PMBSF_emb_WPC)
cts_genes_especif_PMBSF_embWPC<-rbind(cts_PMBSF_PMBSF_emb_WPC,cts_ALBSF_PMBSF_emb_WPC)
dim(cts_genes_especif_PMBSF_embWPC)
controles2<-cts_genes_especif[,1:8]
PMBSF_EMB1<-cts_genes_especif[,13:16]
genes_especificos_data2 <- cbind(controles2,PMBSF_EMB1)
controles_coldata_vs_PMBSF_EMB_coldata<-coldata[c(1:8,13:16),]
dim(genes_especificos_data2)


###dds###################################################################
dds_Spc_PMBSF_EmbWPC <- DESeqDataSetFromMatrix(countData = genes_especificos_data2,
                                               colData = controles_coldata_vs_PMBSF_EMB_coldata,
                                               design = ~ Tratamiento + Region)
###prefiltering############################################################
keep_Spc_PMBSF_EmbWPC <- rowSums(counts(dds_Spc_PMBSF_EmbWPC)) >= 10
dds_Spc_PMBSF_EmbWPC <- dds_Spc_PMBSF_EmbWPC[keep_Spc_PMBSF_EmbWPC,]
##PCA######################################################################
rld_Spc_PMBSF_EmbWPC<- rlog(dds_Spc_PMBSF_EmbWPC)
data_Spc_PMBSF_EmbWPC<-plotPCA(rld_Spc_PMBSF_EmbWPC, intgroup = c("Tratamiento","Region"),returnData=TRUE)
percentVar <- round(100 * attr(data_Spc_PMBSF_EmbWPC, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, shape=Region, data=data_Spc_PMBSF_EmbWPC, label=colnames(rld_Spc_PMBSF_EmbWPC)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p + geom_text(aes(label=colnames(rld_Spc_PMBSF_EmbWPC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                                    alpha = 0.2, 
                                                                                                    
                                                                                                    show.legend = FALSE, 
                                                                                                    
                                                                                                    level = 0.75
) 
dds_Spc_PMBSF_EmbWPC <- DESeq(dds_Spc_PMBSF_EmbWPC)
resSpc_PMBSF_EmbWPC<-results(dds_Spc_PMBSF_EmbWPC,alpha = 0.1, lfcThreshold= 0)
summary(resSpc_PMBSF_EmbWPC)
resSpc_PMBSF_EmbWPC_apeglm<- lfcShrink(dds_Spc_PMBSF_EmbWPC,res = resSpc_PMBSF_EmbWPC, coef = "Region_PMBSF_vs_ALBSF")#Tabla resultado
summary(resSpc_PMBSF_EmbWPC_apeglm)

##HEATMAP########################
rld_Spc_PMBSF_EmbWPC<- rlog(dds_Spc_PMBSF_EmbWPC, blind=FALSE)
mat_Spc_PMBSF_EmbWPC<-assay(rld_Spc_PMBSF_EmbWPC)[rownames(genes_especificos_data2), rownames(controles_coldata_vs_PMBSF_EMB_coldata)]
mat_Spc_PMBSF_EmbWPC.scaled <- t(apply(mat_Spc_PMBSF_EmbWPC, 1, scale)) #center and scale each column (Z-score) then transpose
mat_Spc_PMBSF_EmbWPC.scaled<-(mat_Spc_PMBSF_EmbWPC - rowMeans(mat_Spc_PMBSF_EmbWPC))/rowSds(mat_Spc_PMBSF_EmbWPC)
colnames(mat_Spc_PMBSF_EmbWPC.scaled) <- paste0(rld_Spc_PMBSF_EmbWPC$Tratamiento,"-",rld_Spc_PMBSF_EmbWPC$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld_Spc_PMBSF_EmbWPC$Tratamiento ]
png(filename = "Heatmap-Spc_genes_mod_PMBSF_emb_WPC_L2FC0137_01padj.png",
    width=1000, height = 1000)
heatmap.2(mat_Spc_PMBSF_EmbWPC.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="none", main = "Spcf_vs_PMBSF_emb_WPC_150up_179down")
dev.off()

#Comparacion P0 vs P8 General ALBSFemb vs Spc###################################################
##Datos P0 Lorenzo#################################################################
#Genes específicos de las regiones#
#resTissue1 <- lapply(openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx"),
 #                    openxlsx::read.xlsx,xlsxFile="./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx")
#names(resTissue1) <- openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx")
#head(resTissue1$ShApe_padj0.05_LFC0.322)


#PMBSF_SpC_P0_lor<- resTissue1$ShApe_Fulltable[which(resTissue1$ShApe_Fulltable$padj <= 0.1 &
                                  #                     resTissue1$ShApe_Fulltable$log2FoldChange > 0.137 ),]
#nrow(PMBSF_SpC_P0_lor)
# 377
#ALBSF_SpC_P0_lor <- resTissue1$ShApe_Fulltable[which(resTissue1$ShApe_Fulltable$padj <= 0.1 & 
                            #                           resTissue1$ShApe_Fulltable$log2FoldChange < -0.137 ),]
#nrow(ALBSF_SpC_P0_lor)

#Genes modificados por la cauterización#
#ALBSF_embWPCvsCtrl_P0_lor <- lapply(openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx"),
 #                    openxlsx::read.xlsx,xlsxFile="./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx")
#names(ALBSF_embWPCvsCtrl_P0_lor) <- openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx")
#head(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable)

#ALBSF_embWPC_P0_lor <- ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable[which(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable$padj <= 0.1 & 
 #                                  abs(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable$log2FoldChange) > 0.137 ),]
#nrow(ALBSF_embWPC_P0_lor)


##Mis datos##########
###Genes expresados en las 2 edades#######
genes_overlap_PO<-readRDS(file = "P0/P0_analisis/overlap_genes_STAR.rds")
genes_overlap_P8<-c(genesALBSF_ALBSF_emb_WPC,genesPMBSF_ALBSF_emb_WPC)
ven_overlap_P0_P8 <-list(overlapP0=rownames(genes_overlap_PO),
                         overlapP8=genes_overlap_P8)
diag_venn_overlap_P0_P8<-ggVennDiagram(ven_overlap_P0_P8,label_size = 16, label_alpha = 0.5) +
  theme(text = element_text(size = 16)) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P0vsP8") 

png(filename = "./figuras/diagVenn-overlap_genes_P0vsP8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_venn_overlap_P0_P8
dev.off()
ven_intersect_P0_P8 <- venn(ven_overlap_P0_P8)
genes_comunes_P0_P8 <-attr(ven_intersect_P0_P8,"intersections")
genes_comunes_P0_P8<-genes_comunes_P0_P8$`overlapP0:overlapP8`

###Anotacion
anotacion_genes_comunes_P0_P8<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = genes_comunes_P0_P8,
                        mart = mouse)
write_csv(anotacion_genes_comunes_P0_P8,"./figuras/genes_comunes_overlap_p0_p8.csv")
##PMBSF:ALBSF_emb Genes expresados en las dos edades################
genes_overlap_PMBSF_PO<-readRDS(file = "P0/P0_analisis/overlap_genes_PMBSF:ALBSF_mbWPC_STAR.rds")
genes_overlap_PMBSF_P8<-genesPMBSF_ALBSF_emb_WPC
ven_overlap_PMBSF_P0_P8<-list(overlapP0=genes_overlap_PMBSF_PO,
                              overlapP8=genes_overlap_PMBSF_P8)

diag_ven_overlap_PMBSF_P0_P8<-ggVennDiagram(ven_overlap_PMBSF_P0_P8,label_size = 16, label_alpha = 0.5) + 
                              scale_fill_gradient(low="blue",high = "red") + 
                              ggtitle("Overlapgenes P0vsP8 PMBSFspc:ALBSF_embWPC") + 
                              theme(text = element_text(size = 16))
png(filename = "./figuras/diagVenn-overlap_genes_PMBSF_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_PMBSF_P0_P8
dev.off()

ven_intersect_PMBSF_P0_P8 <- venn(ven_overlap_PMBSF_P0_P8)
genes_comunes_PMBSF_P0_P8 <-attr(ven_intersect_PMBSF_P0_P8,"intersections")
genes_comunes_PMBSF_P0_P8<-genes_comunes_PMBSF_P0_P8$`overlapP0:overlapP8`
###Anotacion
anotacion_genes_comunes_PMBSF_P0_P8<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                                     filters = "ensembl_gene_id",
                                     values = genes_comunes_PMBSF_P0_P8,
                                     mart = mouse)

##ALBSF:ALBSF_emb Genes expresados en las dos edades################
genes_overlap_ALBSF_PO<-readRDS(file = "P0/P0_analisis/overlap_genes_ALBSF:ALBSF_mbWPC_STAR.rds")
genes_overlap_ALBSF_P8<-genesALBSF_ALBSF_emb_WPC
ven_overlap_ALBSF_P0_P8<-list(overlapP0=genes_overlap_ALBSF_PO,
                              overlapP8=genes_overlap_ALBSF_P8)

diag_ven_overlap_ALBSF_P0_P8<-ggVennDiagram(ven_overlap_ALBSF_P0_P8,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P0vsP8 ALBSFspc:ALBSF_embWPC") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/diagVenn-overlap_genes_ALBSF_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_ALBSF_P0_P8
dev.off()

ven_intersect_ALBSF_P0_P8 <- venn(ven_overlap_ALBSF_P0_P8)
genes_comunes_ALBSF_P0_P8 <-attr(ven_intersect_ALBSF_P0_P8,"intersections")
genes_comunes_ALBSF_P0_P8<-genes_comunes_ALBSF_P0_P8$`overlapP0:overlapP8`
###Anotacion
anotacion_genes_comunes_ALBSF_P0_P8<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                                     filters = "ensembl_gene_id",
                                     values = genes_comunes_ALBSF_P0_P8,
                                     mart = mouse)


#Comparacion P0 vs P8 Especificos###################################################
##PMBSFspc Genes expresados en las dos edades################
genes_PMBSFspc_PO<-readRDS(file = "P0/P0_analisis/PMBSF_SpC_P0.rds")
genes_PMBSFspc_P8<-rownames(PMBSF_SpC)
ven_overlap_PMBSFspc_P0_P8<-list(P0=genes_PMBSFspc_PO,
                                 P8=genes_PMBSFspc_P8)

diag_ven_overlap_PMBSFspc_P0_P8<-ggVennDiagram(ven_overlap_PMBSFspc_P0_P8,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P0vsP8 PMBSFspc") + 
  theme(text = element_text(size = 16))

png(filename = "./figuras/diagVenn-overlap_genes_PMBSFspc_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_PMBSFspc_P0_P8
dev.off()

pdf("./figuras/diagVenn-overlap_genes_PMBSFspc_P0_P8_lfc0137_01padj.pdf",
    width=15, height = 15)
diag_ven_overlap_PMBSFspc_P0_P8
dev.off()

pdf("./figuras/diagVenn-overlap_genes_PMBSFspc_P0_P8_lfc0137_01_euler_padj.pdf",
    width=10, height = 10)
plot(euler(ven_overlap_PMBSFspc_P0_P8, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 10,
     legend = F)
dev.off()

#metemos el nombre de los genes en una lista
ven_intersect_PMBSFspc_P0_P8 <- venn(ven_overlap_PMBSFspc_P0_P8)
genes_comunes_PMBSFspc_P0_P8 <-attr(ven_intersect_PMBSFspc_P0_P8,"intersections")
genes_comunes_PMBSFspc_P0_P8<-genes_comunes_PMBSFspc_P0_P8$`overlapP0:overlapP8`

###Anotacion
anotacion_genes_comunes_PMBSFspc_P0_P8<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                                           filters = "ensembl_gene_id",
                                           values = genes_comunes_PMBSFspc_P0_P8,
                                           mart = mouse)

##ALBSFspc Genes expresados en las dos edades################
genes_ALBSFspc_PO<-readRDS(file = "P0/P0_analisis/ALBSF_SpC_P0.rds")
genes_ALBSFspc_P8<-rownames(ALBSF_SpC)
ven_overlap_ALBSFspc_P0_P8<-list(P0=genes_ALBSFspc_PO,
                                 P8=genes_ALBSFspc_P8)

diag_ven_overlap_ALBSFspc_P0_P8<-ggVennDiagram(ven_overlap_ALBSFspc_P0_P8,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P0vsP8 ALBSFspc") + 
  theme(text = element_text(size = 16))

png(filename = "./figuras/diagVenn-overlap_genes_ALBSFspc_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_ALBSFspc_P0_P8
dev.off()

pdf("./figuras/diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P8_lfc0137_01padj_euler.pdf",
    width=1000, height = 1000)
diag_ven_overlap_ALBSFspc_P0_P8
dev.off()

pdf("./figuras/diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P8_lfc0137_01padj_euler.pdf",
    width=10, height = 10)
plot(euler(ven_overlap_ALBSFspc_P0_P8, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 10,
     legend = F)
dev.off()

ven_intersect_ALBSFspc_P0_P8 <- venn(ven_overlap_ALBSFspc_P0_P8)
genes_comunes_ALBSFspc_P0_P8 <-attr(ven_intersect_ALBSFspc_P0_P8,"intersections")
genes_comunes_ALBSFspc_P0_P8<-genes_comunes_ALBSFspc_P0_P8$`overlapP0:overlapP8`

###Anotacion
anotacion_genes_comunes_ALBSFspc_P0_P8<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                                              filters = "ensembl_gene_id",
                                              values = genes_comunes_ALBSFspc_P0_P8,
                                              mart = mouse)

##Juntos###############################
ven_overlap_SPC<-list(P0_PMBSF=genes_PMBSFspc_PO,
                                 P8_PMBSF=genes_PMBSFspc_P8,
                                 P0_ALBSF=genes_ALBSFspc_PO,
                                 P8_ALBSF=genes_ALBSFspc_P8)
pdf("./figuras/diag_ven_SPC_P8yP0_euler.pdf",
    width=10, height = 10)
plot(euler(ven_overlap_SPC, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 10,
     legend = F)
dev.off()

#DiagVenn SPC_P8yP0 y los genes ALBSFembWPC P0############################

ALBSFemb_WPC_P0_genes<-readRDS("./P0/P0_analisis/ALBSF_embWPC_P0.rds")
ven_overlap_SPC_P8yP0_ALBSFembWPC_P0<-list(ALBSFspc_Both=genes_comunes_ALBSFspc_P0_P8,
                                           PMBSFspc_Both=genes_comunes_PMBSFspc_P0_P8,
                                           ALBSFemb_WPC_P0=ALBSFemb_WPC_P0_genes)

diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P0<-ggVennDiagram(ven_overlap_SPC_P8yP0_ALBSFembWPC_P0,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes específicos coincidentes en P0 y P8 expresados en ALBSF_emb_WPC a P0") + 
  theme(text = element_text(size = 16))

png(filename = "./figuras/diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P0_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P0
dev.off()

#DiagVenn SPC_P8yP0 y los genes ALBSFembWPC P8############################
ALBSFemb_WPC_P8_genes<-rownames(ALBSF_mbWPC)
ven_overlap_SPC_P8yP0_ALBSFembWPC_P8<-list(ALBSFspc_Both=genes_comunes_ALBSFspc_P0_P8,
                                           PMBSFspc_Both=genes_comunes_PMBSFspc_P0_P8,
                                           ALBSFemb_WPC_P8=ALBSFemb_WPC_P8_genes)

diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P8<-ggVennDiagram(ven_overlap_SPC_P8yP0_ALBSFembWPC_P8,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes específicos coincidentes en P0 y P8 expresados en ALBSF_emb_WPC a P8") + 
  theme(text = element_text(size = 16))


pdf("./figuras/diag_ven_overlap_SPC_P8yP0_ALBSFembWPC_P8_lfc0137_euler_padj01.pdf",
    width=10, height = 10)
plot(euler(ven_overlap_SPC_P8yP0_ALBSFembWPC_P8, shape = "circle"), 
     # fills = c("#196af7", "#34b31b", "#f71951"), edges = T,
     quantities = list(type = c( "counts")), #, "percent"
     fontsize = 10,
     legend = F)
dev.off()

#DiagVenn ALBSFspc P0vs08#################################################
ALBSFemb_WPC_P8_genes<-rownames(ALBSF_mbWPC)
ALBSFemb_WPC_P0_genes<- readRDS("P0/P0_analisis/ALBSF_embWPC_P0.rds")
ven_overlap_SPC_P8yP0_ALBSFembWPC_P8<-list(ALBSFemb_WPC_P0=ALBSFemb_WPC_P0_genes,
                                           ALBSFemb_WPC_P8=ALBSFemb_WPC_P8_genes)

diag_ven_overlap_P8yP0_ALBSFembWPC<-ggVennDiagram(ven_overlap_SPC_P8yP0_ALBSFembWPC_P8,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes de ALBSF_emb_WPC expresados a P0 y P8") + 
  theme(text = element_text(size = 16))

png(filename = "./figuras/diag_ven_ALBSFembWPC_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_P8yP0_ALBSFembWPC
dev.off()

##Enriquecimiento###########################

###Overlapping P8 genes##########################################################
emb_WPC_LFC<-subset(resA_apeglm,rownames(resA_apeglm) %in% rownames(cts_genes_especif) )
emb_WPC_LFC<-rownames_to_column(as.data.frame(emb_WPC_LFC), var="ensembl_gene_id")
emb_WPC_LFC<- inner_join(emb_WPC_LFC,dfC_anot,"ensembl_gene_id") 
GO_resultsgenes_overlapP8_2<- enrichGO(gene = emb_WPC_LFC$ensembl_gene_id, 
                                   OrgDb = "org.Mm.eg.db", 
                                   keyType = "ENSEMBL", 
                                   ont = "BP",
                                   readable=TRUE,
                                   pvalueCutoff = 0.2,
                                   pAdjustMethod = "BH")

#procesos de interés
go_ids<- c(	"GO:0050432","GO:0051932","GO:0060828","GO:0000910","GO:0016079","GO:0006887","GO:1901632","GO:0044843","GO:0050808","GO:0071634","GO:0034349","GO:0006836","GO:0051937","GO:0098742")
filtered_go <- GO_resultsgenes_overlapP8_2 %>%
  dplyr::filter(ID %in% go_ids)

png(filename = "figuras/dotplot_genes_spc_mod_P8.png", width = 600, height = 600)
dotplot(filtered_go,showCategory = 14, font.size=19)
dev.off()

genes_emb_WPC_LFC<-emb_WPC_LFC$log2FoldChange
names(genes_emb_WPC_LFC) <- emb_WPC_LFC$mgi_symbol
png(filename = "figuras/cnetplot_genes_spc_mod_P8.png", width = 800, height = 800)
cnetplot(filtered_go,color.params = list(foldChange=genes_emb_WPC_LFC),showCategory = 14, categorySize="pvalue",fontsize=19)+scale_color_gradient2(name='log2FoldChange', low='darkgreen', high='firebrick')
dev.off()

#emriquecimiento por región
compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genes<-list(genesPMBSF_ALBSF_emb_WPC,genesALBSF_ALBSF_emb_WPC)
names(compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genes) <- c("PMBSF:ALBSF_embWPC", "ALBSF:ALBSF_embWPC")
ck_PMBSF_VS_ALBSF_genes<- compareCluster(geneCluster = compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genes,
                                         fun = enrichGO,
                                         OrgDb = "org.Mm.eg.db", 
                                         keyType = "ENSEMBL", 
                                         ont = "BP",
                                         readable=TRUE,
                                         pvalueCutoff = 0.2,
                                         pAdjustMethod = "BH")
filtered_go_clster<-ck_PMBSF_VS_ALBSF_genes %>%
  dplyr::filter(ID %in% go_ids)

dotplot(filtered_go_clster)
write_csv(ck_PMBSF_VS_ALBSF_genes@compareClusterResult, "./figuras/GO_results_overlap_genes.csv")
overlap_P8_PMBSF<-data.frame(ensembl_gene_id=compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genes$`PMBSF:ALBSF_embWPC`,Region="PMBSF:ALBSF_embWPC")
overlap_P8_ALBSF<-data.frame(ensembl_gene_id=compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genes$`ALBSF:ALBSF_embWPC`,Region="ALBSF:ALBSF_embWPC")

overlap_P8<-rbind(overlap_P8_PMBSF,overlap_P8_ALBSF)
overlapP8_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                                              filters = "ensembl_gene_id",
                                              values = overlap_P8$ensembl_gene_id,
                                              mart = mouse)

overlapP8_anotacion<-inner_join(overlap_P8,overlapP8_anotacion,"ensembl_gene_id")
write_csv(overlapP8_anotacion,"./figuras/genes_overlapP8.csv")

overlapP8_P0_genes_anotados<-inner_join(overlapP8_anotacion,anotacion_genes_comunes_P0_P8,"ensembl_gene_id")
write_csv(overlapP8_P0_genes_anotados,"./figuras/genes_overlapP0&P8.csv")

# Filtrar valores etiquetados con "PMBSF:ALBSF_embWPC"
PMBSF_ALBSF_embWPC_values <- subset(overlapP8_P0_genes_anotados, Region == "PMBSF:ALBSF_embWPC")

# Filtrar valores etiquetados con "ALBSF:ALBSF_embWPC"
ALBSF_ALBSF_embWPC_values <- subset(overlapP8_P0_genes_anotados, Region == "ALBSF:ALBSF_embWPC")

png(filename = "./figuras/cnetplotBP_overlapping_genes_P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(ck_PMBSF_VS_ALBSF_genes, showCategory=categorys,cex.params=list( category_label = 2, gene_label = 2))
dev.off()

png(filename = "./figuras/dotplotBP_overlapping_genes_P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(ck_PMBSF_VS_ALBSF_genes ,showCategory = 15,by="geneRatio")
dev.off()

genes_overP8<-c(genesPMBSF_ALBSF_emb_WPC,genesALBSF_ALBSF_emb_WPC)

GO_results_PMBSF_p8_overlap <- enrichGO(gene = genesALBSF_ALBSF_emb_WPC, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
dotplot(GO_results_PMBSF_p8_overlap)

###Overlapping P0 genes################################################

compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genesP0<-list(genes_overlap_PMBSF_PO,genes_overlap_ALBSF_PO)
names(compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genesP0) <- c("PMBSF:ALBSF_embWPC", "ALBSF:ALBSF_embWPC")
ck_PMBSF_VS_ALBSF_genesP0<- compareCluster(geneCluster = compare_ALBSF_emb_WPCPMBSF_VS_ALBSF_emb_WPCALBSF_genesP0,
                                         fun = enrichGO,
                                         OrgDb = "org.Mm.eg.db", 
                                         keyType = "ENSEMBL", 
                                         ont = "BP",
                                         readable=TRUE,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH")

png(filename = "./figuras/cnetplotBP_overlapping_genes_P0_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(ck_PMBSF_VS_ALBSF_genesP0, showCatgory=10,cex.params=list( category_label = 2, gene_label = 2))
dev.off()

png(filename = "./figuras/dotplotBP_overlapping_genes_P0_L2FC0137_01padj.png",
    width=600, height = 600)
dotplot(ck_PMBSF_VS_ALBSF_genesP0,  showCategory = 20, font.size=19)
dev.off()
### Overlap P8vsP0#####
genes_comunes_P0_P8
GO_resultsgenes_overlapP8vsP0<- enrichGO(gene = genes_comunes_P0_P8, 
                                       OrgDb = "org.Mm.eg.db", 
                                       keyType = "ENSEMBL", 
                                       ont = "BP",
                                       readable=TRUE,
                                       pvalueCutoff = 0.2,
                                       pAdjustMethod = "BH")
barplot(GO_resultsgenes_overlapP8vsP0,showCategory=8)

###Spc P8 genes################################################
PMBSF_SpC_genes<-rownames(PMBSF_SpC)
ALBSF_SpC_genes<-rownames(ALBSF_SpC)

compare_spc_genesP8<-list(PMBSF_SpC_genes,ALBSF_SpC_genes)
names(compare_spc_genesP8) <- c("PMBSF", "ALBSF")
ck_spc_genesP8<- compareCluster(geneCluster = compare_spc_genesP8,
                                           fun = enrichGO,
                                           OrgDb = "org.Mm.eg.db", 
                                           keyType = "ENSEMBL", 
                                           ont = "BP",
                                           readable=TRUE,
                                           pvalueCutoff = 0.05,
                                           pAdjustMethod = "BH")

png(filename = "./figuras/cnetplotBP_SPC_genes_P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(ck_spc_genesP8, showCatgory=10,cex.params=list( category_label = 2, gene_label = 2))
dev.off()

png(filename = "./figuras/dotplotBP_SPC_genes_P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(ck_spc_genesP8,  showCategory = 10, by="geneRatio")
dev.off()

###spc P8vsP0 overlapping############################
GO_results_overlapP8_P0<- enrichGO(gene = compare_spc_genesP8yP0, 
                                       OrgDb = "org.Mm.eg.db", 
                                       keyType = "ENSEMBL", 
                                       ont = "BP",
                                       readable=TRUE,
                                       pvalueCutoff = 0.2,
                                       pAdjustMethod = "BH")
dotplot(GO_results_overlapP8_P0,showCategory=15)

compare_spc_genesP8yP0<-list(genes_comunes_PMBSFspc_P0_P8,genes_comunes_ALBSFspc_P0_P8)
names(compare_spc_genesP8yP0) <- c("PMBSF:P0&P8", "ALBSF:P0&P8")
ck_spc_genesP8yP0<- compareCluster(geneCluster = compare_spc_genesP8yP0,
                                fun = enrichGO,
                                OrgDb = "org.Mm.eg.db", 
                                keyType = "ENSEMBL", 
                                ont = "BP",
                                readable=TRUE,
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")

png(filename = "./figuras/cnetplotBP_SPC_genes_P0&P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(ck_spc_genesP8yP0, showCatgory=10,cex.params=list( category_label = 2, gene_label = 2))
dev.off()

png(filename = "./figuras/dotplotBP_SPC_genes_P0&P8_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(ck_spc_genesP8yP0,  showCategory = 10, by="geneRatio")
dev.off()
#Scatter plots#########################################################

##SPC:ALBSF_embWPC P8#########################################################
#Datos ALBSF_embWPC
ALBSF_embWPC<-as.data.frame(resA_apeglm)
ALBSF_embWPC<-rownames_to_column(ALBSF_embWPC,var="EnsemblID")

#Datos SPC
Tissue<-as.data.frame(resC_apglm)
Tissue<-rownames_to_column(Tissue,var="EnsemblID")

#tabla conjunta
ALBSF_embWPC_vs_Tissue<- full_join(Tissue,ALBSF_embWPC,"EnsemblID")
ALBSF_embWPC_vs_Tissue_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                        filters = "ensembl_gene_id",
                        values = ALBSF_embWPC_vs_Tissue$EnsemblID,
                        mart = mouse)
ALBSF_embWPC_vs_Tissue_anotacion$EnsemblID<-ALBSF_embWPC_vs_Tissue_anotacion$ensembl_gene_id
ALBSF_embWPC_vs_Tissue_anotacion$ensembl_gene_id<-NULL
ALBSF_embWPC_vs_Tissue_completa<- inner_join(ALBSF_embWPC_vs_Tissue,ALBSF_embWPC_vs_Tissue_anotacion,"EnsemblID")
#etiquetado
ALBSF_embWPC_vs_Tissue_completa <- ALBSF_embWPC_vs_Tissue_completa %>%
  mutate(label = ifelse(EnsemblID %in% genesPMBSF_ALBSF_emb_WPC, "PMBSF:ALBSF_emb",
                        ifelse(EnsemblID %in% genesALBSF_ALBSF_emb_WPC, "ALBSF:ALBSF_emb", "Neg")))
ALBSF_embWPC_vs_Tissue_completa <- replace(ALBSF_embWPC_vs_Tissue_completa, is.na(ALBSF_embWPC_vs_Tissue_completa), 0)

#numero de genes
table(ALBSF_embWPC_vs_Tissue_completa$label)
#grafica
fourway <-ggplot(data=ALBSF_embWPC_vs_Tissue_completa,
                 aes(x = log2FoldChange.x,
                     y = log2FoldChange.y,
                     color= label,
                     label=mgi_symbol,
                     text = paste("GeneSymbol:", mgi_symbol)))  +
  geom_point(data = ALBSF_embWPC_vs_Tissue_completa[ALBSF_embWPC_vs_Tissue_completa$label=="Neg",], aes(col=label), size=3, alpha = (0.1)) +
  geom_point(data = ALBSF_embWPC_vs_Tissue_completa[ALBSF_embWPC_vs_Tissue_completa$label!="Neg",], aes(col=label), size=4, alpha = (0.8)) +
  labs(x = "LFC Tissue1_PMBSF_vs_ALBSF ",
       y = "LFC ALBSF_EmbWPC_vs_Control_ALBSF",
       color = "PMBSF UP/ ALBSF DOWN")+
  #shape = "Log2 Fold Change",
  #size = "Fold Change") +
  scale_color_manual(values=c("#00a4a2", "grey","#b05d95")) +
  geom_vline(xintercept=0, col="black",linetype="dashed") +
  geom_hline(yintercept=0, col="black",linetype="dashed")+
  ylim(-3,3)+
  xlim(-3,3)+
  geom_text_repel(data=ALBSF_embWPC_vs_Tissue_completa[ALBSF_embWPC_vs_Tissue_completa$label!="Neg",],size=7) +
  theme_dark() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "ALBSF_embWPC:SPC P8 ")+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10))
png(filename="./figuras/scatterplot_ALBSF_embWPC:SPC_P8.png",    
    width=800, height = 800)
fourway
dev.off()

pdf("./figuras/scatterplot_ALBSF_embWPC:SPC_P8.pdf",    
    width=20, height = 20)
fourway
dev.off()
fourway_plotly <- ggplotly(fourway)

fourway_plotly
htmlwidgets::saveWidget(config(fourway_plotly, showLink = T), "./figuras/scatterplot_ALBSF_embWPC:SPC_P8.html")

##SPC P8vsP0 : ALBSF_embWPC_P0&P8####################################

#Datos spc_P8
spc_P8<-as.data.frame(resC_apglm)
spc_P8<-rownames_to_column(spc_P8,var="EnsemblID")

#Datos SPC_P0
spcP0<-readRDS(file = "P0/P0_analisis/SPC_resAPEGLM_P0.rds")
spcP0<-as.data.frame(spcP0)
spcP0<-rownames_to_column(spcP0,var="EnsemblID")

#tabla conjunta
spc_P8_vs_spcP0<- inner_join(spc_P8,spcP0,"EnsemblID")
spc_P8_vs_spcP0_anot<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                                        filters = "ensembl_gene_id",
                                        values = spc_P8_vs_spcP0$EnsemblID,
                                        mart = mouse)
spc_P8_vs_spcP0_anot$EnsemblID<-spc_P8_vs_spcP0_anot$ensembl_gene_id
spc_P8_vs_spcP0_anot$ensembl_gene_id<-NULL
spc_P8_vs_spcP0_completa<- inner_join(spc_P8_vs_spcP0,spc_P8_vs_spcP0_anot,"EnsemblID")
#etiquetado
spc_P8_vs_spcP0_completa <- spc_P8_vs_spcP0_completa %>%
  mutate(label = ifelse(EnsemblID %in% genes_comunes_PMBSFspc_P0_P8_lor , "PMBSF_spc",
                        ifelse(EnsemblID %in% genes_comunes_ALBSFspc_P0_P8_lor, "ALBSF_spc", "Neg")))

#numero de genes
table(spc_P8_vs_spcP0_completa$label)
#grafica
fourway <-ggplot(data=spc_P8_vs_spcP0_completa,
                 aes(x = log2FoldChange.x,
                     y = log2FoldChange.y,
                     color= label,
                     #label=mgi_symbol,
                     text = paste("GeneSymbol:", mgi_symbol)))  +
  geom_point(data = spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$label=="Neg",], aes(col=label), size=2, alpha = (0.1)) +
  geom_point(data = spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$label!="Neg",], aes(col=label), size=5, alpha = (0.9)) +
  labs(x = "LFC SPC P8 ",
       y = "LFC SPC P0",
       color = "PMBSF UP/ ALBSF DOWN")+
  #shape = "Log2 Fold Change",
  #size = "Fold Change") +
  scale_color_manual(values=c("red", "black","darkgreen"))+
  geom_vline(xintercept=0, col="black",linetype="dashed") +
  geom_hline(yintercept=0, col="black",linetype="dashed")+
  ylim(-3,3)+
  xlim(-3,3)+
  #geom_text_repel(data=spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$label!="Neg",]) +
  theme_dark() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Genes específicos de región conservados a P0&P8 ")+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10))
png(filename="./figuras/scatterplot_SPC_P8vsP0_sin_modEMBwpc.png",    
    width=1000, height = 1000)
fourway
dev.off()
pdf("./figuras/scatterplot_SPC_P8vsP0_sin_modEMBwpc.pdf",    
    width=10, height = 10)
fourway
dev.off()
fourway_plotly <- ggplotly(fourway)

fourway_plotly
htmlwidgets::saveWidget(config(fourway_plotly, showLink = T), "./figuras/scatterplot_SPC_P8vsP0_ALBSF_embWPC_P0&P8.html")

##SPC P8vsP0####################################

#Datos spc_P8
spc_P8<-as.data.frame(resC_apglm)
spc_P8<-rownames_to_column(spc_P8,var="EnsemblID")

#Datos SPC_P0
spcP0_lor<-resTissue1$ShApe_Fulltable
spcP0<-readRDS(file = "P0/P0_analisis/SPC_resAPEGLM_P0.rds")
spcP0<-as.data.frame(spcP0)
spcP0<-rownames_to_column(spcP0,var="EnsemblID")

#tabla conjunta
spc_P8_vs_spcP0_lor<- inner_join(spc_P8,spcP0,"EnsemblID")
spc_P8_vs_spcP0_lor_anot<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                                filters = "ensembl_gene_id",
                                values = spc_P8_vs_spcP0$EnsemblID,
                                mart = mouse)
spc_P8_vs_spcP0_lor_anot$EnsemblID<-spc_P8_vs_spcP0_lor_anot$ensembl_gene_id
spc_P8_vs_spcP0_lor_anot$ensembl_gene_id<-NULL
spc_P8_vs_spcP0_completa<- inner_join(spc_P8_vs_spcP0_lor,spc_P8_vs_spcP0_lor_anot,"EnsemblID")

#etiquetado
spc_P8_vs_spcP0_completa <- spc_P8_vs_spcP0_completa %>%
  mutate(label = ifelse(EnsemblID %in% genesSPC_ALBSF_emb_WPC_P0_P8 , "EmbWPC mod P8&P0",
                        ifelse(EnsemblID %in% genes_comunes_PMBSFspc_P0_P8_lor , "PMBSFspc P8&P0",
                        ifelse(EnsemblID %in% genes_comunes_ALBSFspc_P0_P8_lor, "ALBSFspc P8&P0", "Neg"))))

spc_P8_vs_spcP0_completa <- replace(spc_P8_vs_spcP0_completa, is.na(spc_P8_vs_spcP0_completa), 0)
spc_P8_vs_spcP0_completa <- spc_P8_vs_spcP0_completa %>%
  mutate(delabel = ifelse(EnsemblID %in% genes_comunes_PMBSF_P0_P8_lor , "PMBSFspc:ALBSFembWPC",
                          ifelse(EnsemblID %in% genes_comunes_ALBSF_P0_P8_lor, "ALBSFspc:ALBSFembWPC","Neg")))
#numero de genes
table(spc_P8_vs_spcP0_completa$label)
#grafica
fourway <-ggplot(data=spc_P8_vs_spcP0_completa,
                 aes(x = log2FoldChange.x,
                     y = log2FoldChange.y,
                     color= label,
                     label=mgi_symbol,
                     text = paste("mgi_symbol:", mgi_symbol)))  +
  geom_point(data = spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$label=="Neg",], aes(col=label), size=2, alpha = (0.1)) +
  geom_point(data = spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$label!="Neg",], aes(col=label), size=5, alpha = (1)) +
  geom_point(data = spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$delabel!="Neg",], aes(col=label), size=5, alpha = (1)) +
  labs(x = "LFC PMBSFctrl_vs_ALBSFctrl  P8 ",
       y = "LFC PMBSFctrl_vs_ALBSFctrl  P0",
       color = "PMBSF UP/ ALBSF DOWN")+
  #shape = "Log2 Fold Change",
  #size = "Fold Change") +
  scale_color_manual(values=c("red","blue", "black","darkgreen"))+
  geom_vline(xintercept=0, col="black",linetype="dashed") +
  geom_hline(yintercept=0, col="black",linetype="dashed")+
  ylim(-3,3)+
  xlim(-3,3)+
  geom_text_repel(data=spc_P8_vs_spcP0_completa[spc_P8_vs_spcP0_completa$delabel!="Neg",]) +
  theme_dark() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Genes especificos de region conservado en P0&P8",
       subtitle = "                                           Los genes etiquetados están modificados por la cauterización a P0&P8")+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10))
png(filename="./figuras/scatterplot_SPC_P8vsP0_mod_EMB_WPC.png",    
    width=1000, height = 1000)
fourway
dev.off()

pdf("./figuras/scatterplot_SPC_P8vsP0_mod_EMB_WPC.pdf",    
    width=8, height = 8)
fourway
dev.off()
fourway_plotly <- ggplotly(fourway)

fourway_plotly
htmlwidgets::saveWidget(config(fourway_plotly, showLink = T), "./figuras/scatterplot_SPC_P8vsP0.html")


#Comprobación datos P0 lore y mios##############################
X85_65genesoverlaping <- read_table("./85_65genesoverlaping.txt")
genesX85_65genesoverlaping<-X85_65genesoverlaping$EnsemblID

ven_comprobacion<-list(Mios=rownames(genes_overlap_PO),
                                Lorenzo=genesX85_65genesoverlaping)

diag_comprobacion<-ggVennDiagram(ven_comprobacion,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Comprobación si los genes de mi análisis P0 son iguales al de Lorenzo") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/diagVenn-comprobacion_analisis_P0.png",
    width=1000, height = 1000)
diag_comprobacion
dev.off()

#TPMS######################################
cts_length<-read.table("./matriz_cuentas/featurecountstotal_STAR.txt",header = TRUE)
rownames(cts_length) <- cts_length$Geneid
head(cts_length)

cts_length2<-cts_length[,-c(1,2,3,4,5,6)]
head(cts_length2)
length_k_2=as.data.frame(cbind(row.names=rownames(cts_length2), 
                               length_k=cts_length$Length[match(rownames(cts_length2),rownames(cts_length)) ]
))

head(length_k_2)
length_k_2<-length_k_2 %>% 
  mutate(length_k = as.numeric(length_k)/1000)

tpm<-function(counts,lengths) {
  rate<-counts/lengths
  rate/sum(rate)*1e6
}

df_tpm <- apply(cts_length2, 2, function(x) tpm(x,length_k_2$length_k))
head(df_tpm)
colSums(df_tpm)
colMeans(df_tpm)

##Boxplots de los tpm###################
goi<-rownames(genes_especificos_data)
coldata_tpm<-as.data.frame(coldata)
coldata_tpm$Names<- rownames(coldata_tpm)

dds_tpm <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata_tpm,
                                  design = ~ Region + Tratamiento)

deseqdataCPAEA<-dds_tpm[,!grepl(paste(c("EP"), collapse = "|"),
                                dds_tpm$Names)]
as.data.frame(colData(deseqdataCPAEA))

deseqdataTPM <- deseqdataCPAEA


#el logaritmo en base 2 de los tpm de los genes que me interesan
tpmscounts <- t(log2((df_tpm[match(goi, rownames(df_tpm)),])+1 )) %>%  #t(log2((RlogctsASTROSwoCXespecificgenes +.5) )) %>%
  merge(colData(deseqdataTPM), ., by="row.names") %>%
  tidyr::gather(geneEnsembl, expression, (ncol(.)-length(goi)+1):ncol(.))
head(tpmscounts)

tcountsannot <- getBM(filters = "ensembl_gene_id",
                      values = tpmscounts$geneEnsembl,
                      attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype"),
                      mart = mouse,
                      useCache = T)

m_cxKD1 <- match(tpmscounts$geneEnsembl, tcountsannot$ensembl_gene_id)
tpmscounts2 <- cbind(tpmscounts[,1:ncol(tpmscounts)], 
                     # EnsemblID=rownames(rld2ScBF),
                     GeneSymbol=tcountsannot$mgi_symbol[m_cxKD1])
head(tpmscounts2)
tpmscounts2$GeneSymbol <- factor(tpmscounts2$GeneSymbol, levels = unique(tpmscounts2$GeneSymbol))

# Generar la columna "Tissue2" con los valores deseados
tpmscounts2$Tissue2 <- ifelse(tpmscounts2$Region == "ALBSF" & tpmscounts2$Tratamiento == "Control", "C_ALBSF",
                              ifelse(tpmscounts2$Region == "ALBSF" & tpmscounts2$Tratamiento == "Cauterizado", "E_ALBSF",
                                     ifelse(tpmscounts2$Region == "PMBSF" & tpmscounts2$Tratamiento == "Control", "C_PMBSF", NA)))

table(tpmscounts2$Tissue2)
dir.create("./ExpressionFigures")
p1<-ggplot(tpmscounts2, aes(Tissue2, expression)) +
  geom_boxplot()+
  geom_jitter(height = 0, width = 0.1, size=1)+
  facet_wrap(~ GeneSymbol, scales="free_y") 
p1 +theme_classic()+ #base_size = 25
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_blank(), strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust=1, colour="black"),
        axis.text.y = element_text(angle = 0, hjust=1, colour="black"),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        plot.title=element_text(size=25, hjust=0.5, face="bold.italic", colour="black", vjust=-1),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        legend.position = "none")+
  ggsave("./ExpressionFigures/TPMs.pdf", 
         device="pdf", dpi = 600, 
         width = 20, height = 10.5, units = "in")



plot_list <- list()
for (i in goi){
  print(i)
  plot_list[[i]]<- ggplot(tpmscounts2[!grepl("EP", tpmscounts2$Name),] %>%
                            # dplyr::filter(tpmscounts2, geneEnsembl == i )
                            dplyr::filter( geneEnsembl == i ),
                          aes(Tissue2, expression)) + #aqui van los valores de cada whisker para cada gene
    stat_boxplot( #aes(Tissue2, expression) ,
      geom='errorbar', linetype=1, width=0.5)+  #whiskers
    geom_boxplot(aes(fill=Tissue2), #dibuja los boxplots
                 alpha=0.85, outlier.size = 0 ) + 
    # geom_jitter(height = 0, width = 0.1, size=5)
    # geom_jitter(aes(color="black"), height = 0, width = 0.1, size=5, stroke=0.2 )  # geom_violin() + geom_jitter(height = 0, width = 0.1)+
    geom_point( aes(color=Tissue2, fill=Tissue2), #, fill=Nuclei
                pch = 21, 
                stroke = 1, size =6, position=position_jitterdodge()) +
    # "Ctrl_PMBSF"="#e2425aff", #rojo "Ctrl_ALBSF"="#e3ab50ff",#amarillo  "EmbWPC_PMBSF"="#00a29aff", #rojo
    
    scale_fill_manual(values=c("#e3ab50ff","#00a29aff","#e2425aff"))+
    scale_color_manual(values=c("#59370eff","#238786ff","#a28b68ff")) +
    # stat_compare_means(comparisons = my_comparisons)+
    labs(title = tpmscounts2$GeneSymbol[tpmscounts2$geneEnsembl == i],
         subtitle = i#,
         # y="log2 (TPM Expression + 1)"
    ) +
    # coord_cartesian(ylim=c(0, max(tpmscounts2$expression[tpmscounts2$geneEnsembl==i]+ 1))) +
    theme_classic(base_size = 25)+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          strip.background = element_blank(), strip.placement = "outside",
          axis.text.x = element_text(angle = 45, hjust=1, colour="black"),
          axis.text.y = element_text(angle = 0, hjust=1, colour="black"),
          axis.title.y =  element_blank(),
          axis.title.x =  element_blank(),
          plot.title=element_text(size=25, hjust=0.5, face="bold.italic", colour="black", vjust=-1),
          plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
          legend.position = "none")
}
for (i in goi) {
  file_name = paste("./ExpressionFigures/TPM_", 
                    unique(tpmscounts2$GeneSymbol[tpmscounts2$geneEnsembl == i]), ".svg", sep="")
  ggsave(file_name, 
         plot_list[[i]], device="svg", dpi = 600, 
         width = 4, height = 5.5, units = "in")
  
}

#Datos Lorenzo P0#############################################

#Genes específicos de las regiones#
resTissue1 <- lapply(openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx"),
                     openxlsx::read.xlsx,xlsxFile="./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx")
names(resTissue1) <- openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322Tissue1.xlsx")
head(resTissue1$ShApe_padj0.05_LFC0.322)


PMBSF_SpC_P0_lor<- resTissue1$ShApe_Fulltable[which(resTissue1$ShApe_Fulltable$padj <= 0.1 &
                                                      resTissue1$ShApe_Fulltable$log2FoldChange > 0.137 ),]
nrow(PMBSF_SpC_P0_lor)
# 377
ALBSF_SpC_P0_lor <- resTissue1$ShApe_Fulltable[which(resTissue1$ShApe_Fulltable$padj <= 0.1 & 
                                                       resTissue1$ShApe_Fulltable$log2FoldChange < -0.137 ),]
nrow(ALBSF_SpC_P0_lor)

#Genes modificados por la cauterización#
ALBSF_embWPCvsCtrl_P0_lor <- lapply(openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx"),
                                    openxlsx::read.xlsx,xlsxFile="./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx")
names(ALBSF_embWPCvsCtrl_P0_lor) <- openxlsx::getSheetNames("./P0/P0_Lorenzo/DEA_padj0.05_LFC0.322State_embWPCvsCtrl_ALBSF.xlsx")
head(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable)

ALBSF_embWPC_P0_lor <- ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable[which(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable$padj <= 0.1 & 
                                                                         abs(ALBSF_embWPCvsCtrl_P0_lor$ShApe_Fulltable$log2FoldChange) > 0.137 ),]
nrow(ALBSF_embWPC_P0_lor)

##Overlapping genes P0 lor##############################################
venP0_lor <-list(ALBSF=ALBSF_SpC_P0_lor$EnsemblID, 
            PMBSF=PMBSF_SpC_P0_lor$EnsemblID, 
            ALBSF_emb_WPC =ALBSF_embWPC_P0_lor$EnsemblID)

diag_vennP0_lor=ggVennDiagram(venP0_lor,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Overlapgenes P0 ALBSF_embWPC") + 
  theme(text = element_text(size = 16)) + scale_fill_gradient(low="blue",high = "red")

png(filename = "./figuras/diagVenn_Overlapping_genes_ALBSF_embWPC_P0_lor.png",
    width=1000, height = 1000)
diag_vennP0_lor
dev.off()
#Lista genes modificados por la cauterizacion en cada region#
ven_intersect_p0_lor <- venn(venP0_lor)
genes_especificos_P0_lor<-attr(ven_intersect_p0_lor,"intersections")
genesPMBSF_ALBSF_emb_WPC_P0_lor<-genes_especificos_P0_lor$`PMBSF:ALBSF_emb_WPC`
genesALBSF_ALBSF_emb_WPC_P0_lor<-genes_especificos_P0_lor$`ALBSF:ALBSF_emb_WPC`

##PMBSF:ALBSF_emb Genes expresados en las dos edades################
genes_overlap_PMBSF_P8<-genesPMBSF_ALBSF_emb_WPC
ven_overlap_PMBSF_P0_P8_lor<-list(P0=genesPMBSF_ALBSF_emb_WPC_P0_lor,
                                  P8=genes_overlap_PMBSF_P8)

diag_ven_overlap_PMBSF_P0_P8_lor<-ggVennDiagram(ven_overlap_PMBSF_P0_P8_lor,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("PMBSFspc:ALBSF_embWPC P0 y P8") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/P0_lorenzo/Lor_diagVenn-overlap_genes_PMBSF_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_PMBSF_P0_P8_lor
dev.off()
pdf("./figuras/P0_lorenzo/Lor_diagVenn-overlap_genes_PMBSF_P0_P8_lfc0137_01padj.pdf",
    width=15, height = 15)
diag_ven_overlap_PMBSF_P0_P8_lor
dev.off()
ven_intersect_PMBSF_P0_P8_lor <- venn(ven_overlap_PMBSF_P0_P8_lor)
genes_comunes_PMBSF_P0_P8_lor <-attr(ven_intersect_PMBSF_P0_P8_lor,"intersections")
genes_comunes_PMBSF_P0_P8_lor<-genes_comunes_PMBSF_P0_P8_lor$`P0:P8`

##ALBSF:ALBSF_emb Genes expresados en las dos edades################
genes_overlap_ALBSF_P8<-genesALBSF_ALBSF_emb_WPC
ven_overlap_ALBSF_P0_P8_lor<-list(P0=genesALBSF_ALBSF_emb_WPC_P0_lor,
                                  P8=genes_overlap_ALBSF_P8)

diag_ven_overlap_ALBSF_P0_P8_lor<-ggVennDiagram(ven_overlap_ALBSF_P0_P8_lor,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("ALBSFspc:ALBSF_embWPC P0 y P8") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_ALBSF_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_ALBSF_P0_P8_lor
dev.off()
pdf("./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_ALBSF_P0_P8_lfc0137_01padj.pdf",
    width=15, height = 15)
diag_ven_overlap_ALBSF_P0_P8_lor
dev.off()
ven_intersect_ALBSF_P0_P8_lor <- venn(ven_overlap_ALBSF_P0_P8_lor)
genes_comunes_ALBSF_P0_P8_lor <-attr(ven_intersect_ALBSF_P0_P8_lor,"intersections")
genes_comunes_ALBSF_P0_P8_lor<-genes_comunes_ALBSF_P0_P8_lor$`P0:P8`

genesSPC_ALBSF_emb_WPC_P0_P8<-c(genes_comunes_PMBSF_P0_P8_lor,genes_comunes_ALBSF_P0_P8_lor)
##Comparacion P0 vs P8 Especificos###################################################
###PMBSFspc Genes expresados en las dos edades################
genes_PMBSFspc_P8<-rownames(PMBSF_SpC)
ven_overlap_PMBSFspc_P0_P8_lor<-list(P0=PMBSF_SpC_P0_lor$EnsemblID,
                                 P8=genes_PMBSFspc_P8)

diag_ven_overlap_PMBSFspc_P0_P8_lor<-ggVennDiagram(ven_overlap_PMBSFspc_P0_P8_lor,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes específicos de PMBSF(PMBSFspc) conservados en P0 y P8") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_PMBSFspc_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_PMBSFspc_P0_P8_lor
dev.off()
pdf( "./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_PMBSFspc_P0_P8_lfc0137_01padj.pdf",
    width=15, height = 15)
diag_ven_overlap_PMBSFspc_P0_P8_lor
dev.off()
ven_intersect_PMBSFspc_P0_P8_lor <- venn(ven_overlap_PMBSFspc_P0_P8_lor)
genes_comunes_PMBSFspc_P0_P8_lor <-attr(ven_intersect_PMBSFspc_P0_P8_lor,"intersections")
genes_comunes_PMBSFspc_P0_P8_lor<-genes_comunes_PMBSFspc_P0_P8_lor$`P0:P8`
###Anotacion
anotacion_genes_comunes_PMBSFspc_P0_P8_lor<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                                              filters = "ensembl_gene_id",
                                              values = genes_comunes_PMBSFspc_P0_P8_lor,
                                              mart = mouse)
colnames(anotacion_genes_comunes_PMBSFspc_P0_P8_lor)[1] <- "EnsemblID"
anotacion_genes_comunes_PMBSFspc_P0_P8_lor<-inner_join(anotacion_genes_comunes_PMBSFspc_P0_P8_lor,spc_P8_vs_spcP0_completa, "EnsemblID")
anotacion_genes_comunes_PMBSFspc_P0_P8_lor<-arrange(anotacion_genes_comunes_PMBSFspc_P0_P8_lor, desc(log2FoldChange))
anotacion_genes_comunes_PMBSFspc_P0_P8_lor_2<-anotacion_genes_comunes_PMBSFspc_P0_P8_lor[,-c(5,6,7,8,9,10,11)]
write.xlsx(anotacion_genes_comunes_PMBSFspc_P0_P8_lor_2, file = "PMBSFspc_P8_P0.xlsx", rowNames = FALSE)
###ALBSFspc Genes expresados en las dos edades################
genes_ALBSFspc_P8<-rownames(ALBSF_SpC)
ven_overlap_ALBSFspc_P0_P8_lor<-list(P0=ALBSF_SpC_P0_lor$EnsemblID,
                                 P8=genes_ALBSFspc_P8)

diag_ven_overlap_ALBSFspc_P0_P8_lor<-ggVennDiagram(ven_overlap_ALBSFspc_P0_P8_lor,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes específicos de ALBSF(ALBSFspc) conservados en P0 y P8") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_ALBSFspc_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_ALBSFspc_P0_P8_lor
dev.off()
pdf( "./figuras/P0_lorenzo/lor_diagVenn-overlap_genes_ALBSFspc_P0_P8_lfc0137_01padj.pdf",
     width=15, height = 15)
diag_ven_overlap_ALBSFspc_P0_P8_lor
dev.off()
ven_intersect_ALBSFspc_P0_P8_lor <- venn(ven_overlap_ALBSFspc_P0_P8_lor)
genes_comunes_ALBSFspc_P0_P8_lor <-attr(ven_intersect_ALBSFspc_P0_P8_lor,"intersections")
genes_comunes_ALBSFspc_P0_P8_lor<-genes_comunes_ALBSFspc_P0_P8_lor$`P0:P8`

#Anotacion
anotacion_genes_comunes_ALBSFspc_P0_P8_lor<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                                                  filters = "ensembl_gene_id",
                                                  values = genes_comunes_ALBSFspc_P0_P8_lor,
                                                  mart = mouse)
colnames(anotacion_genes_comunes_ALBSFspc_P0_P8_lor)[1] <- "EnsemblID"
anotacion_genes_comunes_ALBSFspc_P0_P8_lor<-inner_join(anotacion_genes_comunes_ALBSFspc_P0_P8_lor,spc_P8_vs_spcP0_completa, by="EnsemblID",)
anotacion_genes_comunes_ALBSFspc_P0_P8_lor<-arrange(anotacion_genes_comunes_ALBSFspc_P0_P8_lor, desc(log2FoldChange))
anotacion_genes_comunes_ALBSFspc_P0_P8_lor_2<-anotacion_genes_comunes_ALBSFspc_P0_P8_lor[,-c(5,6,7,8,9,10,11)]
write.xlsx(anotacion_genes_comunes_ALBSFspc_P0_P8_lor_2, file = "ALBSFspc_P8_P0.xlsx", rowNames = FALSE)


##ALBSF_embWPC genes P0 yP8#################################################
ALBSFemb_WPC_P8_genes<-rownames(ALBSF_mbWPC)
ven_overlap_P8yP0_ALBSFembWPC_LOR<-list(ALBSFemb_WPC_P0=ALBSF_embWPC_P0_lor$EnsemblID,
                                    ALBSFemb_WPC_P8=ALBSFemb_WPC_P8_genes)

diag_ven_overlap_P8yP0_ALBSFembWPC_LOR<-ggVennDiagram(ven_overlap_P8yP0_ALBSFembWPC_LOR,label_size = 16, label_alpha = 0.5) + 
  scale_fill_gradient(low="blue",high = "red") + 
  ggtitle("Genes de ALBSF_emb_WPC expresados a P0 y P8") + 
  theme(text = element_text(size = 16))
png(filename = "./figuras/P0_lorenzo/lor_diag_ven_ALBSFembWPC_P0_P8_lfc0137_01padj.png",
    width=1000, height = 1000)
diag_ven_overlap_P8yP0_ALBSFembWPC_LOR
dev.off()

spc_P8_vs_spcP0_completa