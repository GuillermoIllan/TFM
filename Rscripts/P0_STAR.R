#Análisis a PO#################################################################

library("AnnotationDbi")
library("org.Mm.eg.db")
library("clusterProfiler")
library("DESeq2")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library(biomaRt)
library("genefilter")
library("ggrepel")
library(pheatmap)
library(ashr)
library(ggVennDiagram)
library("dplyr")
library("ComplexHeatmap")
library("readr")
library("tibble")
library("apeglm")
library("openxlsx")
library(plotly)
library(GeneOverlap)
setwd("~/Escritorio/datosrecuperados/scripts/Rscripts/AED/P0/P0_analisis")

##Importar datos################################################################
countsP0 <- as.matrix(read.csv("STAR/FeatureCounts/featurecountstotal_STAR.tsv",sep="",row.names="Geneid"))

coldataP0<-readr::read_csv("./metaDataMar_P0.csv")
coldataP0<- coldataP0 %>% column_to_rownames(var="Muestra")
coldataP0 <- coldataP0[,c("Tratamiento","Region")]
coldataP0$Tratamiento <- factor(coldataP0$Tratamiento)
coldataP0$Region <- factor(coldataP0$Region)

#library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countsP0,
                              colData = coldataP0,
                              design = ~ Region + Tratamiento)
##PRE-FILTERING################################################################
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
colSums(counts(dds))
##PCA##########################################################################
rld <- rlog(dds)
data<-plotPCA(rld, intgroup = c("Tratamiento","Region"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, shape=Region, data=data, label=colnames(rld)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "./figuras/PCA_CAUTvsCTRL_general_L2FC0322_005padj.png",
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

##Normalización################################################################
dds$Tratamiento <- relevel(dds$Tratamiento, "Control")
dds <- DESeq(dds,test = "Wald")

res<- results(dds, alpha = 0.1, lfcThreshold= 0)
resShr<- lfcShrink(dds,res = res, coef="Tratamiento_Cauterizado_vs_Control")
#Tabla resultado
dim(resShr)
head(resShr)
summary(resShr)
sigs<-na.omit(resShr)
sigs<- sigs[sigs$padj < 0.1,]
selected<-which(resShr$padj <= 0.1)
resShr_sig <- resShr[selected,]
df <- as.data.frame(resShr_sig)
dim(df)

#0.05padj
selected2<-which(resShr$padj <= 0.05)
resShr_sig2 <- resShr[selected2,]
df2 <- as.data.frame(resShr_sig2)
dim(df2)
##Anotación con biomaRt#######################################################
#library(biomaRt)
# Create mouse mart object
# listDatasets(ensembl_genes)
ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                         host =  'https://may2021.archive.ensembl.org')

mouse <- useDataset("mmusculus_gene_ensembl", ensembl_genes)
listEnsemblArchives()
# ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "104")
# lista_atributos<-listAttributes(ensembl)
tabla_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(df),
                       mart = mouse)

dim(tabla_anotacion)
dim(df)
# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combined <- merge(df, tabla_anotacion %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#m_cx <- match(rownames(df), tabla_anotacion$ensembl_gene_id)
#out_res_anotobj <- cbind(df[,1:ncol(df)], 
#                        EnsemblID=rownames(df),
#                       GeneSymbol=tabla_anotacion$mgi_symbol[m_cx],
#                       EntrezID=tabla_anotacion$entrezgene_id[m_cx], 
#                       Description=tabla_anotacion$description[m_cx],
#                       Gene_biotype=tabla_anotacion$gene_biotype[m_cx])
#dim(out_res_anotobj)

#eliminamos duplicados
df_final<- unique(df_combined)
rownames(df_final)<-df_final$Row.names
df_final$Row.names<-NULL
head(df_final)
dim(df_final)

#Tissue SPC###########################
#Obtención de los DEGs de los territorios específicos
controles<-countsP0[,1:8]
controles_coldata<-coldataP0[1:8,]
##Objeto DESeq #########################################
ddsC <- DESeqDataSetFromMatrix(countData = controles,
                               colData = controles_coldata,
                               design = ~ Region)

ddsC
##Pre-filtering#############################################
keepC <- rowSums(counts(ddsC)) >= 10
ddsC <- ddsC[keepC,]
##Normalización#########################################################################
ddsC <- DESeq(ddsC)
resC<-results(ddsC, alpha = 0.1, lfcThreshold= 0)
summary(resC)
resC_apglm<- lfcShrink(ddsC,res = resC, coef = "Region_PMBSF_vs_ALBSF")#Tabla resultado
summary(resC_apglm)
head(resC_apglm)
saveRDS(resC_apglm, "SPC_resAPEGLM_P0.rds")
rldC<-rlog(ddsC)

##PCA############################
data<-plotPCA(rldC, intgroup = "Region",returnData=TRUE,
              ntop = 500)

percentVar <- round(100 * attr(data, "percentVar"))
p<-qplot(PC1, PC2, color=Region, data=data, label=colnames(rldC)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "figuras/PCA_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0137_0.1padj.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rldC)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Region),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75) 

dev.off()
##Anotación con biomaRt############

rownames_ordenadosC <- sort(rownames(dfC))
datos_ordenados_C <- dfC[rownames_ordenadosC, ]
valuesC <- rownames(datos_ordenados_C)
tabla_anotacionC<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = valuesC,
                        mart = mouse)
#tabla con todos los genes, sin filtrar
resC_apglm_sinNA <- na.omit(resC_apglm)
tabla_anotacionC_completo<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                        filters = "ensembl_gene_id",
                        values = rownames(resC_apglm_sinNA),
                        mart = mouse)

# Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combinedC <- merge(dfC, tabla_anotacionC %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)
resC_apglm_Anotado<- merge(as.data.frame(resC_apglm), tabla_anotacionC_completo %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_finalC<- unique(df_combinedC)
rownames(df_finalC)<-df_finalC$Row.names
df_finalC$Row.names<-NULL

##eleccion LFC y padj##############
##0.1padj dataframe
selected<-which(resC_apglm$padj <= 0.1)
resShr_sig <- resC_apglm[selected,]
df <- as.data.frame(resShr_sig)
dim(df)

#0.05padj
selected2<-which(resC_apglm$padj <= 0.05)
resShr_sig2 <- resC_apglm[selected2,]
df2 <- as.data.frame(resShr_sig2)
dim(df2)

resC2<-results(ddsC, alpha = 0.05, lfcThreshold= 0)

###numero de upregulados y downregulados############################
calculateDiffExpressed <- function(df, threshold) {
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > threshold] <- "UP"
  df$diffexpressed[df$log2FoldChange < -threshold] <- "DOWN"
  return(table(df$diffexpressed))
}

# Uso de la función con los tres casos
summary(resC_apglm)
calculateDiffExpressed(dfC, 0.137)
calculateDiffExpressed(dfC, 0.332)
calculateDiffExpressed(dfC, 0.5)


###0.05padj numero de upregulados y downregulados############################
summary(resC2)
calculateDiffExpressed(df2, 0.137)
calculateDiffExpressed(df2, 0.332)
calculateDiffExpressed(df2, 0.5)

##Heatmap simple####################

#filtramos los genes para quedarnos con los significativos
df.topC <- df_finalC[(abs(df_finalC$log2FoldChange) > 0.137),]
df.topC <- df.topC[order(df.topC$log2FoldChange, decreasing = TRUE),]
df.topC$diffexpressed <- "NO"
df.topC$diffexpressed[df.topC$log2FoldChange > 0.137 ] <- "UP"
df.topC$diffexpressed[df.topC$log2FoldChange < -0.137 ] <- "DOWN"
table(df.topC$diffexpressed)

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
png(filename = "figuras/Heatmap_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
heatmap.2(matC.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()
##Complex Heatmap##############################################################################
#anotación en symbol
dfC_anot<-getBM(attributes =c("ensembl_gene_id","mgi_symbol","gene_biotype") ,
                filters = "ensembl_gene_id",
                values = rownames(dfC),
                mart = mouse)
dfc2<-rownames_to_column(dfC,var = "ensembl_gene_id")
dfC_biotype <- inner_join(dfc2,dfC_anot,"ensembl_gene_id")
dfC_filtrado <- dfC_biotype %>% filter(gene_biotype == "protein_coding")

#50 genes más expresados diferencialmente
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

#Posicionar la anotación en la matriz del heatmap
posiciones_matrizC <- which(!is.na(genes_anotacion_C))

#asigna los valores de top_50_genesC$mgi_symbol al vector label_genesC basándose en los elementos no nulos de genes_anotacion_C.
label_genesC <- vector("character", length(genes_anotacion_C))
for (i in 1:length(genes_anotacion_C)) {
  if (!is.na(genes_anotacion_C[i])) {
    label_genesC[i] <- top_50_genesC$mgi_symbol[genes_anotacion_C[i]]
  }
}
label_genesC<-subset(label_genesC, nzchar(label_genesC))
print(label_genesC)
#anotación comples heatmap
anotaciones_c = rowAnnotation(foo = anno_mark(at = posiciones_matrizC,
                                              labels =label_genesC))

col_fun= colorRamp2(seq(-1.5,1.5, length.out=100),
                    (colorRampPalette(c("darkorchid1","black","yellow1"))(100)),
                    space="RGB")
png(filename="figuras/Heatmap_ALBSF_ctrl_VS_PMBSF_ctrl_01padjP0.png",
    width=600, height = 600)
Heatmap(matC.scaled,
        col=col_fun,
        na_col = "black",
        column_title = "PMBSFspc(391)_vs_ALBSFspc(392)",
        column_title_gp = gpar(fontsize = 18),
        right_annotation = anotaciones_c,
        row_names_gp = gpar(fontsize = 14),
        show_row_names = FALSE,
        cluster_columns = TRUE,)
dev.off()
##Volcano###################################################################
resC_apglm_Anotado$diffexpressed <- "NO"
resC_apglm_Anotado$diffexpressed[resC_apglm_Anotado$log2FoldChange > 0.137 & resC_apglm_Anotado$padj < 0.1] <- "UP"
resC_apglm_Anotado$diffexpressed[resC_apglm_Anotado$log2FoldChange < -0.137 & resC_apglm_Anotado$padj < 0.1] <- "DOWN"
table(resC_apglm_Anotado$diffexpressed)

#etiquetamos los genes diferenciados
resC_apglm_Anotado$delabel <- NA
resC_apglm_Anotado$delabel[resC_apglm_Anotado$diffexpressed != "NO"] <- resC_apglm_Anotado$mgi_symbol[resC_apglm_Anotado$diffexpressed != "NO"]

#representamos con ggplot
png(filename = "figuras/volcano_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0137_01padj.png",
    width=900, height = 900)
ggplot(data=resC_apglm_Anotado, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point(size=3.5) + 
  theme_minimal() +
  geom_text_repel(min.segment.length = 0,size=7.5) +
  scale_color_manual(values=c("#00a4a2", "grey","#b05d95")) +
  geom_vline(xintercept=c(-0.137, 0.137), col="black",linetype="dashed") +
  geom_hline(yintercept=1, col="black",linetype="dashed") +
  theme(axis.title.x = element_text(size=12))+
  theme(axis.title.y = element_text(size=12))+
  theme(legend.title = element_text(size=4))+
  theme(legend.text = element_text(size=4))+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

##MA Plot#######################################################################
png(filename = "figuras/MAplot_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0137_01padj.png",
    width=800, height = 800)
ggplot(data=resC_apglm_Anotado, aes(x=log2(baseMean+1), y=log2FoldChange, col=diffexpressed, label=delabel)) +
  geom_point(size=3.5) + 
  theme_minimal() +
  geom_text_repel(min.segment.length = 0,size=7.5) +
  scale_color_manual(values=c("#00a4a2", "grey","#b05d95")) +
  geom_hline(yintercept=c(0.137, -0.137), col="black",linetype="dashed") +
  theme(axis.title.x = element_text(size=12))+
  theme(axis.title.y = element_text(size=12))+
  theme(legend.title = element_text(size=4))+
  theme(legend.text = element_text(size=4))+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

#ALBSF emb_vs_ctrl##############################################################
ALBSF<-countsP0[,1:4]
cauterizados<-countsP0[,9:12]
ALBSF_cts<-cbind(ALBSF,cauterizados)

ALBSF_coldat<-coldataP0[1:4,]
coldata_cauterizados<-coldataP0[9:12,]
ALBSF_coldata<-rbind(ALBSF_coldat,coldata_cauterizados)

##Objeto DESeq##########################################
ddsA <- DESeqDataSetFromMatrix(countData = ALBSF_cts,
                               colData = ALBSF_coldata,
                               design = ~ Tratamiento)

##Pre-filtering###########################################
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
png(filename = "PCA_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rldA)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75
) 
dev.off()

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
rownames_ordenadosA <- sort(rownames(dfA))
datos_ordenados_A <- dfA[rownames_ordenadosA, ]
valuesA <- rownames(dfA)

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
png(filename = "figuras/Heatmap_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
heatmap.2(matA.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

##Volcano plot con ggplot################################################
#seleccionamos los genes expresados diferencialmente
df_finalA$diffexpressed <- "NO"
df_finalA$diffexpressed[df_finalA$log2FoldChange > 0.322 ] <- "UP"
df_finalA$diffexpressed[df_finalA$log2FoldChange < -0.322] <- "DOWN"
table(df_finalA$diffexpressed)

#etiquetamos los genes diferenciados
df_finalA$delabel <- NA
df_finalA$delabel[df_finalA$diffexpressed != "NO"] <- df_finalA$mgi_symbol[df_finalA$diffexpressed != "NO"]

#representamos con ggplot
png(filename = "figuras/volcano_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalA, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.322, 0.322), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

##MAplot#########################################################
png(filename = "figuras/MAplot_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalA, aes(x=log2(baseMean+1), y=log2FoldChange, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_hline(yintercept=c(0, 0), col="red") 
  #geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

##PMBSF emb_vs_ctrl#############################################################
#PMBSF_EMB_VS_PMBSF_ctrl########################################################
PMBSF_ctrl<-countsP0[,5:8]
PMBSF_emb_WPC<-countsP0[,13:16]
PMBSF_cts<-cbind(PMBSF_ctrl,PMBSF_emb_WPC)

PMBSF_coldat<-coldataP0[5:8,]
PMBSF_ctrl_coldata<-coldataP0[13:16,]
PMBSF_coldata<-rbind(PMBSF_coldat,PMBSF_ctrl_coldata)
##dds#######################################################
ddsP <- DESeqDataSetFromMatrix(countData = PMBSF_cts,
                               colData = PMBSF_coldata,
                               design = ~ Tratamiento)

ddsP

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
png(filename = "figuras/PCA-PMBSF_EMB_VS_PMBSF_ctrlL2FC0322_005padj.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rldP)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75
) 
dev.off()

##Expresión diferencial#########################################################
ddsP$Tratamiento <- relevel(ddsP$Tratamiento, "Control")
ddsP <- DESeq(ddsP)
resP<-results(ddsP,alpha=0.05, lfcThreshold= 0)
resP_apeglm<- lfcShrink(ddsP,res=resP,coef = "Tratamiento_Cauterizado_vs_Control")#Tabla resultado
summary(resP_apeglm)

#Filtrar los datos
sigsP <- na.omit(resP_apeglm)
sigsP<- sigsP[sigsP$padj < 0.05,]
dfP <- as.data.frame(sigsP)
dim(dfP)
summary(sigsP)
#plotMA
plotMA(resP_apeglm, ylim=c(-5,5))

##Anotación con biomaRt########################################################
#library(biomaRt)
#lista_atributos<-listAttributes(ensembl)
tabla_anotacionP<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = rownames(dfP),
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
df.topP <- df_finalP[(abs(df_finalP$log2FoldChange) > 0.322),]
df.topP <- df.topP[order(df.topP$log2FoldChange, decreasing = TRUE),]
df.topP$diffexpressed <- "NO"
df.topP$diffexpressed[df.topP$log2FoldChange > 0.322 ] <- "UP"
df.topP$diffexpressed[df.topP$log2FoldChange < -0.322] <- "DOWN"
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
png(filename = "figuras/Heatmap-PMBSF_EMB_VS_PMBSF_ctrlL2FC0322_005padj.png",
    width=1000, height = 1000)
heatmap.2(matP.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

##Volcano############################################
df_finalP$diffexpressed <- "NO"
df_finalP$diffexpressed[df_finalP$log2FoldChange > 0.322 ] <- "UP"
df_finalP$diffexpressed[df_finalP$log2FoldChange < -0.322] <- "DOWN"

#etiquetamos los genes diferenciados
df_finalP$delabel <- NA
df_finalP$delabel[df_finalP$diffexpressed != "NO"] <- df_finalP$mgi_symbol[df_finalP$diffexpressed != "NO"]
#representamos con ggplot
#library(ggrepel)
png(filename = "figuras/VOLCANO-PMBSF_EMB_VS_PMBSF_ctrlL2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalP, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.322, 0.322), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

png(filename = "figuras/MAplot-PMBSF_EMB_VS_PMBSF_ctrlL2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalP, aes(x=log2(baseMean+1), y=log2FoldChange, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_hline(yintercept=c(0, 0), col="red") 
#geom_hline(xintercept=-log10(0.1), col="red")
dev.off()

#Overlapping/Venn Diag####################################
PMBSF_SpC <- resC_apglm[which(resC_apglm$padj <= 0.1 & 
                                resC_apglm$log2FoldChange > 0.137 ),]
nrow(PMBSF_SpC)
saveRDS(rownames(PMBSF_SpC),"PMBSF_SpC_P0.rds")
ALBSF_SpC <- resC_apglm[which(resC_apglm$padj <= 0.1 & 
                                resC_apglm$log2FoldChange < -0.137 ),]
nrow(ALBSF_SpC)
saveRDS(rownames(ALBSF_SpC),"ALBSF_SpC_P0.rds")
ALBSF_mbWPC <- resA_apeglm[which(resA_apeglm$padj <= 0.1 & 
                                   abs(resA_apeglm$log2FoldChange) > 0.137 ),]
nrow(ALBSF_mbWPC)
saveRDS(rownames(ALBSF_SpC),"ALBSF_embWPC_P0.rds")
venALBSF_emb_WPC <-list(ALBSF=rownames(ALBSF_SpC), 
                        PMBSF=rownames(PMBSF_SpC), 
                        ALBSF_mbWPC = rownames(ALBSF_mbWPC))

diag_vennALBSF_emb_WPC=ggVennDiagram(venALBSF_emb_WPC) + scale_fill_gradient(low="blue",high = "red")
png(filename = "figuras/DiagVenn_ALBSF_EMB_VS_ALBSF_ctrl_vs_PMBSFctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
diag_vennALBSF_emb_WPC
dev.off()

##Fisher's Test##################################################
go.objALBSF_SpC <- newGeneOverlap(rownames(ALBSF_SpC),
                                  rownames(ALBSF_mbWPC),
                                  nrow(dds) )
go.objALBSF_SpC <- testGeneOverlap(go.objALBSF_SpC)
print(go.objALBSF_SpC)

go.objPMBSF_SpC<- newGeneOverlap(rownames(PMBSF_SpC),
                                 rownames(ALBSF_mbWPC),
                                 nrow(dds) )
go.objPMBSF_SpC <- testGeneOverlap(go.objPMBSF_SpC)
print(go.objPMBSF_SpC)
##seleccionamos los genes específicos##########################
ven_intersect <- venn(venALBSF_emb_WPC)
genes_especificos<-attr(ven_intersect,"intersections")
genesPMBSF_ALBSF_emb_WPC<-genes_especificos$`PMBSF:ALBSF_mbWPC`
saveRDS(genesPMBSF_ALBSF_emb_WPC, file = "overlap_genes_PMBSF:ALBSF_mbWPC_STAR.rds")
genesALBSF_ALBSF_emb_WPC<-genes_especificos$`ALBSF:ALBSF_mbWPC`
saveRDS(genesALBSF_ALBSF_emb_WPC, file = "overlap_genes_ALBSF:ALBSF_mbWPC_STAR.rds")
#SPCvsALBSF_EMB_WPC#########################################
#counts
cts_PMBSF_ALBSF_emb_WPC<-countsP0[genesPMBSF_ALBSF_emb_WPC,]
cts_ALBSF_ALBSF_emb_WPC<-countsP0[genesALBSF_ALBSF_emb_WPC,]
dim(cts_PMBSF_ALBSF_emb_WPC)
cts_genes_overlap<-rbind(cts_PMBSF_ALBSF_emb_WPC,cts_ALBSF_ALBSF_emb_WPC)
saveRDS(cts_genes_overlap, file = "overlap_genes_STAR.rds")
controles1<-cts_genes_overlap[,1:8]
ALBSF_EMB1<-cts_genes_overlap[,9:12]
genes_especificos_data <- cbind(controles1,ALBSF_EMB1)
#coldata
controles_coldata_vs_ALBSF_EMB_coldata<-coldataP0[1:12,]
#Creamos el objeto
dds_Spc_EmbWPC <- DESeqDataSetFromMatrix(countData = genes_especificos_data,
                                         colData = controles_coldata_vs_ALBSF_EMB_coldata,
                                         design = ~ Tratamiento + Region)

###prefiltering############################################################
keep_Spc_EmbWPC <- rowSums(counts(dds_Spc_EmbWPC)) >= 10
dds_Spc_EmbWPC <- dds_Spc_EmbWPC[keep_Spc_EmbWPC,]

##PCA########
rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC, blind=FALSE)
data<-plotPCA(rld_Spc_EmbWPC, intgroup = c("Tratamiento","Region"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Tratamiento, shape=Region, data=data, label=colnames(rld_Spc_EmbWPC)) +
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
##3d PCAplot#####################################################
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
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[rownames(genes_especificos_data), rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
mat_Spc_EmbWPC.scaled <- t(apply(mat_Spc_EmbWPC, 1, scale)) #center and scale each column (Z-score) then transpose
mat_Spc_EmbWPC.scaled<-(mat_Spc_EmbWPC - rowMeans(mat_Spc_EmbWPC))/rowSds(mat_Spc_EmbWPC)
colnames(mat_Spc_EmbWPC.scaled) <- paste0(rld_Spc_EmbWPC$Tratamiento,"-",rld_Spc_EmbWPC$Region)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld_Spc_EmbWPC$Tratamiento ]
png(filename = "figuras/Heatmap_ALBSF_EMB_VS_ALBSF_ctrl_vs_PMBSFctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
heatmap.2(mat_Spc_EmbWPC.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="none")
dev.off()
##Complex Heatmap##############################################################################
ALBSFembWPC_SPC<-rownames_to_column(as.data.frame(cts_genes_overlap),"ensembl_gene_id")
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
        column_title = "Genes específicos modificados por la cauterización 66_PMBSF 86_ALBSF ",
        column_title_gp = gpar(fontsize = 18),
        right_annotation = anotaciones_h,
        row_names_gp = gpar(fontsize = 12),
        show_row_names = FALSE,
        cluster_columns = FALSE)
dev.off()


##Anotacion Overlap#############
valuescomplex <- rownames(genes_especificos_data)
#lista_atributos<-listAttributes(ensembl)
tabla_anotacionComplex<-getBM(attributes =c("ensembl_gene_id","mgi_symbol") ,
                              filters = "ensembl_gene_id",
                              values = valuescomplex,
                              mart = mouse)
data_anotadaComplex <- merge(genes_especificos_data, tabla_anotacionComplex %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE, sort=FALSE)
data_anotadaComplex

###genes mas importantes################

genes_mas_expresadosPMBSF <- resP_apeglm[rownames(resP_apeglm) %in% genesPMBSF_ALBSF_emb_WPC, ]
genes_mas_expresadosPMBSF<-genes_mas_expresadosPMBSF[order(-genes_mas_expresadosPMBSF$log2FoldChange),]
genes_mas_expresadosPMBSF<- head(genes_mas_expresadosPMBSF, n = 15)

genes_mas_expresadosALBSF <- resP_apeglm[rownames(resP_apeglm) %in% genesALBSF_ALBSF_emb_WPC, ]
genes_mas_expresadosALBSF<-genes_mas_expresadosALBSF[order(-genes_mas_expresadosALBSF$log2FoldChange),]
genes_mas_expresadosALBSF<- head(genes_mas_expresadosALBSF, n = 15)

genes_mas_expresados<-rbind(genes_mas_expresadosALBSF,genes_mas_expresadosPMBSF)
##Complex Heatmap####################################
rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC, blind=FALSE)
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[rownames(genes_especificos_data), rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
mat_Spc_EmbWPC.scaled <- t(apply(mat_Spc_EmbWPC, 1, scale)) #center and scale each column (Z-score) then transpose
mat_Spc_EmbWPC.scaled<-(mat_Spc_EmbWPC - rowMeans(mat_Spc_EmbWPC))/rowSds(mat_Spc_EmbWPC)
genes_anotacion_m<-match(rownames(mat_Spc_EmbWPC.scaled), rownames(genes_mas_expresados))


rownames(mat_Spc_EmbWPC.scaled)<- data_anotadaComplex$mgi_symbol 
colnames(mat_Spc_EmbWPC.scaled) <- paste0(rld_Spc_EmbWPC$Tratamiento,"-",rld_Spc_EmbWPC$Region)

anotaciones_h<- rowAnnotation(foo=anno_mark(at = genes_anotacion_m, 
                                            labels = rownames(mat_Spc_EmbWPC.scaled)),
                              width = unit(1, "cm") +
                                max_text_width(rownames(mat_Spc_EmbWPC.scaled)))

data_anotadaComplex$mgi_symbol[genes_anotacion_m]
png(filename = "figuras/Heatmap_ALBSF_EMB_VS_ALBSF_ctrl_VS_PMBSF_ctrl_L2FC0322_005padj.png",
    width=600, height = 1000)
Heatmap(mat_Spc_EmbWPC.scaled, 
        column_title = "ALBSF_EMB_VS_ALBSF_ctrl_VS_PMBSF_ctrl_P0",
        column_title_gp = gpar(fill = "white", col = "black", border = "white"),
        right_annotation = anotaciones_h,
        row_names_gp = gpar(fontsize = 0))
dev.off()

#Comparación alineadores################################################
##Salmon#################
overlap_genes_ALBSF_ALBSF_mbWPC_Salmon<-readRDS(file = "genes_overlapALBSF-ALBSF_mbWPC_Salmon.rds")
overlap_genes_ALBSF_ALBSF_mbWPC_Salmon <- gsub("\\..*", "", overlap_genes_ALBSF_ALBSF_mbWPC_Salmon)
overlap_genes_PMBSF_ALBSF_mbWPC_Salmon<-readRDS(file = "genes_overlapPMBSF-ALBSF_mbWPC_Salmon.rds")
overlap_genes_PMBSF_ALBSF_mbWPC_Salmon <- gsub("\\..*", "", overlap_genes_PMBSF_ALBSF_mbWPC_Salmon)
overlap_genesSalmon<-c(overlap_genes_ALBSF_ALBSF_mbWPC_Salmon,overlap_genes_PMBSF_ALBSF_mbWPC_Salmon)
##Kallisto######################################################################
overlap_genes_ALBSF_ALBSF_mbWPC_Kallisto<-readRDS(file = "genes_overlapALBSF-ALBSF_mbWPC_Kallisto.rds")
overlap_genes_ALBSF_ALBSF_mbWPC_Kallisto <- gsub("\\..*", "", overlap_genes_ALBSF_ALBSF_mbWPC_Kallisto)
overlap_genes_PMBSF_ALBSF_mbWPC_Kallisto<-readRDS(file = "genes_overlapPMBSF-ALBSF_mbWPC_Kallisto.rds")
overlap_genes_PMBSF_ALBSF_mbWPC_Kallisto <- gsub("\\..*", "", overlap_genes_PMBSF_ALBSF_mbWPC_Kallisto)
overlap_genesKallisto<-c(overlap_genes_ALBSF_ALBSF_mbWPC_Kallisto,overlap_genes_PMBSF_ALBSF_mbWPC_Kallisto)

##General overlap##################################################
venAlineadores <-list(Kallisto_127=overlap_genesKallisto, 
                        Salmon_137=overlap_genesSalmon, 
                        STAR_152 = rownames(cts_genes_overlap))

diag_venAlineadores=ggVennDiagram(venAlineadores) + scale_fill_gradient(low="blue",high = "red")
png(filename = "figuras/overlap_general_alineadores.png",
    width=1000, height = 1000)
diag_venAlineadores
dev.off()

#Enriquecimiento overlap###############################################
##STAR########################################
emb_WPC_LFC<-subset(resA_apeglm,rownames(resA_apeglm) %in% rownames(cts_genes_overlap) )
emb_WPC_LFC<-rownames_to_column(as.data.frame(emb_WPC_LFC), var="ensembl_gene_id")
emb_WPC_LFC<- inner_join(emb_WPC_LFC,tabla_anotacionComplex,"ensembl_gene_id") 
GO_resultsgenes_overlapP0_2<- enrichGO(gene = emb_WPC_LFC$ensembl_gene_id, 
                                       OrgDb = "org.Mm.eg.db", 
                                       keyType = "ENSEMBL", 
                                       ont = "BP",
                                       readable=TRUE,
                                       pvalueCutoff = 0.2,
                                       pAdjustMethod = "BH")
#dotplot
png(filename = "./figuras/dotplotBP_overlapping_genes_P0_L2FC0137_01padj.png",
    width=600, height = 600)
dotplot(GO_resultsgenes_overlapP0_2,showCategory=14,font.size=15)
dev.off()

#cnetplot
genes_emb_WPC_LFC<-emb_WPC_LFC$log2FoldChange#genes con el LFC anotado
names(genes_emb_WPC_LFC) <- emb_WPC_LFC$mgi_symbol#genes anotados es Symbol
png(filename = "figuras/cnetplotBP_overlapping_L2FC0137_01padj_P0.png",
    width=900, height = 900)
cnetplot(GO_resultsgenes_overlapP0_2,color.params = list(foldChange=genes_emb_WPC_LFC),showCategory = 14)+scale_color_gradient2(name='log2FoldChange', low='darkgreen', high='firebrick')
dev.off()

gç


###Spc P0 genes################################################
PMBSF_SpC_genes<-rownames(PMBSF_SpC)
ALBSF_SpC_genes<-rownames(ALBSF_SpC)

#compareCluster
compare_spc_genesP8<-list(PMBSF_SpC_genes,ALBSF_SpC_genes)
names(compare_spc_genesP8) <- c("PMBSF", "ALBSF")
ck_spc_genesP8<- compareCluster(geneCluster = compare_spc_genesP8,
                                fun = enrichGO,
                                OrgDb = "org.Mm.eg.db", 
                                keyType = "ENSEMBL", 
                                ont = "BP",
                                readable=TRUE,
                                pvalueCutoff = 0.1,
                                pAdjustMethod = "BH")

png(filename = "./figuras/cnetplotBP_SPC_genes_P0_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(ck_spc_genesP8, showCatgory=10,cex.params=list( category_label = 2, gene_label = 2))
dev.off()

png(filename = "./figuras/dotplotBP_SPC_genes_P0_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(ck_spc_genesP8,  showCategory = 10, by="geneRatio")
dev.off()

##Kallisto####################
GO_resultsgenes_overlapKal<- enrichGO(gene = overlap_genesKallisto, 
                                   OrgDb = "org.Mm.eg.db", 
                                   keyType = "ENSEMBL", 
                                   ont = "BP",
                                   readable=TRUE,
                                   pvalueCutoff = 0.1,
                                   pAdjustMethod = "BH")
write.csv(GO_resultsgenes_overlapKal,"GO_resultsgenes_overlap_Kalllisto.csv")
png(filename = "kallisto/figuras/cnetplotBP_overlapping_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(GO_resultsgenes_overlapKal)
dev.off()
png(filename = "kallisto/figuras/dotplotBP_overlapping_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(GO_resultsgenes_overlapKal,  showCategory = 10)
dev.off()

##Salmon####################
GO_resultsgenes_overlapSal<- enrichGO(gene = overlap_genesSalmon, 
                                      OrgDb = "org.Mm.eg.db", 
                                      keyType = "ENSEMBL", 
                                      ont = "BP",
                                      readable=TRUE,
                                      pvalueCutoff = 0.1,
                                      pAdjustMethod = "BH")

write.csv(GO_resultsgenes_overlapSal, "GO_resultsgenes_overlapSalmon.csv")
png(filename = "salmon/figuras/cnetplotBP_overlapping_L2FC0137_01padj.png",
    width=1000, height = 1000)
cnetplot(GO_resultsgenes_overlapSal)
dev.off()

png(filename = "salmon/figuras/dotplotBP_overlapping_L2FC0137_01padj.png",
    width=1000, height = 1000)
dotplot(GO_resultsgenes_overlapSal,  showCategory = 10)
dev.off()

#TPMS######################################
cts_length<-read.table("STAR/FeatureCounts/featurecountstotal_STAR.txt",header = TRUE)
rownames(cts_length) <- cts_length$Geneid
head(cts_length)

cts_length2<-cts_length[,-c(1,2,3,4,5,6)]
head(cts_length2)
# cts_Count_TPM <- merge(cts_Count_TPM, cts_length, by = "row.names", all = FALSE)
# cts_Count_TPM<-column_to_rownames(cts_Count_TPM,var = "Row.names")

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


##Boxplots TPM###################
genes85_65<-read_table("/home/guillermo/Escritorio/datosrecuperados/scripts/Rscripts/AED/85_65genesoverlaping.txt")
coldataP0_tpm<-as.data.frame(coldataP0)
coldataP0_tpm$Names<- rownames(coldataP0_tpm)

dds_tpm <- DESeqDataSetFromMatrix(countData = countsP0,
                              colData = coldataP0_tpm,
                              design = ~ Region + Tratamiento)

deseqdataCPAEA<-dds_tpm[,!grepl(paste(c("EP"), collapse = "|"),
                                dds_tpm$Names)]
as.data.frame(colData(deseqdataCPAEA))

deseqdataTPM <- deseqdataCPAEA
goi <- genes85_65$EnsemblID

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
#Scatter plot###############################################################
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
  scale_color_manual(values=c("red", "black","darkgreen"))+
  geom_vline(xintercept=0, col="black",linetype="dashed") +
  geom_hline(yintercept=0, col="black",linetype="dashed")+
  ylim(-3,3)+
  xlim(-3,3)+
  geom_text_repel(data=ALBSF_embWPC_vs_Tissue_completa[ALBSF_embWPC_vs_Tissue_completa$label!="Neg",]) +
  theme_dark() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "ALBSF_embWPC:SPC P0 ")+
  theme_classic(base_size = 25 ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10))
png(filename="./figuras/scatterplot_ALBSF_embWPC:SPC_P0.png",    
    width=1000, height = 1000)
fourway
dev.off()