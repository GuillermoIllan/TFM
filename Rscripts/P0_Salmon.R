library(dplyr) # data wrangling
library(ggplot2) # plotting
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  
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
setwd("~/Escritorio/datosrecuperados/scripts/Rscripts/AED/P0/P0_analisis/salmon")

#https://nbisweden.github.io/workshop-RNAseq/2011/lab_kallisto.html#21_Get_gene_quantifications

coldata <- read.csv("../metaDataMar_P0.csv",sep="," ,row.names=1)
coldata <- coldata[,c("Tratamiento","Region")]
coldata$Tratamiento <- factor(coldata$Tratamiento)
coldata$Region <- factor(coldata$Region)

files <- paste("salmon_analisis/general" ,  
               list.files(path = "salmon_analisis/general",pattern = "quant.sf", recursive = TRUE),
               sep = "/")
names(files) <- rownames(coldata)

tx2gene <- read_csv("./tx2gene.mm.GRCm39.cdna.csv")
txi.salmon.tsv <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

dds <- DESeqDataSetFromTximport(txi.salmon.tsv, coldata, ~Region +Tratamiento)
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
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rld)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Tratamiento),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) 
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



##Anotación con biomaRt#######################################################

ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                         host =  'https://may2021.archive.ensembl.org')

mouse <- useDataset("mmusculus_gene_ensembl", ensembl_genes)
listEnsemblArchives()

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
head(df_final)
dim(df_final)


#Tissue SPC###########################
controles_coldata<-coldata[1:8,]
filesTissue <- paste("salmon_analisis/Tissue" ,  
                     list.files(path = "salmon_analisis/Tissue",pattern = "quant.sf", recursive = TRUE),
                     sep = "/")
names(filesTissue) <- rownames(controles_coldata)

txi.salmon.tissue.tsv <- tximport(filesTissue, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
##Objeto DESeq #########################################
ddsC <- DESeqDataSetFromTximport(txi.salmon.tissue.tsv, controles_coldata, ~Region)
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
rldC<-rlog(ddsC)

##PCA############################
data<-plotPCA(rldC, intgroup = "Region",returnData=TRUE,
              ntop = 500)


percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
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
##eleccion LFC y padj##############
selected<-which(resC$padj <= 0.1)
resShr_sig <- resC[selected,]
df <- as.data.frame(resShr_sig)
dim(df)

#0.05padj
resC2<-results(ddsC, alpha = 0.05, lfcThreshold= 0)
selected2<-which(resC2$padj <= 0.05)
resShr_sig2 <- resC2[selected2,]
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

#0.5padj
summary(resC2)
calculateDiffExpressed(df2, 0.137)
calculateDiffExpressed(df2, 0.332)
calculateDiffExpressed(df2, 0.5)

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

##Heatmap####################
#library("genefilter")
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

##Volcano###################################################################
df_finalC$diffexpressed <- "NO"
df_finalC$diffexpressed[df_finalC$log2FoldChange > 0.137 & df_finalC$padj < 0.01] <- "UP"
df_finalC$diffexpressed[df_finalC$log2FoldChange < -0.137 & df_finalC$padj < 0.01] <- "DOWN"
table(df_finalC$diffexpressed)
#etiquetamos los genes diferenciados
df_finalC$delabel <- NA
df_finalC$delabel[df_finalC$diffexpressed != "NO"] <- df_finalC$mgi_symbol[df_finalC$diffexpressed != "NO"]
#representamos con ggplot
#library(ggrepel)
png(filename = "figuras/volcano_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalC, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_vline(xintercept=c(-0.322, 0.322), col="red") 
#geom_hline(yintercept=-log10(0.05), col="red")
dev.off()

##MA Plot#######################################################################
png(filename = "figuras/MAplot_ALBSF_ctrl_vs_PMBSF_ctrl_L2FC0322_005padj.png",
    width=1000, height = 1000)
ggplot(data=df_finalC, aes(x=log2(baseMean+1), y=log2FoldChange, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#40E0D0", "grey","#77DD77")) +
  geom_hline(yintercept=c(0, 0), col="red") 
#geom_hline(xintercept=-log10(0.1), col="red")
dev.off()

#ALBSF emb_vs_ctrl##############################################################
ALBSF_coldat<-coldata[1:4,]
coldata_cauterizados<-coldata[9:12,]
ALBSF_coldata<-rbind(ALBSF_coldat,coldata_cauterizados)

filesALBSF <- paste("salmon_analisis/ALBSF_emb_VS_ctrl" ,  
                    list.files(path = "salmon_analisis/ALBSF_emb_VS_ctrl",pattern = "quant.sf", recursive = TRUE),
                    sep = "/")
names(filesALBSF) <- rownames(ALBSF_coldata)

txi.salmon.ALBSF.tsv <- tximport(filesALBSF, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
##Objeto DESeq##########################################
ddsA <- DESeqDataSetFromTximport(txi.salmon.ALBSF.tsv, ALBSF_coldata, ~Tratamiento)

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
png(filename = "figuras/PCA_ALBSF_EMB_VS_ALBSF_ctrl_L2FC0322_005padj.png",
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
#plotMA
plotMA(resA_apeglm, ylim=c(-5,5))

#Anotación con biomaRt

tabla_anotacionA<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id", "gene_biotype","description") ,
                        filters = "ensembl_gene_id",
                        values = rownames(dfA),
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

#PMBSF_EMB_VS_PMBSF_ctrl########################################################
PMBSF_coldat<-coldata[5:8,]
PMBSF_ctrl_coldata<-coldata[13:16,]
PMBSF_coldata<-rbind(PMBSF_coldat,PMBSF_ctrl_coldata)

filesPMBSF <- paste("salmon_analisis/PMBSF_emb_vs_ctrl" ,  
                    list.files(path = "salmon_analisis/PMBSF_emb_vs_ctrl",pattern = "quant.sf", recursive = TRUE),
                    sep = "/")
names(filesPMBSF) <- rownames(PMBSF_coldata)

txi.salmon.PMBSF.tsv <- tximport(filesPMBSF, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
##dds#######################################################
ddsP <- DESeqDataSetFromTximport(txi.salmon.PMBSF.tsv, PMBSF_coldata, ~Tratamiento)
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

ALBSF_SpC <- resC_apglm[which(resC_apglm$padj <= 0.1 & 
                                resC_apglm$log2FoldChange < -0.137 ),]
nrow(ALBSF_SpC)
ALBSF_mbWPC <- resA_apeglm[which(resA_apeglm$padj <= 0.1 & 
                                   abs(resA_apeglm$log2FoldChange) > 0.137 ),]
nrow(ALBSF_mbWPC)

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
saveRDS(genesPMBSF_ALBSF_emb_WPC, file = "../genes_overlapPMBSF-ALBSF_mbWPC_Salmon.rds")
genesALBSF_ALBSF_emb_WPC<-genes_especificos$`ALBSF:ALBSF_mbWPC`
saveRDS(genesALBSF_ALBSF_emb_WPC, file = "../genes_overlapALBSF-ALBSF_mbWPC_Salmon.rds")

##Heatmap representando el overlap############################
#counts
cts_PMBSF_ALBSF_emb_WPC <- tx2gene %>% filter(GENEID %in% genesPMBSF_ALBSF_emb_WPC)

cts_ALBSF_ALBSF_emb_WPC <- tx2gene %>% filter(GENEID %in% genesALBSF_ALBSF_emb_WPC)
cts_genes_overlap<-rbind(cts_PMBSF_ALBSF_emb_WPC,cts_ALBSF_ALBSF_emb_WPC)
dim(cts_genes_overlap)


#coldata
controles_coldata_vs_ALBSF_EMB_coldata<-coldata[1:12,]

filesoverlap <- paste("salmon_analisis/overlap" ,  
                      list.files(path = "salmon_analisis/overlap",pattern = "quant.sf", recursive = TRUE),
                      sep = "/")
names(filesoverlap) <- rownames(controles_coldata_vs_ALBSF_EMB_coldata)

txi.salmon.overlap.tsv <- tximport(filesoverlap, type = "salmon", tx2gene = cts_genes_overlap, ignoreAfterBar = TRUE)
#Creamos el objeto
dds_Spc_EmbWPC <- DESeqDataSetFromTximport(txi.salmon.overlap.tsv, controles_coldata_vs_ALBSF_EMB_coldata, ~ Region + Tratamiento)
###prefiltering############################################################
keep_Spc_EmbWPC <- rowSums(counts(dds_Spc_EmbWPC)) >= 10
dds_Spc_EmbWPC <- dds_Spc_EmbWPC[keep_Spc_EmbWPC,]

genes_especificos_data <- unique(cts_genes_overlap$GENEID)
##Heatmap##########################################################
rld_Spc_EmbWPC<- rlog(dds_Spc_EmbWPC, blind=FALSE)
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[genes_especificos_data, rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
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
##Anotacion Overlap#############
valuescomplex <- genes_especificos_data
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
mat_Spc_EmbWPC<-assay(rld_Spc_EmbWPC)[genes_especificos_data, rownames(controles_coldata_vs_ALBSF_EMB_coldata)]
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
