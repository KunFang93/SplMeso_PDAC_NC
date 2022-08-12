library(Seurat)
library(sctransform)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape)
library(writexl)
library(readxl)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)


# loading the data
load('seurat_anaylsis_combined_clusters_210827.RData')
efg.seu <- readRDS('TotalCells_E9.5_Updated.rds')
efg.seu <- subset(efg.seu, subset = LineageAnnotations %in% c("Spl_Meso"))
nvsp.seu <- combined_cluster
rm(combined_cluster)

# how many common gene
table(rownames(efg.seu) %in% rownames(nvsp.seu))
common_genes <- rownames(efg.seu)[rownames(efg.seu) %in% rownames(nvsp.seu)]

############# 6 vs 0, 1, 13 ##########
nvsp.seu <- subset(nvsp.seu, subset = cluster_name %in% c("inflammatory fibroblasts","myofibroblasts"))
cluster6.spe.markers <- FindMarkers(nvsp.seu, ident.1 = "6", verbose = FALSE)
write.csv(cluster6.spe.markers,'degs_6vs0.1.13_071522.csv')

# data wrangle
nvsp.seu@meta.data$stim <- colsplit(nvsp.seu@meta.data$stim, split = "_", names = c('cell', 'trt','time'))$trt
nvsp.seu@meta.data[,c("integrated_snn_res.0.5","seurat_clusters")] <- list(NULL)

efg.seu@meta.data[,c("nGene","nUMI","Stages","PreviousClusters","PreviousLineages","AaronLineageAnnotations","cluster_name",
                     "S.Score","G2M.Score","Phase","old.ident","CC.Difference","res.0.4","AaronLineageAnnotations1","CellType"
)] <- list(NULL)
colnames(efg.seu@meta.data) = c('orig.ident', 'percent.mt', 'nCount_RNA', 'nFeature_RNA','cluster_name')
efg.seu@meta.data$stim <- 'Spl_Meso'

############# deg in cardiac, othermeso, and splmeso ###########
cluster.genes <- read_excel('markers_all_clusters_072222.xlsx',sheet='Sheet1')
cardiac.genes <- cluster.genes[cluster.genes$cluster=='Cardiac',]$genes
othermeso.genes <- cluster.genes[cluster.genes$cluster=='Other_Meso',]$genes
splmeso.genes <- cluster.genes[cluster.genes$cluster=='Spl_Meso',]$genes

Idents(nvsp.seu) <- as.factor(nvsp.seu@meta.data$stim)
pvsn.markers <- FindMarkers(nvsp.seu,assay = "RNA", ident.1 = "PDAC", ident.2 = "normal",verbose = FALSE)

pvsn.cardiac.genes <- pvsn.markers[cardiac.genes %in% rownames(pvsn.markers),]
pvsn.othermeso.genes <- pvsn.markers[othermeso.genes %in% rownames(pvsn.markers),]
pvsn.splmeso.genes <- pvsn.markers[splmeso.genes %in% rownames(pvsn.markers),]

pvsn.cardiac.genes <- pvsn.cardiac.genes[order(pvsn.cardiac.genes$avg_log2FC),]
pvsn.othermeso.genes <- pvsn.othermeso.genes[order(pvsn.othermeso.genes$avg_log2FC),]
pvsn.splmeso.genes <- pvsn.splmeso.genes[order(pvsn.splmeso.genes$avg_log2FC),]

ave.exp.pvsn <- AverageExpression(nvsp.seu)
ave.exp.pvsn.rna <- ave.exp.pvsn$RNA[,c('normal','PDAC')]
my.breaks <- c(seq(-2, 2, by=0.01))
my.colors <- c(colorRampPalette(colors = c("skyblue","black", "coral1"))(length(my.breaks)))

Heatmap(ave.exp.pvsn.rna.z[rownames(pvsn.splmeso.genes),], name = "z(mat)", 
        col = my.colors, cluster_columns =FALSE,column_split = c('normal','PDAC'),
        show_row_names = TRUE,cluster_rows = FALSE,row_names_gp = gpar(fontsize = 6),
        show_parent_dend_line = FALSE,column_title=NULL)
Heatmap(ave.exp.pvsn.rna.z[rownames(pvsn.cardiac.genes),], name = "z(mat)", 
        col = my.colors, cluster_columns =FALSE,column_split = c('normal','PDAC'),
        show_row_names = TRUE,cluster_rows = FALSE,row_names_gp = gpar(fontsize = 4),
        show_parent_dend_line = FALSE,column_title=NULL)
Heatmap(ave.exp.pvsn.rna.z[rownames(pvsn.othermeso.genes),], name = "z(mat)", 
        col = my.colors, cluster_columns =FALSE,column_split = c('normal','PDAC'),
        show_row_names = TRUE,cluster_rows = FALSE,row_names_gp = gpar(fontsize = 5),
        show_parent_dend_line = FALSE,column_title=NULL)

scatterdiffplot <- function(genes.list,ave.exp.rna){
  library(ggplot2)
  library(ggrepel)
  exp.rna.log <- as.data.frame(log10(ave.exp.rna[rownames(ave.exp.rna) %in% genes.list,]+0.01))
  exp.rna.log$ave_SC <- exp.rna.log$PDAC-exp.rna.log$normal
  exp.rna.log$Gene <- rownames(exp.rna.log)
  sigval <- (abs(quantile(exp.rna.log$ave_SC,0.05))+abs(quantile(exp.rna.log$ave_SC,0.95)))/2
  exp.rna.log$Significant <- ifelse(abs(exp.rna.log$ave_SC) > sigval, "Sig", "Not Sig")
  print(paste('Sig Val',sigval))
  ggplot(exp.rna.log, aes(x = normal, y = PDAC)) + 
    geom_point(aes(color = Significant)) + scale_color_manual(values = c("grey", "red")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = subset(exp.rna.log, abs(ave_SC) > sigval),
      aes(label = Gene),
      size = 2,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps = 20
    ) + geom_abline(slope=1, intercept = 0,linetype=2)
}

scatterdiffplot(splmeso.genes,ave.exp.pvsn.rna)
scatterdiffplot(cardiac.genes,ave.exp.pvsn.rna)
scatterdiffplot(othermeso.genes,ave.exp.pvsn.rna)

######## Eng expression #########
Idents(nvsp.seu) <- nvsp.seu@meta.data$cluster_name
DimPlot(nvsp.seu,split.by = 'stim',label.box=FALSE)

######## slide14 of bioinfo-to-do #########
FeaturePlot(nvsp.seu,features="Eng")


# only keep common genes
nvsp.seu <- nvsp.seu[rownames(nvsp.seu) %in% common_genes,]
efg.seu <- efg.seu[rownames(efg.seu) %in% common_genes,]

nvsp.seu <- SCTransform(nvsp.seu, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

efg.seu <- SCTransform(efg.seu, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

ienp.list <- list(nvsp = nvsp.seu, efg = efg.seu)
features <- SelectIntegrationFeatures(object.list = ienp.list, nfeatures = 3000)
ienp.list <- PrepSCTIntegration(object.list = ienp.list, anchor.features = features)
# require large amount of RAM
ienp.anchors <- FindIntegrationAnchors(object.list = ienp.list, normalization.method = "SCT",
                                       anchor.features = features)

saveRDS(ienp.anchors,'ienp.anchors.rds')
saveRDS(common_genes,'common_genes.rds')

# only can run in server
ienp.combined.sct <- IntegrateData(anchorset = ienp.anchors, 
                                   normalization.method = "SCT",
                                   features.to.integrate = common_genes)
Idents(ienp.combined.sct)<- ienp.combined.sct@meta.data$stim

ienp.combined.sct <- PrepSCTFindMarkers(ienp.combined.sct)
ienp.combined.sct <- readRDS('ienp.combined.sct.rds')

evsn.markers <- FindMarkers(ienp.combined.sct,assay = "SCT", ident.1 = "Spl_Meso", ident.2 = "normal",verbose = FALSE)
evsp.markers <- FindMarkers(ienp.combined.sct,assay = "SCT", ident.1 = "Spl_Meso", ident.2 = "PDAC",verbose = FALSE)
pvsn.markers <- FindMarkers(nvsp.seu,assay = "RNA", ident.1 = "PDAC", ident.2 = "normal",verbose = FALSE)

# read data from server
ave.exp.sct <- AverageExpression(ienp.combined.sct,assays = "SCT",slot='data')$SCT

Eup.Ndw.Pup <- Reduce(intersect, list(rownames(evsn.markers[evsn.markers$avg_log2FC>0,]),
                                      rownames(pvsn.markers[pvsn.markers$avg_log2FC>0,])))
Eup.Nup.Pdw <- Reduce(intersect, list(rownames(evsp.markers[evsp.markers$avg_log2FC>0,]),
                                      rownames(pvsn.markers[pvsn.markers$avg_log2FC<0,])))

plotAveExpHeatmap <- function(ave_exp, features, rowsplit,clusterrows){
  library(ComplexHeatmap)
  ave_exp_sub <- as.data.frame(subset(ave_exp[features,]))
  ave_exp_sub <- as.matrix(ave_exp_sub[c('Spl_Meso','normal','PDAC')])
  my.breaks <- c(seq(-2, 2, by=0.01))
  my.colors <- c(colorRampPalette(colors = c("skyblue","black", "coral1"))(length(my.breaks)))
  ave_exp_sub_z <- t(scale(t(ave_exp_sub)))
  Heatmap(ave_exp_sub_z, name = "mat", col = my.colors, cluster_columns =FALSE,
          show_row_names = TRUE,cluster_rows = clusterrows,
          row_split = rowsplit,
          show_parent_dend_line = FALSE,row_title=NULL,
          row_names_gp = gpar(fontsize = 6),cluster_row_slices = FALSE
  )
}

features.selected <- c(Eup.Ndw.Pup,Eup.Nup.Pdw)
rowsplit.selected <- c(rep("normal",length(Eup.Ndw.Pup)),
                       rep("PDAC",length(Eup.Nup.Pdw)))

Eup.Ndw.Pup.exp <- as.data.frame(ave.exp.sct$RNA[Eup.Ndw.Pup,])
Eup.Nup.Pdw.exp <- as.data.frame(ave.exp.sct$RNA[Eup.Nup.Pdw,])


# heatmap ordered by FC spl_meso.vs.normal (top) and spl_meso.vs.pdac (bottom) 
Eup.Ndw.Pup.exp <- cbind(Eup.Ndw.Pup.exp,evsn.markers[rownames(Eup.Ndw.Pup.exp),])
Eup.Nup.Pdw.exp <- cbind(Eup.Nup.Pdw.exp,evsn.markers[rownames(Eup.Nup.Pdw.exp),])

Eup.Ndw.Pup.FC.srt <- rownames(Eup.Ndw.Pup.exp[order(Eup.Ndw.Pup.exp$avg_log2FC,decreasing = TRUE),])
Eup.Nup.Pdw.FC.srt <- rownames(Eup.Nup.Pdw.exp[order(Eup.Nup.Pdw.exp$avg_log2FC,decreasing = TRUE),])
plotAveExpHeatmap(ave.exp.sct,c(Eup.Ndw.Pup.FC.srt,Eup.Nup.Pdw.FC.srt),rowsplit.selected,FALSE)


