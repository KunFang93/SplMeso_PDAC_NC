library(Seurat)
library(SeuratDisk)
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
library(collections)
library(Matrix)

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  output <- dict()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      # in case cannot find in human
      if (length(human_genes)==0){
        output$set(gene, gene)
      }else{
        output$set(gene, human_genes)
        }
      }
      # for(human_gene in human_genes){
      #   output = append(output,human_gene)
      # }
     else {
      output$set(gene, gene)
    }
  }
  return (output)
}

# load dataset
fibroX <- readRDS('Human_PS_Fibro.RDS')
fibroX.pdac <- subset(fibroX, subset = Tissue %in% c('PDAC'))
fibroX.pdac@project.name <- "FibroX.pdac"
rm(fibroX)

# Convert("/Users/kunfang/Documents/lab/Jin_lab/Collaberations/Leone_lab/PDAC/data/TS_Pancreas.h5ad", 
#         dest = "h5seurat", overwrite = TRUE)
TS.seu <- LoadH5Seurat("TS_Pancreas.h5seurat",assays = "RNA")
TS.seu.fibro <- subset(TS.seu, subset = free_annotation %in% c('fibroblast','pancreatic stellate cell'))
TS.seu.fibro@project.name <- "TS.normal"
rm(TS.seu)

efg.seu <- readRDS('TotalCells_E9.5_Updated.rds')
efg.seu <- subset(efg.seu, subset = LineageAnnotations %in% c("Spl_Meso"))

human.list <- convert_mouse_to_human(rownames(efg.seu@assays$RNA))
human.list <- human.list$as_list()

efg.counts <- GetAssayData(efg.seu,assay = "RNA",slot = "counts") 
human.all.list <- c()
all.names <- c()
for (gene in names(human.list)){
  human_genes = human.list[gene][[1]]
  human.all.list = append(human.all.list,human_genes)
  all.names <- c(all.names, rep(gene,length(human_genes)))
}
names(human.all.list) <- all.names


efg.human.count <- efg.counts[names(human.all.list),]
rownames(efg.human.count) <- human.all.list

# get the mean value for rows with same gene name
efg.human.count.dup = efg.human.count[duplicated(rownames(efg.human.count))|duplicated(rownames(efg.human.count), fromLast = TRUE),]
efg.human.count.dup.mean= aggregate(efg.human.count.dup, list(row.names(efg.human.count.dup)), mean)
rownames(efg.human.count.dup.mean) <- efg.human.count.dup.mean$Group.1
efg.human.count.dup.mean$Group.1 <- NULL
efg.human.count.dup.mean.spar <- Matrix(as.matrix(efg.human.count.dup.mean), sparse=T)

efg.human.count.nodup = efg.human.count[!(duplicated(rownames(efg.human.count))|duplicated(rownames(efg.human.count), fromLast = TRUE)),]

merge.sparse = function(M.list) {
  A = M.list[[1]]
  
  for (i in 2:length(M.list)){ #i indexes of matrices
    # finding what's missing
    misA = colnames(M.list[[i]])[!colnames(M.list[[i]]) %in% colnames(A)]
    misB = colnames(A)[!colnames(A) %in% colnames(M.list[[i]])]
    
    misAl = as.vector(numeric(length(misA)), "list")
    names(misAl) = misA
    misBl = as.vector(numeric(length(misB)), "list")
    names(misBl) = misB
    
    ## adding missing columns to initial matrices
    An = Reduce(cbind, c(A, misAl))
    if (length(misA) > 0)
    {
      lenA <- ncol(An)-length(misA)+1
      colnames(An)[lenA:ncol(An)] = names(misAl)
    }
    
    Bn = Reduce(cbind, c(M.list[[i]], misBl))
    if(length(misB) > 0)
    {
      lenB <- ncol(Bn)-length(misB)+1
      colnames(Bn)[lenB:ncol(Bn)] = names(misBl)
    }
    
    Bn <- Bn[,colnames(An)]
    
    # final bind
    A = rbind(An, Bn, use.names = T)
    print(c(length(M.list), i))
  } 
  A
}

efg.human.count.merge <- merge.sparse(c(efg.human.count.nodup,efg.human.count.dup.mean.spar))

# Initial counts
efg.seu.lift <- CreateSeuratObject(counts = efg.human.count.merge, meta.data = efg.seu@meta.data)


# data wrangle
# "nCount_RNA"   "nFeature_RNA" "Cluster"      "ClustName"    "Tissue" 
efg.seu.lift@meta.data[,c("orig.ident", "percent.mito","nGene","nUMI","Stages","PreviousClusters","PreviousLineages","AaronLineageAnnotations","cluster_name",
                     "S.Score","G2M.Score","Phase","old.ident","CC.Difference","res.0.4","AaronLineageAnnotations1","CellType"
)] <- list(NULL)
colnames(efg.seu.lift@meta.data) = c('nCount_RNA', 'nFeature_RNA','stim')

fibroX.pdac@meta.data[,c("Cluster", "ClustName")] <- list(NULL)
colnames(fibroX.pdac@meta.data) <- c('nCount_RNA', 'nFeature_RNA','stim')

TS.seu.fibro@meta.data[,c("organ_tissue","method","donor",
                          "anatomical_information","n_counts_UMIs",
                          "n_genes","cell_ontology_class",
                          "free_annotation","manually_annotated",
                          "compartment","gender")] <- list(NULL)
TS.seu.fibro@meta.data$stim <- 'normal'


### merge onevsone
### find keep common genes
efg.pdac.common <- rownames(efg.seu.lift)[rownames(efg.seu.lift) %in% rownames(fibroX.pdac)]
efg.normal.common <- rownames(efg.seu.lift)[rownames(efg.seu.lift) %in% rownames(TS.seu.fibro)]
pdac.normal.common <- rownames(fibroX.pdac)[rownames(fibroX.pdac) %in% rownames(TS.seu.fibro)]

common_genes <- Reduce(intersect, list(efg.pdac.common,efg.normal.common,pdac.normal.common))


# SCT integration
efg.seu.lift <- SCTransform(efg.seu.lift, vst.flavor = "v2", verbose = FALSE)
fibroX.pdac <- SCTransform(fibroX.pdac, vst.flavor = "v2", verbose = FALSE) 
TS.seu.fibro <- SCTransform(TS.seu.fibro, vst.flavor = "v2", verbose = FALSE) 

ienp.list <- list(normal=TS.seu.fibro,pdac=fibroX.pdac, efg = efg.seu.lift)
sctinter <- function(i.list,commone_genes){
  i.features <- SelectIntegrationFeatures(object.list = i.list, nfeatures = 3000)
  i.list <- PrepSCTIntegration(object.list = i.list, anchor.features = i.features)
  
  i.anchors <- FindIntegrationAnchors(object.list = i.list, normalization.method = "SCT",
                                        anchor.features = i.features,k.filter=20)
  
  i.combined.sct <- IntegrateData(anchorset = i.anchors, 
                                  normalization.method = "SCT",
                                  features.to.integrate = common_genes,k.weight=0)
  return(i.combined.sct)
}

ienp.combined.sct <- sctinter(ienp.list,common_genes)
Idents(ienp.combined.sct)<- ienp.combined.sct@meta.data$stim
ienp.combined.sct <- PrepSCTFindMarkers(ienp.combined.sct)


i.efg.pdac.markers <- FindMarkers(ienp.combined.sct,assay = "SCT",ident.1 = 'Spl_Meso', ident.2='PDAC',verbose = FALSE)
i.efg.normal.markers <- FindMarkers(ienp.combined.sct,assay = "SCT",ident.1 = 'Spl_Meso', ident.2 = 'normal',verbose = FALSE)
i.pdac.normal.markers <- FindMarkers(ienp.combined.sct,assay = "SCT",ident.1 = 'PDAC',ident.2='normal',verbose = FALSE)

Eup.Ndw.Pup <- Reduce(intersect, list(rownames(i.efg.normal.markers[i.efg.normal.markers$avg_log2FC>0,]),
                                      rownames(i.pdac.normal.markers[i.pdac.normal.markers$avg_log2FC>0,])))
Eup.Nup.Pdw <- Reduce(intersect, list(rownames(i.efg.pdac.markers[i.efg.pdac.markers$avg_log2FC>0,]),
                                      rownames(i.pdac.normal.markers[i.pdac.normal.markers$avg_log2FC<0,])))

ave.exp.sct <- AverageExpression(ienp.combined.sct,assays = "SCT",slot='data')$SCT
plotAveExpHeatmap2 <- function(ave_exp, features, rowsplit,clusterrows,barname){
  library(ComplexHeatmap)
  ave_exp_sub <- as.data.frame(subset(ave_exp[features,]))
  ave_exp_sub <- as.matrix(ave_exp_sub[c('Spl_Meso','normal','PDAC')])
  my.colors <- colorRamp2(c(-1.5,0,1.5),c("skyblue","black", "coral1"))
  ave_exp_sub_z <- t(scale(t(ave_exp_sub)))
  Heatmap(ave_exp_sub_z, name = barname, col = my.colors, cluster_columns =FALSE,
          show_row_names = FALSE,cluster_rows = clusterrows,
          row_split = rowsplit,
          show_parent_dend_line = FALSE,row_title=NULL,
          row_names_gp = gpar(fontsize = 3),cluster_row_slices = FALSE
  )
}

Eup.Ndw.Pup.sct.exp <- as.data.frame(ave.exp.sct[Eup.Ndw.Pup,])
Eup.Nup.Pdw.sct.exp <- as.data.frame(ave.exp.sct[Eup.Nup.Pdw,])

Eup.Ndw.Pup.sct.exp <- cbind(Eup.Ndw.Pup.sct.exp,i.efg.normal.markers[rownames(Eup.Ndw.Pup.sct.exp),])
Eup.Nup.Pdw.sct.exp <- cbind(Eup.Nup.Pdw.sct.exp,i.efg.pdac.markers[rownames(Eup.Nup.Pdw.sct.exp),])

Eup.Ndw.Pup.sct.FC.srt <- rownames(Eup.Ndw.Pup.sct.exp[order(Eup.Ndw.Pup.sct.exp$avg_log2FC,decreasing = TRUE),])
Eup.Nup.Pdw.sct.FC.srt <- rownames(Eup.Nup.Pdw.sct.exp[order(Eup.Nup.Pdw.sct.exp$avg_log2FC,decreasing = TRUE),])

features.selected <- c(Eup.Ndw.Pup.sct.FC.srt,Eup.Nup.Pdw.sct.FC.srt)
rowsplit.selected <- c(rep("normal",length(Eup.Ndw.Pup.sct.FC.srt)),
                       rep("PDAC",length(Eup.Nup.Pdw.sct.FC.srt)))

plotAveExpHeatmap2(ave.exp.sct,features.selected,
                   rowsplit.selected,FALSE,"z(mat)")

sheets_filt <- list("Eup.Ndw.Pup" = cbind(genes=rownames(Eup.Ndw.Pup.sct.exp),Eup.Ndw.Pup.sct.exp[,c("normal","PDAC","Spl_Meso")]), 
                    "Eup.Nup.Pdw" = cbind(genes=rownames(Eup.Nup.Pdw.sct.exp),Eup.Nup.Pdw.sct.exp[,c("normal","PDAC","Spl_Meso")]))
write_xlsx(sheets_filt, 'Spl_Meos_human.PDAC.Normal.xlsx')