## Setup ----
require(Hmisc) # for %nin%
require(Seurat)
require(magrittr)
require(ggplot2)
require(RColorBrewer)
require(tidyverse)
require(fgsea)
require(viridis)
require(ComplexHeatmap)

## ----------------------- -
## UMAPs and tdTomato DimPlot
# load('seurat_anaylsis_combined_clusters_210827.RData') # original dataset; 8 samples
load('new_combined_cluster_8samples_2022-01-14.RData') # re-processed; rescued cells
# load('new_combined_cluster_FibroSamples_2022-01-14.RData') # subsetted to only the 4 fibroblast samples

DefaultAssay(combined_cluster) <- "RNA"
DimPlot(combined_cluster, reduction = "umap", group.by="stim", label = FALSE)
DimPlot(combined_cluster, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE, repel = TRUE)

## Show tdTomato expression
opts <- list( scale_colour_gradientn(colours = c(rep(rgb(.7,.7,.7,.3), 1),
                                                 rep('firebrick'     , 3)) ),
              theme(axis.title = element_text(size=8, face="plain"),
                    plot.title = element_text(size=8, face="plain")) )

# FeaturePlot(combined_cluster, features = "tdTomato", order = TRUE, pt.size=0.3) & opts # for 2021-08-27 dataset
FeaturePlot(combined_cluster, features = "Raw_tdTomato.v2", order = TRUE, pt.size=0.3) & opts # for 2022-01-14 dataset

## Show simplified sample coloring
u <- combined_cluster@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(sample = combined_cluster@meta.data$stim) %>%
  separate(., "sample", "_", into = c("type","orig","time"), remove = FALSE)

han_gray <- '#5b5c5e'
han_red  <- '#ec2028'
sample_colors <- c("epithelial_normal_day4rep1" = han_gray, "epithelial_normal_day4rep2" = han_gray,
                   "epithelial_PDAC_day2"       = han_gray, "epithelial_PDAC_day3"       = han_gray,
                   "fibroblast_normal_day4rep1" = han_red , "fibroblast_normal_day4rep2" = han_red,
                   "fibroblast_PDAC_day2"       = han_red , "fibroblast_PDAC_day3"       = han_red  )
ggplot(u, aes(x=UMAP_1, y=UMAP_2, color=sample)) + geom_point(size=0.3, alpha=1) +
  scale_color_manual(values=sample_colors) +
  theme_bw() + coord_fixed(ratio = 1) + facet_grid(rows=vars(type), cols=vars(orig))


## ----------------------- -
## Combine into one fibroblast object -

##
# load('new_combined_cluster_FibroSamples_2022-01-14.RData') # from GSE200903


## 
# load('merged_fibroblast_dataset.v4.RData') # concordant wiht the 2022-01-14 dataset, but re-integrated with fibro only
# merged.f2$CellGroup <- plyr::mapvalues(merged.f2$CellGroup,
#                                        from = c("fibroblast_PDAC.fibroblasts", "fibroblast_normal.fibroblasts"),
#                                        to   = c('PDAC Fibroblasts'           , 'Normal Fibroblasts'))

##
load('seurat_anaylsis_combined_clusters_210827.RData') # original dataset; 8 samples
merged.f2 <- subset(combined_cluster,
                    subset=(cluster_name %in% c('inflammatory fibroblasts','myofibroblasts')) &
                           (stim         %in% c('fibroblast_normal_day4rep1','fibroblast_normal_day4rep2',
                                                'fibroblast_PDAC_day2', 'fibroblast_PDAC_day3')) )
merged.f2$CellGroup <- plyr::mapvalues(merged.f2$stim,
                                       from = c("fibroblast_PDAC_day2", 'fibroblast_PDAC_day3',
                                                "fibroblast_normal_day4rep1", "fibroblast_normal_day4rep2"),
                                       to   = c('PDAC Fibroblasts'    ,'PDAC Fibroblasts'     ,
                                                'Normal Fibroblasts'        , 'Normal Fibroblasts'))


## ----------------------- -
## Lu's embryonic data -
E9_markers <- read_tsv("E9.5_embryo/Markers_All_Clusters.txt",
                       col_names = c('rowname','p_val','avg_log2FC','pct.1','pct.2',
                                     'p_val_adj','cluster','gene'), skip = 1)
head(E9_markers)
# tapply(E9_markers$gene, E9_markers$cluster, function(x){vt(x, rownames(fibro))})

## read sc data
TotalCells_E9.5_Updated <- readRDS("E9.5_embryo/TotalCells_E9.5_Updated.rds")

lu.meso <- subset(TotalCells_E9.5_Updated, subset=LineageAnnotations %in% c('Cardiac','Other_Meso','Spl_Meso'))
lu.meso@meta.data$LineageAnnotations <- lu.meso@meta.data$LineageAnnotations %>% droplevels
table(lu.meso@meta.data$LineageAnnotations)

# message('intersection of genes from previous study and the current study:')
# vt(rownames(lu.meso), rownames(fibro))
# g <- intersect(rownames(lu.meso), rownames(fibro))
# str(g)


## ------------------------------------------------- - 
## Tabula Sapiens pancreas fibroblasts/stellate ----
load('TabulaSapiens/tabula_sapiens_seurat.panc.fibro.RData') # ts.panc.fibro

table(ts.panc.fibro$free_annotation %>% droplevels) %>% sort %>% rev

ts.panc.fibro$free_annotation2 <- ts.panc.fibro$free_annotation %>% 
  plyr::mapvalues(., from = 'fibroblast', to = 'pancreatic stellate cell') %>% droplevels
table(ts.panc.fibro$free_annotation2) %>% sort %>% rev

ts <- NormalizeData(ts.panc.fibro)
ts.a <- AverageExpression(ts)
ts.s <- AverageExpression(ts, assays="RNA", group.by = "free_annotation")$RNA


## ------------------------------------------------- - 
## DE Genes ----
a   <- AverageExpression(merged.f2, group.by = "CellGroup")
d.r <- dist(t(a$RNA))
dst <- a$RNA[, "Normal Fibroblasts"] - a$RNA[, "PDAC Fibroblasts"]

deg <- FindMarkers(merged.f2, group.by = "CellGroup",
                   ident.1 = "Normal Fibroblasts", ident.2 = "PDAC Fibroblasts",
                   test.use = "MAST")

g.m <- homologene::homologene(rownames(deg), inTax = 10090, outTax = 9606)
m   <- match(rownames(deg), g.m$`10090`)
all(rownames(deg) == g.m$`10090`[m], na.rm = TRUE)
deg$human_symbol <- g.m$`9606`[m]

m <- match(rownames(deg), names(dst))
all(names(dst)[m] == rownames(deg), na.rm = TRUE)

smoothScatter(pseudo_log2(dst[m]), deg$avg_log2FC, nrpoints = 0,
              main = sprintf("Cor(diff, DE FC) = %.2f", cor(pseudo_log2(dst[m]), deg$avg_log2FC)) )
abline(lm(deg$avg_log2FC ~ pseudo_log2(dst[m])), col="gray")


## DE genes - Heatmap ----
col_rat = circlize::colorRamp2(c(-5          , 0      , 5    ),
                               c("dodgerblue", "white", "darkred"))
col_08  = circlize::colorRamp2(c(-5          , 0      , 7    ),
                               c("dodgerblue", "white", "darkred"))
col_02  = circlize::colorRamp2(c(-2          , 0      , 2    ),
                               c("purple4"   , "white", "darkorange"))
col_z   = circlize::colorRamp2(c(-2          , 0      , 2    ),
                               c("dodgerblue", "white", "darkred"))

## FC heatmap, 2 columns sorted by PDAC level
mat.2 <- as.matrix(a$RNA[ rownames(deg), ])
# mat.2 <- as.matrix(a$RNA[ rownames(deg), c(4,5) ])
all(rownames(mat.2) == rownames(deg))

Heatmap(mat.2 %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2),
        show_row_names = FALSE, name = "log2\nGene\nExpression", column_title = "DE Genes, ordered by PDAC level",
        row_order  = order(mat.2[,2]), col = col_rat, column_order = c(2,1), border = TRUE)

## FC heatmap, 2 columns, FC order, with FC shown
Heatmap(mat.2 %>% apply(., 2, log2),
        show_row_names = FALSE, name = "log2\nGene\nExpression", column_title = "DE Genes, ordered by FC",
        row_order  = order(deg$avg_log2FC), col = col_08, column_order = c(2,1), border = TRUE) +
Heatmap(deg$avg_log2FC %>% matrix(., ncol = 1) %>% set_colnames(., c("log2FC")),
        show_row_names = FALSE, name = "log2FC",
        row_order  = order(deg$avg_log2FC), col = col_02, border = TRUE)


## Same, with Tabula Sapiens data added
# ts.a$RNA[deg$human_symbol,] %>% is.na %>% table
ts.m <- match( deg$human_symbol, rownames(ts.a$RNA) )
all(deg$human_symbol == rownames(ts.a$RNA)[ts.m], na.rm = TRUE)

Heatmap(mat.2 %>% apply(., 2, log2),
        show_row_names = FALSE, name = "log2\nGene\nExpression", column_title = "DE Genes",
        row_order  = order(deg$avg_log2FC), col = col_08, column_order = c(2,1), border = TRUE) +
Heatmap(deg$avg_log2FC %>% matrix(., ncol = 1) %>% set_colnames(., c("log2FC")),
        show_row_names = FALSE, name = "log2FC",
        row_order  = order(deg$avg_log2FC), col = col_02, border = TRUE) +
Heatmap(ts.a$RNA[ts.m] %>% matrix(., ncol = 1) %>% log2 %>% set_colnames(., c("Tabula Sapiens\nPancreas\nFibroblasts")),
        show_row_names = FALSE, name = "log2\nGene\nExpression",
        row_order  = order(deg$avg_log2FC), col = col_08, border = TRUE)

cor( cbind(mat.2[,1] %>% log2, # normal
           mat.2[,2] %>% log2, # pdac
           ts.a$RNA[ts.m] %>% log2), use="complete.obs", method = "spearman")

obs <- cor( mat.2[,2], ts.a$RNA[ts.m], use="complete.obs", method = "spearman")
bg  <- sapply(1:1000, function(i){
  ri <- sample(1:dim(ts.a$RNA)[1], length(ts.m), replace = FALSE)
  cor( mat.2[,2], ts.a$RNA[ri], use="complete.obs", method = "spearman")
})
plot(density(bg), xlim=range(c(bg,obs)))
abline(v=obs, col="orange", lwd=2)


## FibroExplorer data ----
fx.hPS <- readRDS('fibroexplorer/Human_PS_Fibro.RDS')
fx.mPS <- readRDS('fibroexplorer/Mouse_PS_Fibro.RDS')
fx.mSS <- readRDS('fibroexplorer/Mouse_SS_Fibro.RDS')

table(fx.hPS$ClustName, fx.hPS$Tissue)
fx.hPDAC <- subset(fx.hPS, subset=(Tissue  %in% c('PDAC')))
fx.mPDAC <- subset(fx.mPS, subset=(Tissue  %in% c('Large_Tumor_PDAC','Small_Tumor_PDAC')))
fx.mPanc <- subset(fx.mSS, subset=(Tissue  %in% c('Pancreas')))
fx.mOther<- subset(fx.mSS, subset=(Tissue %nin% c('Pancreas')))

fx.hc <- AverageExpression(NormalizeData(fx.hPDAC), assays="RNA", group.by = 'Tissue')$RNA
fx.mc <- AverageExpression(NormalizeData(fx.mPDAC), assays="RNA", group.by = 'Tissue')$RNA
fx.mn <- AverageExpression(NormalizeData(fx.mPanc), assays="RNA", group.by = 'Tissue')$RNA

fx.m <- base::merge(fx.mc, fx.mn, by=0, all=TRUE)
colnames(fx.m)[2] <- "Large-Tumor PDAC Fibroblasts"
colnames(fx.m)[3] <- "Small-Tumor PDAC Fibroblasts"
colnames(fx.m)[4] <- "Normal Fibroblasts"
colnames(fx.hc)   <- "PDAC Fibroblasts"

##
fxts <- base::merge(fx.hc, ts.a$RNA, by=0, all=TRUE)
colnames(fxts) <- c("gene","fx","ts")
fxts.m <- homologene::homologene(fxts$gene, inTax = 9606, outTax = 10090)

##
fx.i <- match( rownames(deg), fx.m$Row.names )
all(rownames(deg) == fx.m$Row.names[fx.i], na.rm = TRUE)

fx.hi <- match( deg$human_symbol, rownames(fx.hc) )
all(deg$human_symbol == rownames(fx.hc)[fx.hi], na.rm = TRUE)

Heatmap(mat.2 %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = "DE Genes",
        row_order  = order(deg$avg_log2FC), col = col_z, column_order = c(2,1), border = TRUE) +
Heatmap(deg$avg_log2FC %>% matrix(., ncol = 1) %>% set_colnames(., c("FC")),
        show_row_names = FALSE, name = "log2FC",
        row_order  = order(deg$avg_log2FC), col = col_02, border = TRUE) +
Heatmap(fx.m[fx.i, 2:4] %>% as.matrix %>% plyr::mapvalues(., c(NA,0), rep(1e-6,2)) %>% apply(., c(1,2), log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = "FibroXplorer, Mouse",
        row_order  = order(deg$avg_log2FC), col = col_z, column_order = c(1,2,3), border = TRUE)


##
Heatmap(mat.2 %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = "DE Genes",
        row_order  = order(deg$avg_log2FC), col = col_z, column_order = c(2,1), border = TRUE) +
Heatmap(deg$avg_log2FC %>% matrix(., ncol = 1) %>% set_colnames(., c("FC")),
        show_row_names = FALSE, name = "log2FC",
        row_order  = order(deg$avg_log2FC), col = col_02, border = TRUE) +
Heatmap(fx.m[fx.i, 2:4] %>% as.matrix %>% plyr::mapvalues(., c(NA,0), rep(1e-6,2)) %>% apply(., c(1,2), log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = "FX, Mouse",
        row_order  = order(deg$avg_log2FC), col = col_z, column_order = c(1,2,3), border = TRUE) +
Heatmap(fx.hc[fx.hi, 1] %>% as.matrix %>% plyr::mapvalues(., c(NA,0), rep(1e-6,2)) %>% apply(., c(1,2), log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = "FX, Human",
        row_order  = order(deg$avg_log2FC), col = col_z, border = TRUE)

cor( cbind(mat.2,
           fx.m[fx.i, 2:4] %>% as.matrix,
           fx.hc[fx.hi, 1] %>% as.matrix
           ), use="complete.obs", method = "spearman") %>% round(.,2)

cor( cbind(a$RNA[,4] %>% log2, # normal
           a$RNA[,5] %>% log2, # pdac
           fx.m[ match( rownames(a$RNA), fx.m$Row.names ), 2:4] %>% as.matrix
           ), use="complete.obs", method = "spearman") %>% round(.,2)

##
##
e9i <- 'Spl_Meso'
E9_markers$gene[E9_markers$cluster == e9i] %in% fxts.m$`10090` %>% table

fxts.e9i <- fxts[ fxts$gene %in% fxts.m$`9606`[fxts.m$`10090` %in% E9_markers$gene[E9_markers$cluster == e9i]], ]
i <- match(fxts.e9i$gene, fxts.m$`9606`)
all( fxts.m$`9606`[i] == fxts.e9i$gene )
fxts.e9i$gene.mouse <- fxts.m$`10090`[i]

fxts.e9i$gene.mouse %in% rownames(a$RNA) %>% table
fxts.e9i <- fxts.e9i[ fxts.e9i$gene.mouse %in% rownames(a$RNA), ]
colnames(fxts.e9i)[2:3] <- c("Tabula Sapiens\nNormal Fibroblasts", "FibroXplorer\nhPDAC")
# str(fxts.e9i)

roh <- (fxts.e9i[, 3] / fxts.e9i[, 2]) %>% plyr::mapvalues(., c(NaN, NA, Inf, -Inf), c(0, 0, 8, -8)) %>% log2

##
mat.ts <- as.matrix(a$RNA[ fxts.e9i$gene.mouse, ])
# mat.ts <- as.matrix(a$RNA[ fxts.e9i$gene.mouse, c(4,5) ])
mat.ts.rat <- log2(mat.ts[, 2] / mat.ts[, 1]) %>% matrix(., ncol = 1) %>% set_colnames(., c("log2ratio"))

##
ri <- order(mat.ts.rat, decreasing = TRUE)
Heatmap(mat.ts %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = str_interp("${e9i} Genes"),
        row_order  = ri, col = col_z, column_order = c(2,1), border = TRUE) +
  Heatmap(mat.ts.rat,
          show_row_names = FALSE, name = "log2ratio\n(PDAC / Normal)",
          row_order  = ri, col = col_02, border = TRUE) +
  Heatmap(roh,
          show_row_names = FALSE, name = "log2ratio\n(TabulaSapiens\n / FibroXplorer)",
          row_order  = ri, col = col_02, border = TRUE) +
  Heatmap(fxts.e9i[, 2:3] %>% as.matrix %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
          show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = str_interp("${e9i} Genes"),
          row_order  = ri, col = col_z, column_order = c(2,1), border = TRUE)

cor( cbind(mat.ts, 
           fxts.e9i[, 2:3]
           ), use="complete.obs", method = "spearman") %>% round(.,2) %>% print


##
ri <- order(roh, decreasing = TRUE)
Heatmap(mat.ts %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
        show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = str_interp("${e9i} Genes"),
        row_order  = ri, col = col_z, column_order = c(2,1), border = TRUE) +
  Heatmap(mat.ts.rat,
          show_row_names = FALSE, name = "log2ratio\n(PDAC / Normal)",
          row_order  = ri, col = col_02, border = TRUE) +
  Heatmap(roh,
          show_row_names = FALSE, name = "log2ratio\n(TabulaSapiens\n / FibroXplorer)",
          row_order  = ri, col = col_02, border = TRUE) +
  Heatmap(fxts.e9i[, 2:3] %>% as.matrix %>% plyr::mapvalues(., 0, 1e-6) %>% apply(., 2, log2) %>% scale,
          show_row_names = FALSE, name = "Z(log2 Gene\nExpression)", column_title = str_interp("${e9i} Genes"),
          row_order  = ri, col = col_z, column_order = c(2,1), border = TRUE)

