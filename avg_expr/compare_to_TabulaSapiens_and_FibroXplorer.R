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


## Lu's embryonic data -
E9_markers <- read_tsv("E9.5_embryo/Markers_All_Clusters.txt",
                       col_names = c('rowname','p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene'), skip = 1)
head(E9_markers)
tapply(E9_markers$gene, E9_markers$cluster, function(x){vt(x, rownames(fibro))})

## read sc data
TotalCells_E9.5_Updated <- readRDS("E9.5_embryo/TotalCells_E9.5_Updated.rds")

lu.meso <- subset(TotalCells_E9.5_Updated, subset=LineageAnnotations %in% c('Cardiac','Other_Meso','Spl_Meso'))
lu.meso@meta.data$LineageAnnotations <- lu.meso@meta.data$LineageAnnotations %>% droplevels
table(lu.meso@meta.data$LineageAnnotations)

message('intersection of genes from previous study and the current study:')
vt(rownames(lu.meso), rownames(fibro))
g <- intersect(rownames(lu.meso), rownames(fibro))
str(g)


## ----------------------- -
## Combine into one fibroblast object -
load('merged_fibroblast_dataset.v4.RData') # or, load('new_combined_cluster_FibroSamples_2022-01-14.RData'), from GSE200903

merged.f2$CellGroup <- plyr::mapvalues(merged.f2$CellGroup,
                                       from = c("fibroblast_PDAC.fibroblasts", "fibroblast_normal.fibroblasts"),
                                       to   = c('PDAC Fibroblasts'           , 'Normal Fibroblasts'))

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
mat.2 <- as.matrix(a$RNA[ rownames(deg), c(4,5) ])
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

