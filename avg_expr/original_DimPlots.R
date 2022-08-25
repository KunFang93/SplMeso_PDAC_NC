require(Seurat)
require(magrittr)
require(ggplot2)
require(RColorBrewer)
require(tidyverse)

## Make 4-condition figure ----
load('new_combined_cluster_8samples_2022-01-14.RData')

u <- combined_cluster@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(sample = combined_cluster@meta.data$stim) %>%
  separate(., "sample", "_", into = c("type","orig","time"), remove = FALSE)

han_gray <- '#5b5c5e'
han_red  <- '#ec2028'
sample_colors <- c("epithelial_normal_day4rep1" = han_gray, "epithelial_normal_day4rep2" = han_gray,
                   "epithelial_PDAC_day2"       = han_gray, "epithelial_PDAC_day3"       = han_gray,
                   "fibroblast_normal_day4rep1" = han_red , "fibroblast_normal_day4rep2" = han_red,
                   "fibroblast_PDAC_day2"       = han_red , "fibroblast_PDAC_day3"       = han_red  )

# pdf('scRNA_per-sample-panels.v3.pdf', height=6, width=8, onefile = TRUE)
ggplot(u, aes(x=UMAP_1, y=UMAP_2, color=sample)) + geom_point(size=0.3, alpha=1) +
  scale_color_manual(values=sample_colors) +
  theme_bw() + coord_fixed(ratio = 1) + facet_grid(rows=vars(type), cols=vars(orig)) + NoLegend()
# dev.off()


## DimPlot by cell type ----
load('new_combined_cluster_FibroSamples_2022-01-14.RData')

DimPlot(combined_cluster, reduction = "umap", group.by = "stim", label = TRUE, repel = TRUE)

DimPlot(combined_cluster, split.by = "stim", group.by="CellType",
        cols = c('black', han_gray, han_red, 'orange', 'royalblue', 'violet', 'salmon'), na.value = 'beige')
