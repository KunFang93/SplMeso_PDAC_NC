---

This code describes how we compared data from the current study to data from previous studies.

---

TabulaSapiens data was downloaded and converted into Seurat objects

    SeuratDisk::Convert("TS_Pancreas.h5ad", "TS_Pancreas.h5seurat")
    ts.panc <- LoadH5Seurat("TS_Pancreas.h5seurat", misc = FALSE)
    ts.panc.fibro <- subset(ts.panc, subset=(cell_ontology_class %in% c("pancreatic stellate cell", "fibroblast")))
    save("ts.panc.fibro", file = 'tabula_sapiens_seurat.panc.fibro.RData')

The above RData file was then used herein.

---

DE genes according to the MAST test and ordered by fold-change (FC), shown as standardized (Z) log2(RPKM)
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/avg_expr/EDF7e.png" width="300">

---

This is the spl_meso (splanchnic mesenchyme) geneset, ordered by the log2ratio in the current study (PDAC/normal), and shown as normalized gene expression level; Z(log2(RPKM)). On the right hand size, is FibroXplorer and TabulaSapiens data. TabulaSapiens dataset was converted from hd5 into R objects using SeuratDisk package. Mouse/human gene mappings were done using the homologene R package.
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/avg_expr/EDF7j.png" width="300">

---

Plots that include more details as to the data sources are produced by this script, and will resemble:
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/avg_expr/EDF7e_full.png" width="400">
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/avg_expr/EDF7j_full.png" width="400">

