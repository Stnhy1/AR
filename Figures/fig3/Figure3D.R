rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/subTreg/nPC_30/UMAP_dist_0.01_nneighbor_30/p2Treg_UMAP_dist_0.01_nneighbor_30_CLUSTER_res_0.1/cluster.rds"

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig3/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

seuratObj <- readRDS(dataPath)

Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters
DEGs <- FindAllMarkers(object = seuratObj,
                       logfc.threshold = 0,
                          only.pos = F)

features <- DEGs %>% filter(gene %in% rownames(seuratObj@assays$RNA@scale.data)) %>%  group_by(cluster) %>% top_n(n = 50, wt=avg_logFC)
pdf(file.path(dirname(dataPath), "heatmap_features.pdf"))
DoHeatmap(seuratObj, features = features$gene)
dev.off()


## DEGs <- read_tsv(file.path(dirname(dataPath), "snn-markers.tsv"))
DEGs <- DEGs[DEGs$cluster == 1,]
DEGs$avg_logFC <- as.numeric(DEGs$avg_logFC)
write_tsv(DEGs, file.path(figurePath, paste0('Figure3D.tsv')))

pdf(file.path(figurePath, "Figure3D.pdf"))
EnhancedVolcano(DEGs,
                lab = str_remove_all(rownames(DEGs), ".1$"),
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-1.5, 1.5),
                title = "Cluster 1",
                pointSize = 2,
                subtitle = NULL,
                shape = c(1, 4, 23, 25),
                pCutoff = 10e-16,
                FCcutoff = 0.5)
dev.off()
