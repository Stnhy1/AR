rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig3/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

#' whole data seuratObj  ######################################################
dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/subTreg/nPC_30/UMAP_dist_0.01_nneighbor_30/p2Treg_UMAP_dist_0.01_nneighbor_30_CLUSTER_res_0.1/cluster.rds"

seuratObj <- readRDS(dataPath)

seuratObj@meta.data$newClusterID <- 1
seuratObj@meta.data$newClusterID[seuratObj@meta.data$seurat_clusters == 0] <- 2

coord = Embeddings(object = seuratObj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP_1", "UMAP_2")
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')

labelT <-  as_tibble(meta) %>%
    select(newClusterID, UMAP_1, UMAP_2) %>%
    group_by(newClusterID) %>%
    summarize(labelx = mean(UMAP_1), labely = mean(UMAP_2))


g <- ggplot(data = meta) +
    geom_point(mapping = aes(x = UMAP_1, y = UMAP_2,
                             color = newClusterID),
               size = 0.5)
pdf(file.path(figurePath, "Figure3A.pdf"))
print(g + theme_classic() + theme(legend.position = 'none'))
dev.off()


g <- ggplot(data = meta) +
    geom_point(mapping = aes(x = UMAP_1, y = UMAP_2, color = control), size = 0.5)
pdf(file.path(figurePath, "Figure3B.pdf"))
print(g + theme_classic() + theme(legend.position = 'none'))
dev.off()
