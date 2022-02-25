rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
require(plyr)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig2/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

#' whole data seuratObj  ######################################################
dataPath <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/dataForpaper.rds'

seuratObj <- readRDS(dataPath)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID

coord = Embeddings(object = seuratObj, reduction = "umap")
colnames(coord) = c('UMAP1','UMAP2')
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
label <- meta %>%
    group_by(newClusterID) %>%
    select(UMAP1, UMAP2) %>%
    summarize_all(mean)
meta$control <- factor(meta$control, levels = c("Blood", "SF"))
g <- ggplot() +
    geom_point(data = meta,
               mapping = aes(x = UMAP1, y = UMAP2, color = control),
               size = 0.1) +
    theme_void() +
    theme(legend.position = "none")
png(file.path(figurePath, "Figure2C.png"), width = 210/2, height = 210/2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
print(g)
dev.off()
