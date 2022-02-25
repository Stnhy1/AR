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

pdf(file.path(figurePath, "Figure2A.pdf"))
DimPlot(seuratObj, label = F) + theme(legend.position = 'none')
dev.off()

png(file.path(figurePath, "Figure2A.png"), width = 210/2, height = 210/2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
DimPlot(seuratObj, label = F) +
    theme_void() +
    theme(legend.position = "none")
dev.off()

