rm(list=ls())
library(Seurat)
library(patchwork)
library(readr)
library(tidyverse)
library(rjson)
library(rlist)
library(SeuratData)
library(harmony)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig2/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

dataPath <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/dataForpaper.rds'

seuratObj <- readRDS(dataPath)

seuratObj$newClusterID <- factor(seuratObj$newClusterID, levels = 17:1)
Idents(seuratObj) <- seuratObj$newClusterID
features <- c("CD3D", "CD4", "CD8A", "IFNG", "IL4", "IL10",
              "IL17A", "IL21", "TNFA", "GZMB", "GZMM", "GZMH", "GZMK",
              "GNLY", "CSF2", "PRF1")
pdf(file.path(figurePath, "Figure2E.pdf"))
DotPlot(seuratObj, features = features) +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
