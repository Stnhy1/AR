rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
require(plyr)

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig1")
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}
setwd(figurePath)


dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/dataForpaper.rds"
seuratObj <- readRDS(dataPath)
Idents(seuratObj) <- seuratObj@meta.data$orig.ident
seuratObj <- subset(seuratObj, idents = c("A13", "A50", "A13B", "A50B"), invert = T)
seuratObj@meta.data$newClusterID <- mapvalues(seuratObj@meta.data$newClusterID,
                                              from=c(1:16,18:31),
                                              to=c(1:30))
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
seuratObj <- subset(seuratObj, idents = c(30), invert = T)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID

png(file.path(figurePath, "Figure1B.png"), width = 210/2, height = 210/2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
DimPlot(seuratObj, label = F) +
    theme_void() +
    theme(legend.position = "none")
dev.off()

