rm(list=ls())
library(Seurat)
library(patchwork)
library(readr)
library(tidyverse)
library(rjson)
library(rlist)
library(SeuratData)
library(harmony)
library(Rmagic)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)
require(plyr)

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
        xlab("") + ylab(feature) + ggtitle("") +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.1),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.margin = plot.margin )
    return(p)
}
extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
}
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
        theme(axis.text.x=element_text(), axis.ticks.x = element_line())
    ymaxs<- purrr::map_dbl(plot_list, extract_max)
    plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
}

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig2/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

dataPath <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/dataForpaper.rds'

seuratObj<- readRDS(dataPath)


features <- c(
    "CD3D",
    "CD4",
    "CCR7",
    "CD69",
    "CXCR3",
    "CXCR5",
    "CXCR6",
    "CCR4",
    "CCR6",
    "KLRB1",
    "IL2RA",
    "FOXP3",
    "CD8A",
    "IL7R",
    "KLRG1",
    "CXCL13",
    "TOX",
    "LAG3",
    "PDCD1",
    "HAVCR2",
    "TIGIT",
    "EOMES",
    "CX3CR1",
    "MKI67",
    "TRGV9",
    "TRGV2",
    "FCGR3A",
    "NCAM1"
)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
levels(seuratObj) <- 1:length(levels(seuratObj))
filePath <- file.path(figurePath, paste0("Figure2B.pdf"))
pdf(filePath, height = 6, width = 6)
StackedVlnPlot(obj = seuratObj, features = features)
dev.off()
