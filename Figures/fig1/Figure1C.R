###############################################################################
#'                        stacked violine plot slide 14                      '#
###############################################################################

rm(list=ls())
library(Seurat)
library(patchwork)
library(tidyverse)
library(rjson)
library(rlist)
library(Rmagic)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)

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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig1/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
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
seuratObj@meta.data$newClusterID <- factor(seuratObj@meta.data$newClusterID, levels = 1:29)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID

features <- c(
    "PTPRC",
    "CD3D",
    "MKI67",
    "CD4",
    "IL2RA",
    "FOXP3",
    "CD8A",
    "TRDV2",
    "CCR7",
    "IL7R",
    "KLRG1",
    "CX3CR1",
    "NCAM1",
    "CD19",
    "MS4A1",
    "CD27",
    "CD38",
    "SDC1",
    "CD14",
    "FCGR3A",
    "S100A8",
    "S100A9",
    "HLA-DRB1",
    "CD1C",
    "PPBP")
## ,
##     "PLA2G2A")
## levels(seuratObj) <- 1:29
filePath <- file.path(figurePath, paste0("Figure1C.pdf"))
pdf(filePath, height = 6, width = 8)
StackedVlnPlot(obj = seuratObj, features = features)
dev.off()

