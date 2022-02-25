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

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig3/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/subTreg/nPC_30/UMAP_dist_0.01_nneighbor_30/p2Treg_UMAP_dist_0.01_nneighbor_30_CLUSTER_res_0.1/cluster.rds"

seuratObj <- readRDS(dataPath)
seuratObj@meta.data$newClusterID <- 1
seuratObj@meta.data$newClusterID[seuratObj@meta.data$seurat_clusters == 0] <- 2

Idents(seuratObj) <- seuratObj@meta.data$orig.ident
seuratObj <- subset(seuratObj, idents = c("A13", "A50", "A13B", "A50B"), invert = T)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
detach("package:plyr", unload=TRUE)
library(dplyr)

md <- as_tibble(seuratObj@meta.data)
totalCellPerSample <- md %>%
    group_by(control, orig.ident) %>%
    summarise(cells = n())
seuratObjMD <- md %>%
    group_by(newClusterID , GROUP, control, orig.ident) %>%
    summarise(cells = n())
seuratObjMD$freq <- seuratObjMD$cells / totalCellPerSample$cells[match(seuratObjMD$orig.ident, totalCellPerSample$orig.ident)]
allStatFreq <- seuratObjMD %>%
    filter(GROUP != "OA")
allStat <- allStatFreq %>%
    group_by(newClusterID, control) %>%
    summarise(avg.freq = mean(freq),
              sd.freq = sd(freq),
              se.freq = sd(freq)/sqrt(n()))
allStat_g3 <- as_tibble(allStat)
allStat_g3$newClusterID <- as.factor(allStat_g3$newClusterID)
allStatFreq$newClusterID <- as.factor(allStatFreq$newClusterID)

g <- ggplot(data = allStat_g3, mapping = aes(x = newClusterID, y = avg.freq, color = control, fill = control)) +
    theme_bw() +
    geom_boxplot() +
    geom_bar(stat = "identity", position=position_dodge()) +
    geom_errorbar(data = allStat_g3,
                  aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),
                  width=.8,
                  size = 0.2,
                  color = "black",
                  position=position_dodge(.9)) +
    geom_jitter(data = allStatFreq,
                mapping = aes(x = newClusterID, y = freq, fill = control),
                color = "black",
                position = position_dodge(width = 0.9),
                size = 0.2) +
    ## scale_y_continuous(expand=c(0,0),limits=c(0,0.3)) +
    ## theme_classic() +
    theme(legend.position = 'none')
filePath <- file.path(figurePath, paste0("figure3c.pdf"))
ggsave(filePath, g, height=4, width=4)

write_tsv(allStat_g3, file.path(figurePath, paste0('figure3c.tsv')))
