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
library(plyr)

dataPath <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/dataForpaper.rds'
figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig1/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

fileName <- "Figure1F"
seuratObj <- readRDS(dataPath)
Idents(seuratObj) <- seuratObj@meta.data$orig.ident
seuratObj <- subset(seuratObj, idents = c("A13", "A50", "A13B", "A50B"), invert = T)
seuratObj@meta.data$newClusterID <- mapvalues(seuratObj@meta.data$newClusterID,
                                              from=c(1:16,18:31),
                                              to=c(1:30))
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
seuratObj <- subset(seuratObj, idents = c(30), invert = T)
## seuratObj@meta.data$newClusterID <- factor(seuratObj@meta.data$newClusterID, levels = 1:29)
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
detach("package:plyr", unload=TRUE)
library(dplyr)

md <- as_tibble(seuratObj@meta.data)
totalCellPerSample <- md %>%
    group_by(control, orig.ident) %>%
    summarise(cells = n())
seuratObjMD <- md %>%
    group_by(newClusterID, GROUP, control, orig.ident) %>%
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
allStat_g3 <- allStat_g3 %>% add_row(newClusterID = 23, control = "SF", avg.freq = 0, sd.freq = 0, se.freq = 0)
allStat_g3 <- allStat_g3 %>% add_row(newClusterID = 24, control = "SF", avg.freq = 0, sd.freq = 0, se.freq = 0)
allStat_g3 <- allStat_g3 %>% add_row(newClusterID = 28, control = "SF", avg.freq = 0, sd.freq = 0, se.freq = 0)
allStat_g3 <- allStat_g3 %>% add_row(newClusterID = 29, control = "SF", avg.freq = 0, sd.freq = 0, se.freq = 0)
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
filePath <- file.path(figurePath, paste0("figure1e.pdf"))
ggsave(filePath, g, height=8, width=10)


write_tsv(allStat_g3, file.path(figurePath, paste0('figure1e.tsv')))


g <- allStatFreq %>%
    ggstatsplot::grouped_ggbetweenstats(
                     data = .,
                     x = control,
                     y = freq,
                     grouping.var = newClusterID,
                     xlab = "Cancer type",
                     ylab = "Sample fraction",
                     ## pairwise.display = "all", # display only significant pairwise comparisons
                     p.adjust.method = "fdr", # adjust p-values for multiple tests using this method
                     ggtheme = theme_classic(),
                     package = "ggsci",
                     palette = "default_jco",
                     plotgrid.args = list(ncol = 3))
ggsave(file.path(figurePath, paste0("figure1e_test.pdf")), g, width = 410, height = 900, units =  "mm")
