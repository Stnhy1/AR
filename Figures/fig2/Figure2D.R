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

md <- as_tibble(seuratObj@meta.data)
md$newGROUP <- md$GROUP
md$newGROUP[md$GROUP %in% c("Mono", "Combo")] <- "Mono+Combo"

totalCellPerSample <- md %>%
    filter(newGROUP != "OA") %>% 
    group_by(control, orig.ident) %>%
    count() %>%
    rename(cells = n)



sampleStat <- md %>%
    group_by(newClusterID, newGROUP, control, orig.ident) %>%
    count() %>%
    rename(cells = n) %>%
    filter(newGROUP == "Mono+Combo") %>%
    mutate(group=paste0(control, "_", newGROUP))

sampleStat$freq <- sampleStat$cells / totalCellPerSample$cells[match(sampleStat$orig.ident, totalCellPerSample$orig.ident)]


allStatFreq <- sampleStat %>%
    filter(newGROUP != "OA")

allStatNew <- sampleStat %>%
    group_by(newClusterID, group) %>%
    dplyr::summarise(avg.freq = mean(freq), sd.freq = sd(freq), se.freq = sd(freq)/sqrt(n()))

allStatNew$newClusterID <- as.factor(allStatNew$newClusterID)
allStat_g3 <- allStatNew %>% filter(group != "OA") %>%  mutate(newGROUP = group)
g <- ggplot(data = allStat_g3,
            mapping = aes(x = newClusterID,
                          y = avg.freq,
                          color = group,
                          fill = group)) +
    theme_bw() +
    geom_boxplot() +
    geom_bar(stat = "identity", position=position_dodge()) +
    geom_errorbar(data = allStat_g3,
                  aes(ymin=avg.freq-sd.freq, ymax=avg.freq+se.freq),
                  width=.8,
                  size = 0.2,
                  color = "black",
                  position=position_dodge(.9))+
    geom_jitter(data = allStatFreq,
                mapping = aes(x = newClusterID, y = freq, fill = group),
                color = "black",
                position = position_dodge(width = 0.9),
                size = 0.2) +
    coord_flip() +
    scale_y_continuous(expand=c(0,0), limits = c(0, 0.35)) +
    theme(legend.position = "none")
filePath <- file.path(figurePath, paste0("Figure2D.pdf"))
ggsave(filePath, g, height=10, width=5)

write_tsv(allStat_g3, file.path(figurePath, paste0('Figure2D.tsv')))
