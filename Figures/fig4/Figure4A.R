###############################################################################
#'                    TCR shared analysis correct version                    '#
###############################################################################

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
library('tcR')

figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig4/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

dataPath <- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/rmDoublets_noTRMTRSgenes/nPC_30/UMAP_dist_0.01/CLUSTER_res_0.8/dataForpaper.rds'
seuratObj <- readRDS(dataPath)

path<- '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/TCR_CC'
print(paste0("Loading data from ", path))
samples = list.dirs(path = path, recursive = F)
file.id = basename(samples)
sample.id <- file.id
if(length(file.id) < 1) stop(paste0('Failed to find data in ', path))
`%notin%` <- Negate(`%in%`)
tab=c()
for(i in 1:length(file.id)){
    fcap_i <-  file.path(path, file.id[i], 'outs/filtered_contig_annotations.csv')
    if(!file.exists(fcap_i)){next}
    ctp_i <-  file.path(path, file.id[i], 'outs/clonotypes.csv')
    if(!file.exists(ctp_i)){next}

    fca_i <- read.csv(fcap_i) %>% select(barcode, raw_clonotype_id)
    ct_i <- read.csv(ctp_i) %>% select(clonotype_id, cdr3s_nt)
    fca_i$cdr3s_nt <-  ""
    fca_i$cdr3s_nt <- ct_i$cdr3s_nt[match(fca_i$raw_clonotype_id, ct_i$clonotype_id)]
    fca_i$barcode <- paste0(file.id[i], "_", str_remove_all(fca_i$barcode, "-[0-9]+$"))

    fca_i <- fca_i %>% select(barcode, cdr3s_nt)
    tab=rbind(tab, fca_i)
}
seuratObj@meta.data$cdr3s_nt <- NA
seuratObj@meta.data$cdr3s_nt <- tab$cdr3s_nt[match(rownames(seuratObj@meta.data), tab$barcode)]

###############################################################################
#'                        draw top 100 cclonotype bar                        '#
###############################################################################

cn <- as_tibble(seuratObj@meta.data) %>%
    filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "Blood") %>%
    group_by(cdr3s_nt) %>% dplyr::count() %>% arrange(desc(n))
cnTop100 <- cn[1:100,]
gdata <- as_tibble(seuratObj@meta.data) %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "Blood") %>%
    group_by(newClusterID, orig.ident, cdr3s_nt) %>% dplyr::count() %>%
    filter(cdr3s_nt %in% cnTop100$cdr3s_nt) %>% group_by(newClusterID, orig.ident) %>%
    dplyr::summarise(n = sum(n))
gdata$orig.ident <- as.factor(gdata$orig.ident)
gdata$orig.ident <- factor(gdata$orig.ident, levels=c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"))
gdata$orig.ident <- plyr::mapvalues(gdata$orig.ident,
          from = c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"),
          to =   c("Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23",
                   "Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23"))
gdata <- bind_rows(gdata, tibble(newClusterID = 4, orig.ident = "Pt_21", n = 0))
gdata <- bind_rows(gdata, tibble(newClusterID = 5, orig.ident = "Pt_21", n = 0))
gdata$newClusterID <- factor(gdata$newClusterID, levels = 1:17)
g <- ggplot(data = gdata) +
    geom_bar(stat = "identity", mapping = aes(x = newClusterID, y = n, fill = orig.ident)) +
    theme_classic()
pdf(file.path(figurePath, "Blood_patient_wise_barplot_100TopTCRClonotype.pdf"))
print(g)
dev.off()

cn <- as_tibble(seuratObj@meta.data) %>%
    filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "SF") %>%
    group_by(cdr3s_nt) %>% dplyr::count() %>% arrange(desc(n))
cnTop100 <- cn[1:100,]
gdata <- as_tibble(seuratObj@meta.data) %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "SF") %>%
    group_by(newClusterID, orig.ident, cdr3s_nt) %>% dplyr::count() %>%
    filter(cdr3s_nt %in% cnTop100$cdr3s_nt) %>% group_by(newClusterID, orig.ident) %>%
    dplyr::summarise(n = sum(n))
gdata$orig.ident <- as.factor(gdata$orig.ident)
gdata$orig.ident <- factor(gdata$orig.ident, levels=c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"))
gdata$orig.ident <- plyr::mapvalues(gdata$orig.ident,
          from = c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"),
          to =   c("Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23",
                   "Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23"))
gdata <- bind_rows(gdata, tibble(newClusterID = 1, orig.ident = "Pt_21", n = 0))
gdata <- bind_rows(gdata, tibble(newClusterID = 2, orig.ident = "Pt_21", n = 0))
gdata <- bind_rows(gdata, tibble(newClusterID = 6, orig.ident = "Pt_21", n = 0))
gdata <- bind_rows(gdata, tibble(newClusterID = 16, orig.ident = "Pt_21", n = 0))
gdata$newClusterID <- factor(gdata$newClusterID, levels = 1:17)
g <- ggplot(data = gdata) + geom_bar(stat = "identity", mapping = aes(x = newClusterID, y = n, fill = orig.ident)) + theme_classic()
pdf(file.path(figurePath, "SF_patient_wise_barplot_100TopTCRClonotype.pdf"))
print(g)
dev.off()


cn <- as_tibble(seuratObj@meta.data) %>%
    filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "Blood") %>%
    group_by(cdr3s_nt) %>% dplyr::count() %>% arrange(desc(n))
cnTop100 <- cn[1:100,]
seuratObj@meta.data$cnTop100Clones <- FALSE
seuratObj@meta.data$cnTop100Clones[seuratObj$cdr3s_nt %in% cnTop100$cdr3s_nt] <- TRUE
coord = Embeddings(object = seuratObj, reduction = "umap")
colnames(coord) = c('UMAP1','UMAP2')
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
colors <- scales::hue_pal()(length(unique(meta$newClusterID)))
meta$color <- ""
meta$color <- colors[match(meta$newClusterID, unique(seuratObj$newClusterID))]
metaWithTCROverlap = meta %>% filter(cnTop100Clones == TRUE)
metaWithoutTCROverlap = meta %>% filter(cnTop100Clones == FALSE)
metaWithTCROverlap$newClusterID <- as.factor(metaWithTCROverlap$newClusterID)
metaWithoutTCROverlap$newClusterID <- as.factor(metaWithoutTCROverlap$newClusterID)
label <- meta %>%
    group_by(newClusterID) %>%
    select(UMAP1, UMAP2) %>%
    summarize_all(mean)
g <- ggplot() +
    geom_point(data = metaWithoutTCROverlap %>% filter(control =="Blood") ,
               mapping = aes(x = UMAP1, y = UMAP2),
               color = "#EEEEEE",
               size = 0.1) +
    geom_point(data = metaWithTCROverlap%>% filter(control =="Blood"),
               mapping = aes(x = UMAP1, y = UMAP2, color = color), size = 0.2) +
    scale_color_identity() +
    ggrepel::geom_text_repel(data = label, aes(x = UMAP1, y = UMAP2, label = newClusterID)) +
    theme_classic()
png(file.path(figurePath, "T_DimPlot_mapped_TCR_PB.png"), width = 210/2.2, height = 297/3.2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
print(g)
dev.off()


cn <- as_tibble(seuratObj@meta.data) %>%
    filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "SF") %>%
    group_by(cdr3s_nt) %>% dplyr::count() %>% arrange(desc(n))
cnTop100 <- cn[1:100,]
seuratObj@meta.data$cnTop100Clones <- FALSE
seuratObj@meta.data$cnTop100Clones[seuratObj$cdr3s_nt %in% cnTop100$cdr3s_nt] <- TRUE
coord = Embeddings(object = seuratObj, reduction = "umap")
colnames(coord) = c('UMAP1','UMAP2')
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
colors <- scales::hue_pal()(length(unique(meta$newClusterID)))
meta$color <- ""
meta$color <- colors[match(meta$newClusterID, unique(seuratObj$newClusterID))]
metaWithTCROverlap = meta %>% filter(cnTop100Clones == TRUE)
metaWithoutTCROverlap = meta %>% filter(cnTop100Clones == FALSE)
metaWithTCROverlap$newClusterID <- as.factor(metaWithTCROverlap$newClusterID)
metaWithoutTCROverlap$newClusterID <- as.factor(metaWithoutTCROverlap$newClusterID)
g <- ggplot() +
    geom_point(data = metaWithoutTCROverlap %>% filter(control =="SF") ,
               mapping = aes(x = UMAP1, y = UMAP2),
               color = "#EEEEEE",
               size = 0.1) +
    geom_point(data = metaWithTCROverlap%>% filter(control =="SF"),
               mapping = aes(x = UMAP1, y = UMAP2, color = color), size = 0.2) +
    scale_color_identity() +
    ggrepel::geom_text_repel(data = label, aes(x = UMAP1, y = UMAP2, label = newClusterID)) +
    theme_classic()
png(file.path(figurePath, "T_DimPlot_mapped_TCR_SF.png"), width = 210/2.2, height = 297/3.2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
print(g)
dev.off()

meta$newClusterID <- factor(meta$newClusterID, levels = unique(seuratObj$newClusterID))
g <- ggplot() +
    geom_point(data = meta,
               mapping = aes(x = UMAP1, y = UMAP2, color = newClusterID), size = 0.2) +
    theme_classic()
png(file.path(figurePath, "meta.png"), width = 210/2.2, height = 297/3.2,
    units = "mm", pointsize = 8, bg = "white", res = 600)
print(g)
dev.off()
