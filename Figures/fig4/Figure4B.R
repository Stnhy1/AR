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

cell.type.order <- c("CD69lo Naive CD4", "CD69hi Naive CD4", "Th1 and Th17-like", "Treg", "CXCL13+ T", "Naive CD8", "Central Memory CD8", "CXCR6+ Effector CD8", "CXCR6- Effector CD8", "CX3CR1+ Effector CD8", "CX3CR1- Effector CD8", "Cycling T", "MAIT T", "CX3CR1- Gamma-Delta T", "CX3CR1+ Gamma-Delta T", "CD16+ NK and NKT", "CD16- NK and NKT")

blood.id <- grep("B$", file.id, value = T)
SF.id <- setdiff(file.id, blood.id)
mdt <- as_tibble(seuratObj@meta.data)
stat <- mdt %>% group_by(orig.ident, control, GROUP) %>% count
Mono.id <- stat$orig.ident[stat$GROUP == "Mono"]
Combo.id <- stat$orig.ident[stat$GROUP == "Combo"]
OA.id <- stat$orig.ident[stat$GROUP == "OA"]
SF.id <- setdiff(SF.id, OA.id)
blood.id <- setdiff(blood.id, OA.id)
file.id <- setdiff(file.id, OA.id)

allIDsList <- list(all = file.id, blood = blood.id, SF = SF.id,
                   Mono = Mono.id, Combo = Combo.id, OA = OA.id)

for(idi in 1:length(allIDsList)){
    groupNames <- names(allIDsList[idi])
    groupSamples <- allIDsList[[idi]]

    md <- as.tibble(seuratObj@meta.data, rownames = NA) %>% mutate(rownames = rownames(seuratObj@meta.data))
    md <- md %>% arrange(factor(cell.types, levels = cell.type.order))

    allClonotypes <- unique(md$cdr3s_nt)

    ctn <- md %>% filter(!is.na(cdr3s_nt)) %>% group_by(cell.types, cdr3s_nt) %>% count
    ctn <- ctn %>% arrange(factor(cell.types, levels = cell.type.order))

    ctn.shared <- matrix(
        rep(NA, length(unique(ctn$cell.types)) ^ 2),
        length(unique(ctn$cell.types)),
        length(unique(ctn$cell.types)))

    p.ctn.shared <- matrix(
        rep(-1, length(unique(ctn$cell.types)) ^ 2),
        length(unique(ctn$cell.types)),
        length(unique(ctn$cell.types)))

    colnames(ctn.shared) <- unique(ctn$cell.types)
    rownames(ctn.shared) <- unique(ctn$cell.types)

    p.ctn.shared.tibble <- tibble(from = character(), to = character(), pvalue = double())
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            cti <- colnames(ctn.shared)[i]
            ctj <- colnames(ctn.shared)[j]
            cti_cts <- ctn %>% filter(cell.types == cti)
            ctj_cts <- ctn %>% filter(cell.types == ctj)

            ctij_cts <- intersect(cti_cts$cdr3s_nt, ctj_cts$cdr3s_nt)
###############################################################################
            #'                         draw pie chart                        '#
            tempPieData <- md %>% filter(cdr3s_nt %in% ctij_cts) %>% group_by(orig.ident) %>% count()
            if(dim(tempPieData)[1] > 0){

                g <- ggplot(tempPieData, aes(x = factor(""), y = n, fill = orig.ident)) +
                    geom_bar(position="fill", stat="identity") +
                    coord_polar("y") + ggtitle(paste0(cti, "_", ctj))

                tempPieDir <- file.path(figurePath, paste0(groupNames, "_matrix_pie"))
                if(!dir.exists(tempPieDir)){
                    dir.create(tempPieDir)
                }
                pdf(file.path(tempPieDir, paste0("matrix_", i, "_", j, "_cluster_pie.pdf")), width = 20, height = 20)
                print(g)
                dev.off()
            }
###############################################################################
            ctn.shared[i, j] <- length(ctij_cts)
            ctn.shared[j, i] <- ctn.shared[i, j]

            ijTest <- matrix(c(
                    length(ctij_cts),
                    length(ctj_cts$cdr3s_nt) - length(ctij_cts),
                    length(cti_cts$cdr3s_nt) - length(ctij_cts),
                    length(setdiff(allClonotypes, union(cti_cts$cdr3s_nt, ctj_cts$cdr3s_nt)))),
                    nrow = 2,
                    dimnames = list(iCluster = c("In", "NotIn"),
                                    jCluster = c("In", "NotIn")))

            ft <- fisher.test(ijTest, alternative = "greater")
            p.ctn.shared[i, j] <- ft[[1]]
            p.ctn.shared[j, i] <- ft[[1]]
            p.ctn.shared.tibble <- p.ctn.shared.tibble %>%
                add_row(from = colnames(ctn.shared)[i], to = colnames(ctn.shared)[j], pvalue = ft[[1]])
        }
    }

    p_values <- c()
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            p_values <- c(p_values, p.ctn.shared[i,j])
        }
    }
    p_values_adjust <- p.adjust(p_values, method = "BH")

    k = 1
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            p.ctn.shared[i,j] <- p_values_adjust[k]
            p.ctn.shared[j,i] <- p_values_adjust[k]

            p.ctn.shared.tibble$pvalue[p.ctn.shared.tibble$from == colnames(ctn.shared)[i] & p.ctn.shared.tibble$to == colnames(ctn.shared)[j]] == p_values_adjust[k]

            k <- k + 1
        }
    }

    p.ctn.shared.tibble <- p.ctn.shared.tibble %>% mutate(isAdjPSig = pvalue < 0.05, nlog10adjp=-log10(pvalue))

    p1<-vis.heatmap(ctn.shared, .title = "shared clonotypes",
                    .labs = c("", ""), .legend = "# clonotypes")
    tcrPlotFolder <- file.path(dirname(dataPath), "tcrPlot")
    if(!dir.exists(tcrPlotFolder)){
        dir.create(tcrPlotFolder)
    }

    write_tsv(as_tibble(p.ctn.shared), file.path(paste0("T_", groupNames, '_adj_p_value_sharedClonotypes.tsv')), col_names = F)
    write_tsv(as_tibble(p.ctn.shared.tibble), file.path(paste0("Edge_", groupNames, '_adj_p_value_sharedClonotypes.tsv')), col_names = T)

    pdf(file.path(paste0("T_", groupNames, "_sharedClonotypes_matrix.pdf")), width = 10, height = 10)
    print(p1)
    dev.off()
}
