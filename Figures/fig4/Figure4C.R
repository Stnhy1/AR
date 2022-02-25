###############################################################################
#'                                 circus plot                               '#
###############################################################################

rm(list=ls())
library(Seurat)
library(tidyverse)
library(iTALK)
library(plyr)
library(scales)
library(RColorBrewer)
library(MASS)

pal = colorRampPalette(c("blue", "red"))
figurePath <- file.path("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/manuscript/fig4/outs")
if(!dir.exists(figurePath)){
    dir.create(figurePath, recursive = T)
}
setwd(figurePath)

dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/dataForpaper.rds"
seuratObj<- readRDS(dataPath)

Idents(seuratObj) <- seuratObj@meta.data$orig.ident
seuratObj <- subset(seuratObj, idents = c("A13", "A50", "A13B", "A50B"), invert = T)
seuratObj@meta.data$newClusterID <- mapvalues(seuratObj@meta.data$newClusterID,
                                              from=c(1:16,18:31),
                                              to=c(1:30))
Idents(seuratObj) <- seuratObj@meta.data$newClusterID
seuratObj <- subset(seuratObj, idents = c(30), invert = T)

#' newClusterID  ##############################################################
cellTypeL <- list(
    T = c("Activated CD4", "Activated CD8", "CD38- CD8 T", "CD38+ CD8 T",
          "CX3CR1hi Effector CD8", "CX3CR1lo Effector CD8", "Cycling T",
          "Effector CD8", "Naive T", "Treg"),
    GDT = c("Gamma-Delta T"),
    NK = c("CD16- NK and NKT",
           "CD16+ NK and NKT"),
    B = c("Memory B",
          "Naive B",
          "Plasmablasts"),
    Myeloid = c(
        "Classical monocytes",
        "mDC",
        "Megakaryocytes",
        "Neutrophil",
        "Non-classical monocytes",
        "pDC"))

seuratObj$new.cell.types <- ""

for(i in 1:length(cellTypeL)){
    newCellTypeNames <- names(cellTypeL[i])
    oldCellTypeNames <- cellTypeL[[i]]

    seuratObj$new.cell.types[seuratObj$cell.types %in% oldCellTypeNames] <- newCellTypeNames
}

Idents(seuratObj) <- seuratObj@meta.data$control
BloodSeuratObj <- subset(seuratObj, ident = "Blood")
SFSeuratObj <- subset(seuratObj, ident = "SF")
topGeneNum <- 100
topLRPNum <- 100
targetList <- list(Whole = seuratObj, Blood = BloodSeuratObj, SF = SFSeuratObj)
markers <- list(one = c( "CXCR3", "CXCL9", "CXCL10", "CXCL11" ),
                two = c( "CXCR6", "CXCL16" ),
                three = c( "CCR1", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8" ))



getPValue <- function(res_cat, groupName){
    res_cat$p <- res_cat$cell_from_mean_exprs * res_cat$cell_to_mean_exprs
    res_cat$key <- paste0(res_cat$ligand, res_cat$receptor, res_cat$cell_from, res_cat$cell_to)

    disMatrix <- matrix(0, dim(res_cat)[1], ncol = 1000)

    for(i in 1:1000){
        wdatapath <- paste0("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/randomShuffleCells/result/", groupName, "_highly_exprs_genes_", i, ".tsv")
        tpt <- read_tsv(wdatapath)
        tpt$p <- tpt$cell_from_mean_exprs * tpt$cell_to_mean_exprs
        tpt$key <- paste0(tpt$ligand, tpt$receptor, tpt$cell_from, tpt$cell_to)
        tpt <- tpt[match(res_cat$key, tpt$key),]
        disMatrix[,i] <- tpt$p
    }

    res_cat$p_value = 0

    for(i in 1:dim(res_cat)[1]){
        fit <- fitdistr(disMatrix[i,], "normal")
        res_cat$p_value[i] <- pnorm(res_cat$p[i],
                                    mean = fit$estimate[1],
                                    sd = fit$estimate[2],
                                    lower.tail = FALSE)
    }
    return(res_cat)
}

pal = colorRampPalette(c("red", "blue"))
crs = pal(1001)
wrs = seq(1, 6, by = 0.05)

for(idi in 1: length(targetList)){
    groupNames <- names(targetList[idi])
    tempSeuratObj <-targetList[[idi]]

    hvg <- VariableFeatures(tempSeuratObj)
    data <- as.data.frame(t(as.matrix(tempSeuratObj@assays$RNA@data[hvg,])))

    data$cell_type <- ""
    data$cell_type <- tempSeuratObj@meta.data$new.cell.types[match(rownames(data), rownames(tempSeuratObj@meta.data))]

    highly_exprs_genes <- rawParse(data, top_genes=2000, stats='mean')

    for(idmi in 1: length(markers)){
        tempHEG <- highly_exprs_genes %>% filter(gene %in% markers[[idmi]])
        ## tempHEG <- highly_exprs_genes
        print(table(tempHEG$cell_type))

        pdf(file.path(figurePath, paste0("signicant_circusplot_wholeData_", groupNames, "_", idmi, ".pdf")), width = 12, height = 10)

        comm_list<-c('cytokine')
        cell_col<-structure(hue_pal()(length(unique(tempSeuratObj@meta.data$new.cell.types))), names= c("T", "Myeloid", "GDT", "NK", "B"))

        for(comm_type in comm_list){
            res_cat<-FindLR(tempHEG, datatype='mean count', comm_type=comm_type)
            res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
            res_cat <- getPValue(res_cat, groupNames)
            res_cat <- res_cat[res_cat$p > 0,]
            ## res_cat <- res_cat[res_cat$p_value < 0.05,]

            colorindex <- round(res_cat$p_value * 1000) + 1
            maxW <- max(res_cat$p)
            minW <- min(res_cat$p)
            normalizedW <- (res_cat$p - minW)/(maxW-minW)

            windex <- round(normalizedW * 100) + 1

            LRPlot(res_cat,
                   datatype='mean count',
                   cell_col=cell_col,
                   link.arr.lwd=wrs[windex],
                   link.arr.width=wrs[windex],
                   ## link.arr.col=pal(dim(res_cat)[1])
                   link.arr.col=crs[colorindex])

            lgd_ = rep(NA, 200)
            lgd_[c(1, 200)] = c(1, 0)
            legend(x = 1.2,
                   y = 0.5,
                   xjust = 1,
                   yjust = 1,
                   legend = lgd_,
                   fill = colorRampPalette(colors = c('blue','red'))(200),
                   border = NA, bty = "n",
                   y.intersp = 0.01,
                   cex = 1, text.font = 1,
                   title = NULL,
                   title.adj = 0)

            legend("topleft",
                   legend = round(quantile(res_cat$p), digits = 5),
                   lwd = quantile(1:6), cex = 0.8)
        }
        dev.off()
    }
}



getPValue <- function(res_cat, groupName){
    res_cat$p <- res_cat$cell_from_mean_exprs * res_cat$cell_to_mean_exprs
    res_cat$key <- paste0(res_cat$ligand, res_cat$receptor, res_cat$cell_from, res_cat$cell_to)

    disMatrix <- matrix(0, dim(res_cat)[1], ncol = 1000)

    for(i in 1:1000){
        wdatapath <- paste0("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/plotscriptsforpaper/randomShuffleCells/resultPatient/", groupName, "_highly_exprs_genes_", i, ".tsv")
        tpt <- read_tsv(wdatapath)
        tpt$p <- tpt$cell_from_mean_exprs * tpt$cell_to_mean_exprs
        tpt$key <- paste0(tpt$ligand, tpt$receptor, tpt$cell_from, tpt$cell_to)
        tpt <- tpt[match(res_cat$key, tpt$key),]
        disMatrix[,i] <- tpt$p
    }

    res_cat$p_value = 0

    for(i in 1:dim(res_cat)[1]){
        fit <- fitdistr(disMatrix[i,], "normal")
        res_cat$p_value[i] <- pnorm(res_cat$p[i],
                                    mean = fit$estimate[1],
                                    sd = fit$estimate[2],
                                    lower.tail = FALSE)
    }
    return(res_cat)
}


totalResCat <- c()
for(oin in unique(seuratObj$orig.ident)){

    Idents(seuratObj) <- seuratObj$orig.ident
    tempSeuratObj <- subset(seuratObj, idents = oin)
    Idents(tempSeuratObj) <- tempSeuratObj$new.cell.types

    hvg <- VariableFeatures(tempSeuratObj)
    data <- as.data.frame(t(as.matrix(tempSeuratObj@assays$RNA@data[hvg,])))

    data$cell_type <- ""
    data$cell_type <- tempSeuratObj@meta.data$new.cell.types[match(rownames(data), rownames(tempSeuratObj@meta.data))]

    highly_exprs_genes <- rawParse(data, top_genes=2000, stats='mean')

    for(idmi in 1: length(markers)){
        tempHEG <- highly_exprs_genes %>% filter(gene %in% markers[[idmi]])
        ## tempHEG <- highly_exprs_genes
        print(table(tempHEG$cell_type))

        pdf(file.path(figurePath, paste0("signicant_circusplot_patient_", oin, "_", idmi, ".pdf")), width = 12, height = 10)

        comm_list <- c('cytokine')
        cell_col <- structure(hue_pal()(5), names= c("T", "Myeloid", "GDT", "NK", "B"))
        cell_col <- cell_col[unique(tempSeuratObj$new.cell.types)]

        for(comm_type in comm_list){
            res_cat <- FindLR(tempHEG, datatype='mean count', comm_type=comm_type)
            if(dim(res_cat)[1] == 0) {
                next
            }

            res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]

            res_cat <- getPValue(res_cat, oin)

            res_cat <- res_cat[res_cat$p > 0,]

            if(dim(res_cat)[1] == 0) {
                next
            }

            ## res_cat <- res_cat[res_cat$p_value < 0.05,]

            colorindex <- round(res_cat$p_value * 1000) + 1
            maxW <- max(res_cat$p)
            minW <- min(res_cat$p)
            normalizedW <- (res_cat$p - minW)/(maxW-minW)

            windex <- round(normalizedW * 100) + 1

            LRPlot(res_cat,
                   datatype='mean count',
                   cell_col=cell_col,
                   link.arr.lwd=wrs[windex],
                   link.arr.width=wrs[windex],
                   ## link.arr.col=pal(dim(res_cat)[1])
                   link.arr.col=crs[colorindex])

            lgd_ = rep(NA, 200)
            lgd_[c(1, 200)] = c(1, 0)
            legend(x = 1.2,
                   y = 0.5,
                   xjust = 1,
                   yjust = 1,
                   legend = lgd_,
                   fill = colorRampPalette(colors = c('blue','red'))(200),
                   border = NA, bty = "n",
                   y.intersp = 0.01,
                   cex = 1, text.font = 1,
                   title = NULL,
                   title.adj = 0)

            legend("topleft", legend = round(quantile(res_cat$p), digits = 5), lwd = quantile(1:6), cex = 0.8)

            res_cat$patient <- oin
            totalResCat <- rbind(totalResCat, res_cat)
        }
        dev.off()
    }
}

write_tsv(totalResCat, file.path(figurePath, paste0('totalResCat', "_", Sys.Date(), '.tsv')))

