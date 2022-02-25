#' filename : manuallySplitCD4AND8.R
#' Date : 2020-06-21
#' contributor : Yanshuo Chu
#' function: manually split CD4 and CD8

suppressMessages({library(optparse)
    library(rjson)
    library(rlist)
    library(dplyr)
    library(readr)
    library(Seurat)
    library(ggplot2)
    library(RColorBrewer)
    library(ggpubr)
})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratObj <- readRDS(opt$data)

## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_TCR_doublets/wholedata/data.rds")

exp.cut <- 1

cellNames <- rownames(GetAssayData(seuratObj))
CD3Markers <- grep("^CD3([^0-9]|$)", cellNames, perl = T, value = T)

assayData <- GetAssayData(seuratObj)
data <- as.matrix(assayData[rownames(assayData) %in% CD3Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
TCells <- colnames(seuratObj)[colSums(data>exp.cut)>0]


assayData <- GetAssayData(seuratObj)
CD4Markers <- grep("^CD4([^0-9]|$)", cellNames, perl = T, value = T)
data<-as.matrix(assayData[rownames(assayData) %in% CD4Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
CD4Cells <-colnames(seuratObj)[colSums(data>exp.cut)>0]


assayData <- GetAssayData(seuratObj)
CD8Markers <- grep("^CD8([^0-9]|$)", cellNames, perl = T, value = T)
data<-as.matrix(assayData[rownames(assayData) %in% CD8Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
CD8Cells <-colnames(seuratObj)[colSums(data>exp.cut)>0]


CD4TCells <- intersect(TCells, CD4Cells)
CD8TCells <- intersect(TCells, CD8Cells)

CD4TandCD8T <- intersect(CD4TCells, CD8TCells)

CD4TCells <- setdiff(CD4TCells, CD4TandCD8T)
CD4TObj <- subset(seuratObj, cells = CD4TCells)
saveRDS(CD4TObj, paste0(dirname(opt$data), '/', 'CD4TObj.rds'))


CD8TCells <- setdiff(CD8TCells, CD4TandCD8T)
CD8TObj <- subset(seuratObj, cells = CD8TCells)
saveRDS(CD8TObj, paste0(dirname(opt$data), '/', 'CD8TObj.rds'))


assayData <- GetAssayData(seuratObj)
CD16Markers <- grep("^FCGR3([^0-9]|$)", cellNames, perl = T, value = T)
data<-as.matrix(assayData[rownames(assayData) %in% CD16Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
CD16Cells <-colnames(seuratObj)[colSums(data>exp.cut)>0]
CD16Obj <- subset(seuratObj, cells = CD16Cells)
saveRDS(CD16Obj, paste0(dirname(opt$data), '/', 'CD16Obj.rds'))


assayData <- GetAssayData(seuratObj)
CD19Markers <- grep("^CD19([^0-9]|$)", cellNames, perl = T, value = T)
data<-as.matrix(assayData[rownames(assayData) %in% CD19Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
CD19Cells <-colnames(seuratObj)[colSums(data>exp.cut)>0]
CD19Obj <- subset(seuratObj, cells = CD19Cells)
saveRDS(CD19Obj, paste0(dirname(opt$data), '/', 'CD19Obj.rds'))



assayData <- GetAssayData(seuratObj)
CD56Markers <- grep("^NCAM1([^0-9]|$)", cellNames, perl = T, value = T)
data<-as.matrix(assayData[rownames(assayData) %in% CD56Markers,])
if(dim(data)[2] == 1){
    data <- t(data)
}
CD56Cells <-colnames(seuratObj)[colSums(data>exp.cut)>0]
CD56Obj <- subset(seuratObj, cells = CD56Cells)
saveRDS(CD56Obj, paste0(dirname(opt$data), '/', 'CD56Obj.rds'))

