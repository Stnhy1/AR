#' filename : mergeTCR.R
#' Date : 2020-06-17
#' contributor : Yanshuo Chu
#' function: merge TCR
suppressMessages({library(optparse)
library(ggplot2)
library(readr)
library(tidyverse)
library(rjson)
library(Seurat)})
print('----loading cellranger outputs----')
##CLI parsing
option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'dataFolder',
                metavar = 'character'),
    make_option(c("-s","--seuratobj"),
                type = 'character',
                help = 'seurat object',
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character'),
    make_option(c("-f", '--sfout'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character'),
    make_option(c("-b",'--bloodout'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character'),
    make_option(c("-m",'--myeloidout'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

## Get all input samples

seuratObj <- readRDS(opt$seuratobj)


seuratObj@meta.data$hasTCR = FALSE

print(paste0("Loading data from ", opt$data))

samples = list.dirs(path = opt$data,recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ', opt$data))

obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)
    targetFile <- paste0(file.path(opt$data, ds, "outs"), "/filtered_contig_annotations.csv")
    barcodeS <- read_csv(targetFile) %>% select(barcode)
    barcodeS$barcode <- paste0(ds,"_", str_remove_all(barcodeS$barcode, "[-1]"))
    seuratObj@meta.data$hasTCR[rownames(seuratObj@meta.data) %in% barcodeS$barcode] = TRUE
}
print(table(seuratObj@meta.data$hasTCR))

pdf(paste0(opt$out, "_hasTCR.pdf"))
DimPlot(seuratObj, group.by='hasTCR', label=TRUE)
dev.off()


nonTcellClusters <- c(28, 29, 19, 6, 15, 30, 8, 23, 31, 18, 7, 25, 24, 22, 27, 26, 16, 5, 20)

md <- seuratObj@meta.data
mdt<- as_tibble(md) %>% filter(seurat_clusters %in% nonTcellClusters) %>% group_by(orig.ident) %>% count(hasTCR)
g <- ggplot(data=mdt, aes(x=orig.ident, y=n, fill=hasTCR)) + geom_bar(stat="identity")
pdf(paste0(opt$out, "_TCR.pdf"))
print(g)
dev.off()

seuratObj <- subset(seuratObj, cells = rownames(md)[!(as.logical(md$seurat_clusters %in% nonTcellClusters) & as.logical(md$hasTCR))])
saveRDS(seuratObj, opt$out)


SFObj <- subset(seuratObj, cells = rownames(seuratObj@meta.data)[as.logical(seuratObj@meta.data$control == "SF")])
saveRDS(SFObj, opt$sfout)

BloodObj <- subset(seuratObj, cells = rownames(seuratObj@meta.data)[as.logical(seuratObj@meta.data$control == "Blood")])
saveRDS(BloodObj, opt$bloodout)

myeloidClusters <- c(15, 19, 6, 30, 8, 23, 31, 18, 7, 25)
myeloidObj <- subset(seuratObj, cells = rownames(seuratObj@meta.data)[as.logical(seuratObj@meta.data$seurat_clusters %in% myeloidClusters)])
saveRDS(myeloidObj, opt$myeloidout)
