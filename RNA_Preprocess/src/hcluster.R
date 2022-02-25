#'--------------------------------------------------------------
#' filename : hcluster.R
#' Date : 2020-11-05
#' contributor : Yanshuo Chu
#' function: hcluster
#'--------------------------------------------------------------

print('<==== hcluster ====>')

suppressMessages({
    library(optparse)
    library(tidyverse)
    library(Seurat)
    library(ggplot2)
    library(Matrix)
})

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

seuratObj <- readRDS(opt$data)

Idents(seuratObj) <- seuratObj$seurat_clusters

featureClusterMatrix <- NULL
for(i in levels(Idents(seuratObj))){
    cellsClusteri <- rownames(seuratObj@meta.data)[seuratObj@meta.data$seurat_clusters == i]
    featurei <- rowMeans(seuratObj@assays$RNA@data[,cellsClusteri])
    featureClusterMatrix <- rbind(featureClusterMatrix, featurei)
    print(dim(featureClusterMatrix))
}
rownames(featureClusterMatrix) <- paste("cluster", levels(Idents(seuratObj)))

hc <- hclust(dist(featureClusterMatrix))
pdf(file.path(dirname(opt$data), "cluster_hclust.pdf"))
plot(hc)
dev.off()

