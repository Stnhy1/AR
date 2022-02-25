#'--------------------------------------------------------------
#' filename : findmarker.r
#' Date : 2020-07-30
#' contributor : Yanshuo Chu
#' function: findmarker
#'--------------------------------------------------------------

print('<==== findmarker ====>')



##libraries
suppressMessages({
    library(optparse)
    library(Seurat)
    library(tidyverse)
})
print('---calculate snn cluster makers---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after snn clustering)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'markers.txt',
                help = 'output file name for the markers [default = %default]',
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
## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/Tcellproject/analysis/rPCAoutsModifyPC/nPC_50/UMAP_dist_0.1_nneighbor_30/p1RRPCA_UMAP_dist_0.1_nneighbor_30_CLUSTER_res_0.8/CD4_V3/nPC_50/UMAP_dist_0.1_nneighbor_50/p1CD4_V3_50_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.3/cluster.rds")

Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters
markerList <- list()
metaData <- as_tibble(seuratObj@meta.data)
all_clusters <- unique(seuratObj@meta.data$seurat_clusters)
for(i in all_clusters){
    if(dim(filter(metaData, seurat_clusters == i))[1] < 20){
        next
    }else{
        iMarkers <- FindMarkers(seuratObj, ident.1 = i, ident.2 = NULL)
        iMarkers$cluster = i
        iMarkers$gene = rownames(iMarkers)
        markerList[[length(markerList) + 1]] <- iMarkers
        print(paste0("Finish cluster ", i))
    }
}

markers <- bind_rows(markerList)

write_tsv(data.frame(markers), opt$out)
print('---end---')
