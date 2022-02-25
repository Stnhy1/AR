#' filename : annotateMarkers.r
#' Date : 2020-04-13
#' contributor : Yanshuo Chu

##libraries
suppressMessages({
    library(optparse)
    library(rjson)
    library(Seurat)

    library(org.Hs.eg.db)
    library(mygene)
})

print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'visualization.pdf',
                help = 'output file name for the plot [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character'),
    make_option(c("-m","--marker"),
                type = 'character',
                help = 'marker celltype data frame',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}
if(is.null(opt$param)) {
    print_help(opt_parser)
    stop("json file name (containing parameters) must be provided", call. = F)
}


##load param
param <- fromJSON(file = opt$param)

## markersTable <- read.table(opt$marker, sep = "\t", head = T)


markersTable <- read.table("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/original/markers.top.tsv", sep = "\t", head = T)

markersTable$gene <- toupper(markersTable$gene)

mid <- mapIds(org.Hs.eg.db, as.vector(markersTable$gene), 'ENTREZID', 'SYMBOL')

clusters <- unique(markersTable$cluster)

pdf(opt$out, width = param$width, height = param$height)

for(clt in clusters){

    print("cluster: ")
    print(clt)

    tempMarkerList <- as.vector(markersTable$marker[markersTable$cluster == clt])

    print("marker list: ")
    print(tempMarkerList)

    ## print(FeaturePlot(seuratobj, reduction = param$reduction, features = tempMarkerList))
}
dev.off()
