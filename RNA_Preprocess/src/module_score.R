#'--------------------------------------------------------------
#' filename : module_score.R
#' Date : 2020-08-25
#' contributor : Yanshuo Chu
#' function: module_score
#'--------------------------------------------------------------

print('<==== module_score ====>')



suppressMessages({
    library(optparse)
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
})


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


moduleScorePlotFolder <- file.path(dirname(opt$data), "moduleScorePlot")
if(!dir.exists(moduleScorePlotFolder)){
    dir.create(moduleScorePlotFolder)
}

seuratObj <- readRDS(opt$data)
Idents(seuratObj) <- seuratObj@meta.data$seurat_clusters


memoryMarkers <- c("SELL", "CCR7", "CD28")
exhaustionMarkers <- c('CD274', 'PDCD1', 'LAG3', 'HAVCR2', 'CTLA4', 'TIGIT')
cytotoxicMarkers <- c('GZMB','GZMA', 'GNLY', 'PRF1')
naiveMarkers <- c('LEF1', "CCR7", "IL7R")

seuratObj <-  AddModuleScore(seuratObj,
                             features = list(
                                 memory = memoryMarkers,
                                 exhaustion = exhaustionMarkers,
                                 cytotoxicMarkers = cytotoxicMarkers,
                                 naiveMarkers = naiveMarkers
                             ),
                             ctrl = 5,
                             name = "ModuleScore")



coord = Embeddings(object = seuratObj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP1", "UMAP2")
coord = data.frame(ID = rownames(coord), coord)
meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')

colnames(meta)[which(names(meta) == "ModuleScore1")] <- "memoryScore"
colnames(meta)[which(names(meta) == "ModuleScore2")] <- "exhaustionScore"
colnames(meta)[which(names(meta) == "ModuleScore3")] <- "cytotoxicScore"
colnames(meta)[which(names(meta) == "ModuleScore4")] <- "naiveScore"

for(tempScore in c("memoryScore", "exhaustionScore", "cytotoxicScore", "naiveScore")){
    g <- ggplot() +
        geom_point(data = meta ,
                   mapping = aes_string(x = "UMAP1", y = "UMAP2", color = tempScore),
                   size = 0.1) + scale_colour_gradientn(colours = c("blue", "yellow", "red"))
        ggtitle(paste0(tempScore)) + theme_classic()

    pdf(file.path(moduleScorePlotFolder, paste0(tempScore, ".pdf")))
    print(g)
    dev.off()
}
