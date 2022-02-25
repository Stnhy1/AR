#'--------------------------------------------------------------
#' filename : AlluvialPlot.R
#' Date : 2021-03-09
#' contributor : Yanshuo Chu
#' function: AlluvialPlot
#'--------------------------------------------------------------

print('<==== AlluvialPlot ====>')

##libraries
suppressMessages({library(optparse)
    library(tidyverse)
    library(Seurat)
    library(ggplot2)
    library(ggalluvial)
})
print('---snn clustering---')
##CLI parsing
option_list = list(
    make_option(c("-l", "--olddata"),
                type = "character",
                default = NULL,
                help = "old seuratObj",
                metavar = 'character'),
    make_option(c("-r", "--newdata"),
                type = "character",
                default = NULL,
                help = "new seuratObj",
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

##Load data

oldSeuratObj <- readRDS(opt$olddata)
newSeuratObj <- readRDS(opt$newdata)


metaOrigin <- oldSeuratObj@meta.data
metaValidate <- newSeuratObj@meta.data

metaValidate$seurat_clusters_old <- -1
metaValidate$seurat_clusters_old <- metaOrigin[rownames(metaValidate), "seurat_clusters"]

sn = 5000
for(ti in 1:2){
    if (ti == 1){
        if(min(table(metaValidate$seurat_clusters)) < sn){
            sn = min(table(metaValidate$seurat_clusters))
        }
        mv <- metaValidate %>% group_by(seurat_clusters) %>% dplyr::sample_n(sn)
        mv$Freq <- 1
    }else{
        if(min(table(metaValidate$seurat_clusters_old)) < sn){
            sn = min(table(metaValidate$seurat_clusters_old))
        }
        mv <- metaValidate %>% group_by(seurat_clusters_old) %>% dplyr::sample_n(sn)
        mv$Freq <- 1
    }
    g <- ggplot(data = mv,
                aes(axis1 = seurat_clusters_old, axis2 = seurat_clusters, y = Freq)) +
        scale_x_discrete(limits = c("seurat_clusters_old", "seurat_clusters"), expand = c(.2, .05)) +
        xlab("clusters") +
        ylab("cell#") +
        geom_alluvium(aes(fill = seurat_clusters)) +
        geom_stratum() +
        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_minimal() 
    if(ti == 1){
        pdf(file.path(dirname(opt$newdata), paste0("new_AlluviumPlot.pdf")))
        print(g)
        dev.off()
    }else{
        pdf(file.path(dirname(opt$newdata), paste0("old_AlluviumPlot.pdf")))
        print(g)
        dev.off()
    }
}
