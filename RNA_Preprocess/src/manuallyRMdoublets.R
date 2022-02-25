#' filename : manuallyRMdoublets.R
#' Date : 2020-06-03
#' contributor : Yanshuo Chu
#' function: remove doublets by markers

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
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'visualization.pdf',
                help = 'output file name for the plot [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
seuratobj <- readRDS(opt$data)

## seuratobj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/combined.data.rds")

B.markers <- c("CD79A", "CD79B");
T.markers <- c("CD3D", "CD3E");

gene<-c('EPCAM','PTPRC','COL1A1','CD79A','CD3D','S100A8','S100A9','CD8A');
com<-combn(gene,2);
com<-com[,-c(2,9:13,25,26)]
exp.cut=1
CELL.del<-NULL
for(i in 1:dim(com)[2]){
    print(i)
    data<-as.matrix(GetAssayData(seuratobj)[com[,i],])
    index<-colSums(data>exp.cut)==2
    cell.del<-colnames(seuratobj)[index];
    CELL.del<-c(CELL.del,cell.del)
}
CELL.del<-unique(CELL.del);

seuratobj@meta.data$manually_is_doublet <- rownames(seuratobj@meta.data) %in% CELL.del

print("finish manually find doublet")
print(colnames(seuratobj@meta.data))

md <- seuratobj@meta.data
mdt<- as_tibble(md) %>% group_by(orig.ident) %>% count(manually_is_doublet)
g <- ggplot(data=mdt, aes(x=orig.ident, y=n, fill=manually_is_doublet)) + geom_bar(stat="identity")
pdf(paste0(opt$out, "_manually_mdt.pdf"))
print(g)
dev.off()

md <- seuratobj@meta.data
seuratobj <- subset(seuratobj, cells = rownames(md)[!as.logical(md$sbt_is_doublet) & !as.logical(md$manually_is_doublet)] )

## write_tsv(CELL.del, paste0(opt$out, 'CELL.del.tsv'))
saveRDS(seuratobj,  opt$out)
