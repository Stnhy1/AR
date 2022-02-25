#' filename : addDoubletInfo.R
#' Date : 2020-05-11
#' contributor : Yanshuo Chu
#' function: add doublet infor

##libraries
suppressMessages({library(optparse)
    library(readr)
    library(stringr)
    library(rjson)
    library(Seurat)
    library(tidyverse)
})
print('---snn clustering---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("--tenXfolder"),
                type = "character",
                default = NULL,
                help = "tenXfolder to get the doublet infor for each sample",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'seuratobjwithbatch.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

## md$sct_doublet_scores=-1
## md$sct_predicted_doublets="False"

seuobj <- readRDS(opt$data)
## seuobj <- readRDS('../../outs/5GE_VI/original/umap.rds')

md <- seuobj@meta.data
md$sbt_is_doublet="False"
md$sbt_doublet_score=-1

###############################################################################
                                        #           load data info            #
###############################################################################
## DataD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_CC'

tenXfolder=opt$tenXfolder
## tenXfolder=DataD

print(paste0("Loading data from ", tenXfolder))

samples = list.dirs(path = tenXfolder, recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ',tenXfolder))

for (ids in seq_along(samples)) {

    ## ids=1

    ds = samples[ids]
    print(ds)

    sbtPath <- file.path(tenXfolder, ds, "outs", "filtered_feature_bc_matrix", "scrubletTable.txt")
    sbt <- read.table(sbtPath, header = F, sep = "\t")
    colnames(sbt) <- c("barcode", "sbt_is_doublet", "sbt_doublet_score")
    rownames(sbt) <- paste0(ds,"_", str_remove_all(sbt$barcode, "[-1]"))

    md$sbt_is_doublet[rownames(md) %in% rownames(sbt)] = sbt$sbt_is_doublet=="True"
    md$sbt_doublet_score[rownames(md) %in% rownames(sbt)] = sbt$sbt_doublet_score
}

md$sbt_is_doublet[md$orig.ident=="62" & md$sbt_doublet_score > 0.22] = T
md$sbt_is_doublet[md$orig.ident=="62B" & md$sbt_doublet_score > 0.4] = T
md$sbt_is_doublet[md$orig.ident=="76" & md$sbt_doublet_score > 0.3] = T

seuobj@meta.data <- md

mdt<- as_tibble(md) %>% group_by(orig.ident) %>% count(sbt_is_doublet)
g <- ggplot(data=mdt, aes(x=orig.ident, y=n, fill=sbt_is_doublet)) + geom_bar(stat="identity")
pdf(paste0(opt$out, "_sbt_mdt.pdf"))
print(g)
dev.off()

saveRDS(seuobj, paste0(opt$out))
