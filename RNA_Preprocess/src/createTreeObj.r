#' filename : createTreeObj.r
#' Date : 2020-06-16
#' contributor : Yanshuo Chu
#' function: create tree object

suppressMessages({
    library(optparse)
    library(readr)
    library(tidyverse)
    library(rjson)
    library(Seurat)
    library(dplyr)
    library(alakazam)
    library(igraph)
    library(visNetwork)
})
option_list = list(
    make_option(
        c("-o", "--outputDir"), type="character", default=NULL,
        help="output directory", metavar="character"),
    make_option(
        c("-g", "--germPassTable"), type="character", default=NULL,
        help="germline pass table", metavar="character"),
    make_option(
        c("-c", "--cloneOutputFile"), type="character", default=NULL,
        help="clone pass table", metavar="character"),
    make_option(
        c("-d", "--dnapars_exec"), type="character", default=NULL,
        help="dnapars runnable file", metavar="character"),
    make_option(
        c("-p", "--subDFOutputFile"), type="character", default=NULL,
        help="Output: subclone data frame data file", metavar="character"),
    make_option(
        c("-t", "--treeOutputFile"), type="character", default=NULL,
        help="Output: Tree data file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(!dir.exists(opt$outputDir)) dir.create(opt$outputDir)

## setwd('/Users/ychu2/Projects/immuTree/TENX_076/outs/immcantation')
## germPassTable = read_tsv('/Users/ychu2/Projects/immuTree/TENX_076/outs/immcantation/s076_heavy_FUNCTIONAL-T.clonal_germ-pass.tab')

setwd(opt$outputDir)
germPassTable = read_tsv(opt$germPassTable)
colnames(germPassTable)[1] <- tolower(colnames(germPassTable)[1])

cloneSize <- germPassTable %>% count(CLONE) %>% arrange(desc(n))
print(cloneSize)
print(cloneSize[1,1])

germPassTable.clone1 = filter(germPassTable, CLONE == cloneSize$CLONE[1])

print(dim(germPassTable.clone1))

clone = makeChangeoClone(
    germPassTable.clone1,
    pad_end =F,
    text_fields = c('C_CALL','JUNCTION_10X_AA'),
    num_fields = c('UMICOUNT'))

## saveRDS(clone, '/Users/ychu2/Projects/immuTree/TENX_076/outs/immcantation/s013_clone1-data.rds')
## tree <- buildPhylipLineage(clone, '/Users/ychu2/miniconda2/bin/dnapars', rm_temp=TRUE)

saveRDS(clone, opt$cloneOutputFile)
tree <- buildPhylipLineage(clone, opt$dnapars_exec, rm_temp=TRUE)

V(tree)$color <- "steelblue"
V(tree)$color[V(tree)$name == "Germline"] <- "black"
V(tree)$color[grepl("Inferred", V(tree)$name)] <- "gray"

vdj <- germPassTable.clone1[germPassTable.clone1$SEQUENCE_ID %in% V(tree)$name, c("SEQUENCE_ID","SEQUENCE_VDJ")]
vdj$ID <- paste("n_",seq.int(nrow(vdj)),sep="")
subDF <- germPassTable.clone1[, c("SEQUENCE_ID","SEQUENCE_VDJ")]
subDF$ID <- vdj$ID[match(subDF$SEQUENCE_VDJ, vdj$SEQUENCE_VDJ)]
subDF$represent_ID <- vdj$SEQUENCE_ID[match(subDF$SEQUENCE_VDJ, vdj$SEQUENCE_VDJ)]
subDF <- subDF[,c("ID","represent_ID","SEQUENCE_ID")]
colnames(subDF) <- c("NodeID","RepresentCell","CellGroup")
subDF <- subDF[complete.cases(subDF),]

saveRDS(subDF, opt$subDFOutputFile)

size = clone@data$COLLAPSE_COUNT[match(V(tree)$name, clone@data$SEQUENCE_ID)] + 5
size[is.na(size)] = 2
size[size > 20]=25
V(tree)$size = size

vlabel <- names(V(tree))
for(i in 1:length(vlabel)){
    if(vlabel[i] %in% subDF$RepresentCell){
        vlabel[i] <- subDF$NodeID[subDF$RepresentCell == names(V(tree))[i]]
    }
}

pdf(file.path(dirname(opt$treeOutputFile), "tree.pdf"))
print(plot(tree))
dev.off()

treeCV <- tree %>% set_vertex_attr("label", value = vlabel)
treeVIS <- visIgraph(treeCV, layout="layout_as_tree", idToLabel = FALSE)

saveRDS(treeVIS, opt$treeOutputFile)
