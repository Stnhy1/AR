#'--------------------------------------------------------------
#' filename : slingshot.R
#' Date : 2020-11-10
#' contributor : Yanshuo Chu
#' function: slingshot
#'--------------------------------------------------------------

print('<==== slingshot ====>')


suppressMessages({
    library(optparse)
    library(slingshot)
    library(Seurat)
    library(stringr)
    library(slingshot)
    library(mclust)
    library(tidyverse)
    library(tidymodels)
    library(viridis)
    library(ggpubr)
    library(scales)
})
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "seurat r data file after normalization",
                metavar = 'character'),
    make_option(c("-r", "--reduction"),
                type = "character",
                default = "umap",
                help = "reduction used by slingshot, from seurat obj",
                metavar = 'character'),
    make_option(c("-m", "--dim"),
                type = "integer",
                default = 2,
                help = "maximum dimension used",
                metavar = 'integer'),
    make_option(c("-g", "--groupby"),
                type = "character",
                default = "seurat_clusters",
                help = "cell clusters",
                metavar = 'character'),
    make_option(c("-s", "--startcluster"),
                type = "integer",
                default = 3,
                help = "start cell clusters",
                metavar = 'integer'),
    make_option(c("-o", "--outDir"),
                type = "character",
                help = "output directory",
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}


seuratObj<- readRDS(opt$data)

seuratObj = subset(seuratObj, cells = sample(Cells(seuratObj), 3000))

RunSlingshot <- function(seuratObj,
                         reduction = opt$reduction,
                         dims = 1:opt$dim,
                         group.by = opt$groupby,
                         output = "seurat",
                         start.cluster = opt$startcluster,
                         stretch = 0,
                         thresh = 0.001) {
    # Extract counts from seurat object. Only works for integrated data.
    count.data <-as.matrix(seuratObj@assays[[seuratObj@active.assay]]@data)

    # Create a single cell experiment object from count data
    sce.object <-
      SingleCellExperiment::SingleCellExperiment(assays = base::list(counts = count.data))

    # Add reduced dimensions to single cell experiment
    rd1 <- as.matrix(seuratObj@reductions[[reduction]]@cell.embeddings[,dims])

    SingleCellExperiment::reducedDims(sce.object) <- S4Vectors::SimpleList(reduced_data = rd1)

    # Cluster cells (Seurat requires clustering)
    if (group.by == "recluster") {
      SummarizedExperiment::colData(sce.object)[["label"]] <-
        mclust::Mclust(rd1)$classification
    } else {
      SummarizedExperiment::colData(sce.object)[["label"]] <-
        as.factor(as.matrix(seuratObj[[group.by]]))
    }
    # Run slingshot
    sce.slingshot <- slingshot::slingshot(sce.object,
                               clusterLabels = "label",
                               reducedDim = 'reduced_data',
                               start.clus = start.cluster,
                               stretch = stretch,
                               thresh = thresh)

    # Output
    if (output == "seurat") {
        # Add cell data
        seuratObj <- AddMetaDataMatrix(seuratObj,
                                       t(as.matrix(SummarizedExperiment::colData(sce.slingshot))))

        location <- paste("slingshot", reduction, sep = "_")
        seuratObj@tools[[location]] <-
          slingshot::SlingshotDataSet(sce.slingshot)
        return(seuratObj)
    } else if (output == "sce") {
        return(sce.slingshot)
    } else if (output == "slingshot") {
        return(slingshot::SlingshotDataSet(sce.slingshot))
    } else {
        stop("Select a valid output type: 'seurat', 'sce', 'slingshot'")
    }
}

AddMetaDataMatrix <- function(seuratObj, data_matrix) {
    for (i in 1:nrow(data_matrix)) {
        seuratObj[[rownames(data_matrix)[i]]] <- data_matrix[i,]
    }
    return(seuratObj)
}

sceObj <- RunSlingshot(seuratObj, output="sce")

seuratObj <- AddMetaDataMatrix(seuratObj,
              t(as.matrix(SummarizedExperiment::colData(sceObj))))
sds <-  slingshot::SlingshotDataSet(sceObj)


coord = Embeddings(object = seuratObj, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP_1", "UMAP_2")
coord = data.frame(ID = rownames(coord), coord)

meta = seuratObj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')

nms <- colnames(meta)[str_detect(colnames(meta), 'slingPseudotime')]

pal <- viridis(100, end = 0.95)
for (i in 1:length(nms)) {
    meta$point_color <- pal[cut(meta[,nms[i]], breaks = 100)]
    meta$point_color[is.na(meta$point_color)] <- "#EEEEEE"

    curveName <- paste0("curve", i)

    g <- ggplot() +
        geom_point(data = meta[meta$point_color!="#EEEEEE",], mapping = aes(x = UMAP_1, y = UMAP_2, color = point_color), size = 0.1) +
        geom_point(data = meta[meta$point_color=="#EEEEEE",], mapping = aes(x = UMAP_1, y = UMAP_2, color = point_color), size = 0.1) +
        geom_line(data = as.data.frame(sds@curves[[curveName]]$s, row.names = T), mapping = aes(x = UMAP_1, y = UMAP_2)) +
        scale_colour_identity() + ggtitle(nms[i]) + theme_classic()

    fileName <- file.path(opt$outDir, paste0("slingshot_", nms[i], ".pdf"))
    ggsave(g, filename = fileName)
}

Idents(seuratObj) <- seuratObj$seurat_clusters

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'vst', nfeatures = 2200)

top_hvg <- HVFInfo(seuratObj) %>% add_rownames(var = "bc") %>% arrange(desc(variance.standardized)) %>% top_n(300, variance.standardized) %>% pull(bc)

dat_use <- t(GetAssayData(seuratObj, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use)
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)

model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
    set_engine("ranger", importance = "impurity", num.threads = 3) %>%
    fit(pseudotime ~ ., data = dat_train)

val_results <- dat_val %>%
    mutate(estimate = predict(model, .[,-1]) %>% pull()) %>%
    select(truth = pseudotime, estimate)

var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:30]

pal <- viridis(100, end = 0.95)

gp.list=list()
for (i in seq_along(top_genes)) {
    meta$point_color <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]

    g <- ggplot() + geom_point(data = meta,
                               mapping = aes(x = UMAP_1,
                                             y = UMAP_2,
                                             color = point_color),
                               size = 0.1)
    for (i in 1:length(nms)) {
        curveName <- paste0("curve", i)
        g <- g + geom_line(data = as.data.frame(sds@curves[[curveName]]$s, row.names = T),
                           mapping = aes(x = UMAP_1, y = UMAP_2))
    }
    g <- g + scale_colour_identity() + ggtitle(top_genes[i]) + theme_classic()
    gp.list <- list.append(gp.list, g)
}

fileName <- file.path(opt$out, paste0("slingshot_DEGs_lineages.pdf"))

combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
pdf(fileName, width = 10, height = 10)
print(combined.gp)
dev.off()
