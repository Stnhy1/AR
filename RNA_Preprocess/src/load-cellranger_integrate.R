suppressMessages({library(optparse)
library(ggplot2)
library(readr)
library(rjson)
library(Seurat)})
print('----loading cellranger outputs----')
##CLI parsing
option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'dataFolder',
                metavar = 'character'),

    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character'),

    make_option(c("-c",'--param'),
                type = 'character',
                help = 'json parameter file for plot figure size',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

param <- fromJSON(file = opt$param)

## Get all input samples
print(paste0("Loading data from ", opt$data))
samples = list.dirs(path = opt$data,recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ',opt$data))

## Get obj list
obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]

    print(ds)

    tenx.data = Read10X(file.path(opt$data, ds, "outs", "filtered_feature_bc_matrix"))

    tenx0 = CreateSeuratObject(counts = tenx.data,
                               min.cells = param$min.cells,
                               min.features = param$min.features,
                               project = ds)

    ###############################################################################
    #              QC subset              #
    ###############################################################################

    mito.features = grep(pattern = '^MT-|^mt-', x = rownames(x = tenx0), value = T)
    percent.mito <- Matrix::colSums(x = GetAssayData(object = tenx0, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = tenx0, slot = 'counts'))
    tenx0[['percent.mito']] = percent.mito
    tenx0 <- subset(x = tenx0,
                    subset = nFeature_RNA > param$nfeature_RNA_min & nFeature_RNA < param$nfeature_RNA_max & percent.mito < param$percent.mito)

    ###############################################################################
    #            Normalize Data           #
    ###############################################################################
    print('log normalization')
    tenx0 <- NormalizeData(object = tenx0)

    ###############################################################################
    #        Find Variable Features       #
    ###############################################################################
    tenx0 <- FindVariableFeatures(tenx0,
                                  selection.method = "vst")

    obj.list[[length(obj.list)+1]] = tenx0
}

names(obj.list) <- samples
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:param$intMaxDim)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:param$intMaxDim)

DefaultAssay(obj.integrated) <- "integrated"

print("saving output")
saveRDS(obj.integrated, file = opt$out)
print('----end----')
