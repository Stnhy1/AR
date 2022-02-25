#'--------------------------------------------------------------
#' filename : pathway-enrichment.R
#' Date : 2020-11-21
#' contributor : Yanshuo Chu
#' function: pathway-enrichment
#'--------------------------------------------------------------

print('<==== pathway-enrichment ====>')

suppressMessages({
    library(optparse)
    library(fgsea)
    library(data.table)
    library(ggplot2)
    library(Seurat)
    library(tidyverse)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

option_list = list(
    make_option(c("-t","--DEGtable"),
                type = 'character',
                help = 'DEG table',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


DEG_T <- read_tsv(opt$DEGtable)
## p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene
figurePath <- file.path(dirname(opt$DEGtable), "DEG-PathwayEnrichment")
if(!dir.exists(figurePath)){
    dir.create(figurePath)
}

## DEG_T <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/project6/result/DEG/featureCounts/human/protein_coding_gene/sampleInfofeatureCounts_protein_coding_gene_countMatrixDESeq2_result.tsv")
## Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	log10fdr
## figurePath <- opt$out

#' read hall mark gene pathways ###############################################

HallMarkGenePath <- "/rsrch3/home/genomic_med/ychu2/data/hallmarkGeneSets/h.all.v7.1.entrez.gmt"
hpw <- gmtPathways(HallMarkGenePath)

#' read DEGs table ############################################################

IDMap <- mapIds(org.Hs.eg.db, keys = DEG_T$gene, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")

DEG_T$ENTREZID <- IDMap[DEG_T$gene]

for(temp_c in unique(DEG_T$cluster)){
    upGenes = DEG_T %>% filter(cluster == temp_c, avg_logFC > 0) %>% arrange(desc(avg_logFC)) %>%  pull(gene)
    downGenes = DEG_T %>% filter(cluster == temp_c, avg_logFC < 0) %>% arrange(avg_logFC) %>% pull(gene)
    allGeneList <- list()
    if(length(upGenes) > 0){
        allGeneList[["up"]] <- upGenes
    }
    if(length(downGenes) > 0){
        allGeneList[["down"]] <- downGenes
    }

    for(i in 1:length(allGeneList)){
        tempName <- names(allGeneList[i])
        tempGenes <- allGeneList[[i]]

        temp_DEG_T <- DEG_T %>% filter((!is.na(ENTREZID))) %>%
            filter(gene %in% tempGenes) %>% mutate(log2FoldChange = abs(avg_logFC)) %>% arrange(desc(log2FoldChange))

        markerRank <- temp_DEG_T$log2FoldChange
        names(markerRank) <- temp_DEG_T$ENTREZID

        fgseaRes <- fgsea(pathways = hpw,
                          stats    = markerRank,
                          minSize  = 15,
                          maxSize  = 1000,
                          nperm = 1000)

        topPathwaysUp <- fgseaRes[ES > 0][head(order(NES, decreasing = T), n=10), pathway]

        ## topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        ## topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        topPathways <- c(topPathwaysUp)
        pdf(file.path(figurePath, paste0("cluster", temp_c, "_", tempName, "_pathwayAnalysis.pdf")), width  = 15)
        plotGseaTable(hpw[topPathways], markerRank, fgseaRes, gseaParam=0.5)
        dev.off()
    }
}
