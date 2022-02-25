#'--------------------------------------------------------------
#' filename : tcR.R
#' Date : 2021-01-15
#' contributor : Yanshuo Chu
#' function: tcR
#'--------------------------------------------------------------

print('<==== tcR ====>')

suppressMessages({
    library(optparse)
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
    library('tcR')
})

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'data.rds',
                metavar = 'character'),
    make_option(c("-t","--tcrFolder"),
                type = 'character',
                help = 'tcr folder',
                metavar = 'character'),
    make_option(c("-o","--out"),
                type = 'character',
                help = 'out',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


seuratObj <- readRDS(opt$data)
setwd(opt$out)

print(paste0("Loading data from ", opt$tcrFolder))
samples = list.dirs(path = opt$tcrFolder, recursive = F)


file.id = basename(samples)
sample.id <- file.id

if(length(file.id) < 1) stop(paste0('Failed to find data in ', path))


`%notin%` <- Negate(`%in%`)

tab=c()
tcr.data.list=list()
for(i in 1:length(file.id)){
    file = file.path(path, file.id[i], 'outs/clonotypes.csv')
    if(!file.exists(file)){next}
    tabi=read.csv(file)
    tcr.i = data.frame(Read.count = tabi$frequency,
                       Read.proportion = tabi$proportion,
                       CDR3.nucleotide.sequence = tabi$cdr3s_nt,
                       Umi.count = tabi$frequency,
                       Umi.proportion = tabi$proportion)

    tcr.i$Read.proportion=tcr.i$Read.count/sum(tcr.i$Read.count)
    tcr.i$Umi.proportion=tcr.i$Read.count/sum(tcr.i$Read.count)
    tcr.i$sample=sample.id[i]
    tcr.data.list[[length(tcr.data.list)+1]]=tcr.i

    tabi$proportion=tabi$frequency/sum(tabi$frequency)
    tabi$sample=sample.id[i]
    tab=rbind(tab,tabi)
}

names(tcr.data.list)=sample.id
# Productive Rearrangements (Observed Richness): The number of unique nucleotide
# rearrangements in the sample.
Richness = tab %>% group_by(sample) %>% summarise(richness = n())
# Shannon entropy: measure of entropy level of the system, defined as H =
# -sum(pi*logpi).
# Larger values indicates higher level of chaos (more even clonotypes in the sample)
Hp = tab %>% group_by(sample) %>% summarize(Hp = -sum(proportion*log2(proportion)))
Hpmax = tab %>% group_by(sample) %>% summarize(Hpmax = log2(length(unique(clonotype_id))))
Hp = inner_join(Hp, Hpmax , by = 'sample')
# Pielou Evenness( normalized Shannon entropy): Pielou’s Evenness (also known as
# the Shannon equitability) is a measure of how uniformly distributed the
# repertoire is , and it is computed as normalized Shannon’s Entropy. Values
# approaching 0 indicate a very skewed distribution of frequencies (i.e more
# variation in abundance) and values approaching 1 indicate that every
# rearrangement is present at nearly identical frequency (i.e. less variation in
# abundance). Pielou Clonality index reported in the Analyzer Sample Overview is
# defined as 1-Pielou’s evenness.
metrics = Hp %>% mutate(eveness = Hp/Hpmax, clonality = 1- eveness) #%>% arrange(desc(clonality))
# Simpson's D: Simpson’s D (also known as Simpson’s dominance index) is the sum
# over all observed rearrangements of the square fractional abundances of each
# rearrangement. This version of Simpson’s D is for an infinite population,
# termed lambda in the original reference. In this Analyzer tool, Simpson’s D is
# calculated using productive templates
simpsonD = tab %>%
    group_by(sample) %>%
    summarize(simpsonD = sum(proportion^2))
# Simpson Clonality: Simpson clonality is a method of quantifying the shape of a
# repertoire, ranging between 0 and 1, where values approaching 1 indicate a
# nearly monoclonal population. Simpson clonality is the square root of the sum
# over all observed rearrangements of the square fractional abundances of each
# rearrangement. Simpson clonality is also the square root of Simpson's D, and
# is robust across differences in sampling depths.
simpsonClonality = tab %>%
    group_by(sample) %>%
    summarize(simpsonClonality = sqrt(sum(proportion^2)))
metrics = left_join(metrics, Richness, by = 'sample')
metrics = left_join(metrics, simpsonD, by = 'sample')
metrics = left_join(metrics, simpsonClonality, by = 'sample')
##calculate other diverisy indices
# Simpson’s Diversity: 1)Simpson’s Diversity can be defined as the complement (1
# - D). Values of (1-D) Diversity range from 1 to 0, where values approaching 1
# correspond to a polyclonal, very diverse sample and values approaching 0
# correspond to a nearly monoclonal, non-diverse sample. 2)Simpsons Diversity
# can also be defined as the reciprocal (1/D). Values of (1/D) Diversity range
# from a minimum of 1 to a maximum of the richness (the number of unique
# nucleotide rearrangements in the sample). A diversity just above 1 corresponds
# to a nearly monoclonal sample and a diversity >>1, approaching the richness
# value, indicates a maximally diverse, polyclonal sample.
simpsonDiversity = repDiversity(tcr.data.list, 'inv.simp','read.prop')
shannonEntropy = sapply(tcr.data.list, function(x) tcR::entropy(x[['Read.proportion']], .norm = T))
# iChao1: A lower bound on repertoire richness, based on an improved
# nonparametric model of sample richness and abundance
chao1 = repDiversity(tcr.data.list, 'chao1')
chao1 = chao1[1,]
metrics$simpsonDiversity=simpsonDiversity[metrics$sample]
metrics$shannonEntropy = shannonEntropy[metrics$sample]
metrics$chao1 = chao1[metrics$sample]
metrics$simpsonEveness = metrics$simpsonD/metrics$richness
write.table(metrics,'TCR.diversity.txt',row.names = F,
            col.names = T,sep='\t',quote=F)
# overlap
twb.shared <- repOverlap(tcr.data.list, "exact", .norm = F, .verbose = F)
twb.shared=round(twb.shared,4)
p1<-vis.heatmap(twb.shared, .title = "shared clonotypes",
                .labs = c("", ""), .legend = "# clonotypes")
jaccard = repOverlap(tcr.data.list, .method = 'jaccard',.seq = 'nuc')
jaccard=round(jaccard,2)
p2<-vis.heatmap(jaccard, .title = "jaccard index",
                .labs = c("", ""), .legend = "")
moi = repOverlap(tcr.data.list, .method = 'morisita',.seq = 'nuc',
                 .quant = 'read.prop')
moi=round(moi,2)
p3<-vis.heatmap(moi, .title = "MOI index",
                .labs = c("", ""), .legend = "")

pdf('TCR.overlap.pdf',height=7,width=21)
grid.arrange(p1,p2,p3, nrow = 1, ncol = 3)
dev.off()


tab=c()
for(i in 1:length(file.id)){
    fcap_i <-  file.path(path, file.id[i], 'outs/filtered_contig_annotations.csv')
    if(!file.exists(fcap_i)){next}
    ctp_i <-  file.path(path, file.id[i], 'outs/clonotypes.csv')
    if(!file.exists(ctp_i)){next}

    fca_i <- read.csv(fcap_i) %>% select(barcode, raw_clonotype_id)
    ct_i <- read.csv(ctp_i) %>% select(clonotype_id, cdr3s_nt)
    fca_i$cdr3s_nt <-  ""
    fca_i$cdr3s_nt <- ct_i$cdr3s_nt[match(fca_i$raw_clonotype_id, ct_i$clonotype_id)]
    fca_i$barcode <- paste0(file.id[i], "_", str_remove_all(fca_i$barcode, "-[0-9]+$"))

    fca_i <- fca_i %>% select(barcode, cdr3s_nt)
    tab=rbind(tab, fca_i)
}

seuratObj@meta.data$cdr3s_nt <- NA
seuratObj@meta.data$cdr3s_nt <- tab$cdr3s_nt[match(rownames(seuratObj@meta.data), tab$barcode)]

###############################################################################
#'                        draw top 100 cclonotype bar                        '#

md <- seuratObj@meta.data
cn <- md %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "Blood") %>% group_by(cdr3s_nt) %>% dplyr::count() %>% arrange(desc(n))
cnTop100 <- cn[1:100,]

gdata <- md %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "Blood") %>%
    group_by(newClusterID, orig.ident, cdr3s_nt) %>% dplyr::count() %>%
    filter(cdr3s_nt %in% cnTop100$cdr3s_nt) %>% group_by(newClusterID, orig.ident) %>%
    dplyr::summarise(n = sum(n))


gdata$newClusterID <- as.factor(gdata$newClusterID)
gdata$orig.ident <- as.factor(gdata$orig.ident)

gdata$orig.ident <- factor(gdata$orig.ident, levels=c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"))


gdata$orig.ident <- mapvalues(gdata$orig.ident,
          from = c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"),
          to =   c("Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23",
                   "Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23"))

g <- ggplot(data = gdata) + geom_bar(stat = "identity", mapping = aes(x = newClusterID, y = n, fill = orig.ident)) + theme_classic()

pdf(file.path(figurePath, "Blood_patient_wise_barplot_100TopTCRClonotype.pdf"))
print(g)
dev.off()


## #' SF #########################################################################
## md <- seuratObj@meta.data
## cn <- md %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "SF") %>% group_by(cdr3s_nt) %>% count() %>% arrange(desc(n))

## cnTop100 <- cn[1:100,]


## gdata <- md %>% filter(!is.na(cdr3s_nt) & GROUP != "OA" & control == "SF") %>%
##     group_by(newClusterID, orig.ident, cdr3s_nt) %>% count() %>%
##     filter(cdr3s_nt %in% cnTop100$cdr3s_nt) %>% group_by(newClusterID, orig.ident) %>%
##     summarise(n = sum(n)) %>% ungroup() %>% add_row(newClusterID = 16, orig.ident = "76", n = 0 )

## gdata$newClusterID <- as.factor(gdata$newClusterID)
## gdata$orig.ident <- as.factor(gdata$orig.ident)

## gdata$orig.ident <- factor(gdata$orig.ident, levels=c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"))


## library(plyr)
## gdata$orig.ident <- mapvalues(gdata$orig.ident,
##           from = c("A7", "A40", "45", "A56", "62", "63", "65", "76", "A7B", "A40B", "45B", "A56B", "62B", "63B", "65B", "76B"),
##           to =   c("Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23",
##                    "Pt_3", "Pt_13", "Pt_16", "Pt_19", "Pt_20", "Pt_21", "Pt_22", "Pt_23"))

## g <- ggplot(data = gdata) + geom_bar(stat = "identity", mapping = aes(x = newClusterID, y = n, fill = orig.ident)) + theme_classic()

## pdf(file.path(figurePath, "SF_patient_wise_barplot_100TopTCRClonotype.pdf"))
## print(g)
## dev.off()

###############################################################################


cell.type.order <- c("CD69lo Naive CD4", "CD69hi Naive CD4", "Th1 and Th17-like", "Treg", "CXCL13+ T", "Naive CD8", "Central Memory CD8", "CXCR6+ Effector CD8", "CXCR6- Effector CD8", "CX3CR1+ Effector CD8", "CX3CR1- Effector CD8", "Cycling T", "MAIT T", "CX3CR1- Gamma-Delta T", "CX3CR1+ Gamma-Delta T", "CD16+ NK and NKT", "CD16- NK and NKT")

blood.id <- grep("B$", file.id, value = T)
SF.id <- setdiff(file.id, blood.id)
mdt <- as_tibble(seuratObj@meta.data)
stat <- mdt %>% group_by(orig.ident, control, GROUP) %>% count
Mono.id <- stat$orig.ident[stat$GROUP == "Mono"]
Combo.id <- stat$orig.ident[stat$GROUP == "Combo"]
OA.id <- stat$orig.ident[stat$GROUP == "OA"]
SF.id <- setdiff(SF.id, OA.id)
blood.id <- setdiff(blood.id, OA.id)
file.id <- setdiff(file.id, OA.id)

allIDsList <- list(all = file.id, blood = blood.id, SF = SF.id,
                   Mono = Mono.id, Combo = Combo.id, OA = OA.id)

for(idi in 1:length(allIDsList)){
    groupNames <- names(allIDsList[idi])
    groupSamples <- allIDsList[[idi]]

    md <- as.tibble(seuratObj@meta.data, rownames = NA) %>% mutate(rownames = rownames(seuratObj@meta.data))
    md <- md %>% arrange(factor(cell.types, levels = cell.type.order))

    allClonotypes <- unique(md$cdr3s_nt)

    ctn <- md %>% filter(!is.na(cdr3s_nt)) %>% group_by(cell.types, cdr3s_nt) %>% count
    ctn <- ctn %>% arrange(factor(cell.types, levels = cell.type.order))

    ctn.shared <- matrix(
        rep(NA, length(unique(ctn$cell.types)) ^ 2),
        length(unique(ctn$cell.types)),
        length(unique(ctn$cell.types)))

    p.ctn.shared <- matrix(
        rep(-1, length(unique(ctn$cell.types)) ^ 2),
        length(unique(ctn$cell.types)),
        length(unique(ctn$cell.types)))

    colnames(ctn.shared) <- unique(ctn$cell.types)
    rownames(ctn.shared) <- unique(ctn$cell.types)

    p.ctn.shared.tibble <- tibble(from = character(), to = character(), pvalue = double())
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            cti <- colnames(ctn.shared)[i]
            ctj <- colnames(ctn.shared)[j]
            cti_cts <- ctn %>% filter(cell.types == cti)
            ctj_cts <- ctn %>% filter(cell.types == ctj)

            ctij_cts <- intersect(cti_cts$cdr3s_nt, ctj_cts$cdr3s_nt)
###############################################################################
            #'                         draw pie chart                        '#
            tempPieData <- md %>% filter(cdr3s_nt %in% ctij_cts) %>% group_by(orig.ident) %>% count()
            if(dim(tempPieData)[1] > 0){

                g <- ggplot(tempPieData, aes(x = factor(""), y = n, fill = orig.ident)) +
                    geom_bar(position="fill", stat="identity") +
                    coord_polar("y") + ggtitle(paste0(cti, "_", ctj))

                tempPieDir <- file.path(figurePath, paste0(groupNames, "_matrix_pie"))
                if(!dir.exists(tempPieDir)){
                    dir.create(tempPieDir)
                }
                pdf(file.path(tempPieDir, paste0("matrix_", i, "_", j, "_cluster_pie.pdf")), width = 20, height = 20)
                print(g)
                dev.off()
            }
###############################################################################
            ctn.shared[i, j] <- length(ctij_cts)
            ctn.shared[j, i] <- ctn.shared[i, j]

            ijTest <- matrix(c(
                    length(ctij_cts),
                    length(ctj_cts$cdr3s_nt) - length(ctij_cts),
                    length(cti_cts$cdr3s_nt) - length(ctij_cts),
                    length(setdiff(allClonotypes, union(cti_cts$cdr3s_nt, ctj_cts$cdr3s_nt)))),
                    nrow = 2,
                    dimnames = list(iCluster = c("In", "NotIn"),
                                    jCluster = c("In", "NotIn")))

            ft <- fisher.test(ijTest, alternative = "greater")
            p.ctn.shared[i, j] <- ft[[1]]
            p.ctn.shared[j, i] <- ft[[1]]
            p.ctn.shared.tibble <- p.ctn.shared.tibble %>%
                add_row(from = colnames(ctn.shared)[i], to = colnames(ctn.shared)[j], pvalue = ft[[1]])
        }
    }

    p_values <- c()
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            p_values <- c(p_values, p.ctn.shared[i,j])
        }
    }
    p_values_adjust <- p.adjust(p_values, method = "BH")

    k = 1
    for(i in 1:(length(colnames(ctn.shared)) - 1)){
        for(j in (i + 1):length(colnames(ctn.shared))){
            p.ctn.shared[i,j] <- p_values_adjust[k]
            p.ctn.shared[j,i] <- p_values_adjust[k]

            p.ctn.shared.tibble$pvalue[p.ctn.shared.tibble$from == colnames(ctn.shared)[i] & p.ctn.shared.tibble$to == colnames(ctn.shared)[j]] == p_values_adjust[k]

            k <- k + 1
        }
    }

    p.ctn.shared.tibble <- p.ctn.shared.tibble %>% mutate(isAdjPSig = pvalue < 0.05, nlog10adjp=-log10(pvalue))

    p1<-vis.heatmap(ctn.shared, .title = "shared clonotypes",
                    .labs = c("", ""), .legend = "# clonotypes")
    tcrPlotFolder <- file.path(dirname(dataPath), "tcrPlot")
    if(!dir.exists(tcrPlotFolder)){
        dir.create(tcrPlotFolder)
    }

    write_tsv(as_tibble(p.ctn.shared), file.path(paste0("T_", groupNames, '_adj_p_value_sharedClonotypes.tsv')), col_names = F)
    write_tsv(as_tibble(p.ctn.shared.tibble), file.path(paste0("Edge_", groupNames, '_adj_p_value_sharedClonotypes.tsv')), col_names = T)

    pdf(file.path(paste0("T_", groupNames, "_sharedClonotypes_matrix.pdf")), width = 10, height = 10)
    print(p1)
    dev.off()
}


###############################################################################
#'                       mergeTCR infor into meta.data                       '#
###############################################################################

sampleNames <- c("S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58")
origIdents <- c("1053324_CTRNeg", "1053324_PEP3Pos", "1053324_PEP4Neg", "1053324_PEP6Pos", "1053324_PEP7Pos", "1053324_CEFPos", "1072015_CTRNeg", "1072015_PEP3Neg", "1072015_PEP6Pos")

sample2origT <- data.frame(sampleName = sampleNames, origIdent = origIdents)

seuratObj@meta.data$hasTCR = FALSE

TCRDataFolder <- "/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/RNASeqOut/TCR"
print(paste0("Loading data from ", TCRDataFolder))

samples = list.dirs(path = TCRDataFolder, recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ', TCRDataFolder))

obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)

    tempOrigIdent <- as.vector(sample2origT$origIdent[sample2origT$sampleName == ds])[1]

    sampleFolder = basename(list.dirs(path = file.path(TCRDataFolder, ds), recursive = F))[1]

    targetFile <- paste0(file.path(TCRDataFolder, ds, sampleFolder, "outs"), "/filtered_contig_annotations.csv")
    barcodeS <- read_csv(targetFile) %>% select(barcode)

    ## barcodeS$barcode <- paste0(ds,"_", str_remove_all(barcodeS$barcode, "[-1]"))

    seuratObj@meta.data$hasTCR[seuratObj@meta.data$orig.ident == tempOrigIdent &
                               str_remove(rownames(seuratObj@meta.data), "_[0-9]+$") %in% barcodeS$barcode] = TRUE
}
print(table(seuratObj@meta.data$hasTCR))

pdf(file.path(figurePath, "Dimplot_hasTCR.pdf"))
DimPlot(seuratObj, group.by='hasTCR', label=TRUE)
dev.off()



###############################################################################
#'                                TCR analysis                               '#
###############################################################################

rm(list=ls())
library(Seurat)
library(patchwork)
library(readr)
library(tidyverse)
library(rjson)
library(rlist)
library(SeuratData)
library(harmony)
library(Rmagic)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)

setwd('/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/scripts/figuresforpaper')
path='/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/RNASeqOut/TCR'
print(paste0("Loading data from ", path))
samples = list.dirs(path = path, recursive = F)
file.id = basename(samples)
sample.id <- file.id

if(length(file.id) < 1) stop(paste0('Failed to find data in ', path))

output=matrix(0,nrow=7,ncol=length(sample.id))
colnames(output)=sample.id
rownames(output)=c('N cells with TCR data',
                   'N cells overlapped with RNA-seq',
                   'Productive clonotype',
                   'TRA-TRB paired',
                   'TRA only',
                   'TRB only',
                   'Clonotype proportion range'
                   )

TCR.all.barcode.clonotype=c()

gp.list = list()
for (i in 1:length(file.id)) {
    print(file.id[i])
    file = file.path(path, file.id[i], 'outs/clonotypes.csv')
    if(!file.exists(file)){next}

    clonotype <- read.csv(file, stringsAsFactors = F)
    clonotype$clonotype_id2 = paste0(clonotype$clonotype_id, '_', file.id[i])

                                        #RNAseq=FetchData(subset(pbmc,subset=orig.ident==names(sample.id)[i]),vars=c('orig.ident'))
                                        #RNAseq$barcode=gsub('([ATCG]+_[NT0-9]+).*','\\1',rownames(RNAseq))

    filtered.contig <- read.csv(file.path(path, file.id[i], 'outs/filtered_contig_annotations.csv'), stringsAsFactors = F)

    filtered.contig$barcode2 = gsub('(.*)-1', paste0('\\1', '_', sample.id[i]), filtered.contig$barcode)
    filtered.contig$raw_clonotype_id2=paste0(filtered.contig$raw_clonotype_id,'_',sample.id[i])
    filtered.contig <- filtered.contig[,c('barcode2','raw_clonotype_id2')] %>% unique
                                        #filtered.contig<-left_join(filtered.contig,RNAseq,by=c("barcode2"="barcode"))
    filtered.contig <- left_join(filtered.contig,clonotype[,c('clonotype_id2','frequency','proportion','cdr3s_aa')],by=c('raw_clonotype_id2'='clonotype_id2'))
    filtered.contig$orig.ident <- sample.id[i]

    TCR.all.barcode.clonotype=rbind(TCR.all.barcode.clonotype, filtered.contig)


    output[1,i] = sum(clonotype$frequency) # cells with TCR data
    output[2,i] = NA  #sum(!is.na(filtered.contig$orig.ident)) # cells overlapped with RNA-seq
    output[3,i] = dim(clonotype)[1] # productive clonotype
    output[4,i] = sum(str_count(clonotype$cdr3s_aa,c('TRB:'))>0 & str_count(clonotype$cdr3s_aa,c('TRA:'))>0) # TRA and TRB paired
    output[5,i] = sum(str_count(clonotype$cdr3s_aa,c('TRB:'))==0) # TRA only
    output[6,i] = sum(str_count(clonotype$cdr3s_aa,c('TRA:'))==0) # TRB only
    output[7,i] = paste0(round(min(clonotype$proportion),4)*100,'%-',round(max(clonotype$proportion),4)*100,'%') # clonotype proportion range

                                        #bar plot
    data2<-as.data.frame(table(filtered.contig[,'raw_clonotype_id2']))
    colnames(data2)<-c('raw_clonotype_id2','Freq')
    data2=data2[rev(order(data2[,2])),]
    data2$proportion=data2$Freq/sum(data2$Freq)
    data2=data2[1:10,]

    plotx<-ggplot(data2, aes(x=reorder(raw_clonotype_id2,-proportion),y=proportion)) +
        geom_bar(stat="identity",position = 'stack',fill='#5777c4') +
        geom_text(aes(label=Freq),vjust=-1) +
        labs(y = 'proportion',x='',title=paste('Top 10 clonotypes in',sample.id[i])) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5),
              panel.background = element_blank(),axis.line = element_line(colour = "black"));
                                        #  ggsave(paste0('/Users/mdang1/Box/HN/figure/SingleCell/TCR/',sample.id[i],'.Top10.clonotype.pdf'),useDingbats=F,plotx,height=6,width = 6)
    gp.list[[length(gp.list) + 1]] = plotx
}

do.call(ggarrange, c(gp.list,ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('TCR.Top10.clonotype.',Sys.Date(),'.pdf'),height=7,width = 7)
print(combined.gp)
dev.off()

output=t(output)
output=data.frame('sample'=rownames(output),output,check.names = F)
write.table(output,'TCR_stat.txt',sep='\t',row.names = F,col.names = T,quote=F)

write.table(TCR.all.barcode.clonotype,'TCR.all.barcode.clonotype.txt',sep='\t',row.names = F,col.names = T,quote=F)


                                        # RNAseq=FetchData(pbmc,vars=c('orig.ident'))
                                        # RNAseq$barcode=gsub('([ATCG]+_[NT0-9]+).*','\\1',rownames(RNAseq))
                                        #
                                        # a<-left_join(RNAseq,TCR.all.barcode.clonotype.celltype,by=c('barcode'='barcode2'))
                                        #
                                        # pbmc[['TCR_clonotype']]<-factor(a$raw_clonotype_id2)
                                        # pbmc$TCR_clonotype_id=gsub('clonotype[0-9]+(_[NT0-9]+)','clonotype\\1',as.character(pbmc$TCR_clonotype))
                                        # pbmc[['TCR_clonotype_proportion']]<-a$proportion
                                        #
                                        # saveRDS(pbmc,'new.scs.rds')

########## clonotype deversity ##########

rm(list=ls())
library(Seurat)
library(patchwork)
library(readr)
library(tidyverse)
library(rjson)
library(rlist)
library(SeuratData)
library(harmony)
library(Rmagic)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)

setwd('/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/scripts/figuresforpaper')
path='/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/RNASeqOut/TCR'

print(paste0("Loading data from ", path))
samples = list.dirs(path = path, recursive = F)

file.id = basename(samples)
sample.id <- file.id

if(length(file.id) < 1) stop(paste0('Failed to find data in ', path))

library('tcR')

`%notin%` <- Negate(`%in%`)

tab=c()
tcr.data.list=list()
for(i in 1:length(file.id)){
    file = file.path(path, file.id[i], 'outs/clonotypes.csv')
    if(!file.exists(file)){next}
    tabi=read.csv(file)
    tcr.i = data.frame(Read.count = tabi$frequency,
                       Read.proportion = tabi$proportion,
                       CDR3.nucleotide.sequence = tabi$cdr3s_nt,
                       Umi.count = tabi$frequency,
                       Umi.proportion = tabi$proportion)

    tcr.i$Read.proportion=tcr.i$Read.count/sum(tcr.i$Read.count)
    tcr.i$Umi.proportion=tcr.i$Read.count/sum(tcr.i$Read.count)
    tcr.i$sample=sample.id[i]
    tcr.data.list[[length(tcr.data.list)+1]]=tcr.i

    tabi$proportion=tabi$frequency/sum(tabi$frequency)
    tabi$sample=sample.id[i]
    tab=rbind(tab,tabi)
}

names(tcr.data.list)=sample.id
# Productive Rearrangements (Observed Richness): The number of unique nucleotide
# rearrangements in the sample.
Richness = tab %>% group_by(sample) %>% summarise(richness = n())
# Shannon entropy: measure of entropy level of the system, defined as H =
# -sum(pi*logpi).
# Larger values indicates higher level of chaos (more even clonotypes in the sample)
Hp = tab %>% group_by(sample) %>% summarize(Hp = -sum(proportion*log2(proportion)))
Hpmax = tab %>% group_by(sample) %>% summarize(Hpmax = log2(length(unique(clonotype_id))))
Hp = inner_join(Hp, Hpmax , by = 'sample')
# Pielou Evenness( normalized Shannon entropy): Pielou’s Evenness (also known as
# the Shannon equitability) is a measure of how uniformly distributed the
# repertoire is , and it is computed as normalized Shannon’s Entropy. Values
# approaching 0 indicate a very skewed distribution of frequencies (i.e more
# variation in abundance) and values approaching 1 indicate that every
# rearrangement is present at nearly identical frequency (i.e. less variation in
# abundance). Pielou Clonality index reported in the Analyzer Sample Overview is
# defined as 1-Pielou’s evenness.
metrics = Hp %>% mutate(eveness = Hp/Hpmax, clonality = 1- eveness) #%>% arrange(desc(clonality))
# Simpson's D: Simpson’s D (also known as Simpson’s dominance index) is the sum
# over all observed rearrangements of the square fractional abundances of each
# rearrangement. This version of Simpson’s D is for an infinite population,
# termed lambda in the original reference. In this Analyzer tool, Simpson’s D is
# calculated using productive templates
simpsonD = tab %>%
    group_by(sample) %>%
    summarize(simpsonD = sum(proportion^2))
# Simpson Clonality: Simpson clonality is a method of quantifying the shape of a
# repertoire, ranging between 0 and 1, where values approaching 1 indicate a
# nearly monoclonal population. Simpson clonality is the square root of the sum
# over all observed rearrangements of the square fractional abundances of each
# rearrangement. Simpson clonality is also the square root of Simpson's D, and
# is robust across differences in sampling depths.
simpsonClonality = tab %>%
    group_by(sample) %>%
    summarize(simpsonClonality = sqrt(sum(proportion^2)))
metrics = left_join(metrics, Richness, by = 'sample')
metrics = left_join(metrics, simpsonD, by = 'sample')
metrics = left_join(metrics, simpsonClonality, by = 'sample')
##calculate other diverisy indices
# Simpson’s Diversity: 1)Simpson’s Diversity can be defined as the complement (1
# - D). Values of (1-D) Diversity range from 1 to 0, where values approaching 1
# correspond to a polyclonal, very diverse sample and values approaching 0
# correspond to a nearly monoclonal, non-diverse sample. 2)Simpsons Diversity
# can also be defined as the reciprocal (1/D). Values of (1/D) Diversity range
# from a minimum of 1 to a maximum of the richness (the number of unique
# nucleotide rearrangements in the sample). A diversity just above 1 corresponds
# to a nearly monoclonal sample and a diversity >>1, approaching the richness
# value, indicates a maximally diverse, polyclonal sample.
simpsonDiversity = repDiversity(tcr.data.list, 'inv.simp','read.prop')
shannonEntropy = sapply(tcr.data.list, function(x) tcR::entropy(x[['Read.proportion']], .norm = T))
# iChao1: A lower bound on repertoire richness, based on an improved
# nonparametric model of sample richness and abundance
chao1 = repDiversity(tcr.data.list, 'chao1')
chao1 = chao1[1,]
metrics$simpsonDiversity=simpsonDiversity[metrics$sample]
metrics$shannonEntropy = shannonEntropy[metrics$sample]
metrics$chao1 = chao1[metrics$sample]
metrics$simpsonEveness = metrics$simpsonD/metrics$richness
write.table(metrics,'TCR.diversity.txt',row.names = F,
            col.names = T,sep='\t',quote=F)
# overlap
twb.shared <- repOverlap(tcr.data.list, "exact", .norm = F, .verbose = F)
twb.shared=round(twb.shared,4)
p1<-vis.heatmap(twb.shared, .title = "shared clonotypes",
                .labs = c("", ""), .legend = "# clonotypes")
jaccard = repOverlap(tcr.data.list, .method = 'jaccard',.seq = 'nuc')
jaccard=round(jaccard,2)
p2<-vis.heatmap(jaccard, .title = "jaccard index",
                .labs = c("", ""), .legend = "")
moi = repOverlap(tcr.data.list, .method = 'morisita',.seq = 'nuc',
                 .quant = 'read.prop')
moi=round(moi,2)
p3<-vis.heatmap(moi, .title = "MOI index",
                .labs = c("", ""), .legend = "")

pdf('TCR.overlap.pdf',height=7,width=21)
grid.arrange(p1,p2,p3, nrow = 1, ncol = 3)
dev.off()

#' Get the shared clone's cells then get DEGs #############################


dataPath <- "/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/RNAanalysisTrim/nPC_30/UMAP_dist_0.1_nneighbor_50/p5_UMAP_dist_0.1_nneighbor_50_CLUSTER_res_0.5/cluster.rds"

seuratObj <- readRDS(dataPath)
Idents(seuratObj) <- seuratObj$orig.ident
metaData <- seuratObj@meta.data
metaData$barcode <- rownames(metaData)
metaData$barcode_origin <- str_remove_all(rownames(metaData), "_[0-9]+$")

TCRPath='/rsrch3/scratch/genomic_med/ychu2/data/P5-2001372/analysis/RNASeqOut/TCR'
for(si in 1:(length(sample.id) - 1)){
    for(sj in (si + 1):length(sample.id)){
        sample_i = sample.id[si]
        sample_j = sample.id[sj]
        tdl_i_l = intersectLogic(tcr.data.list[[sample_i]]$CDR3.nucleotide.sequence, tcr.data.list[[sample_j]]$CDR3.nucleotide.sequence)
        tdl_j_l = intersectLogic(tcr.data.list[[sample_j]]$CDR3.nucleotide.sequence, tcr.data.list[[sample_i]]$CDR3.nucleotide.sequence)
        shared_i_CDR3 = tcr.data.list[[sample_i]]$CDR3.nucleotide.sequence[tdl_i_l]
        shared_j_CDR3 = tcr.data.list[[sample_j]]$CDR3.nucleotide.sequence[tdl_j_l]
        shared_i_clonotypes = tab$clonotype_id[tab$sample == sample_i & tab$cdr3s_nt %in% shared_i_CDR3]
        shared_j_clonotypes = tab$clonotype_id[tab$sample == sample_j & tab$cdr3s_nt %in% shared_j_CDR3]

        filtered.contig_i = read_csv(file.path(TCRPath, sample_i, 'outs/filtered_contig_annotations.csv'))
        filtered.contig_j = read_csv(file.path(TCRPath, sample_j, 'outs/filtered_contig_annotations.csv'))

        shared_barcode_origin_i = filtered.contig_i$barcode[filtered.contig_i$raw_clonotype_id %in% shared_i_clonotypes]
        shared_barcode_origin_j = filtered.contig_j$barcode[filtered.contig_j$raw_clonotype_id %in% shared_j_clonotypes]

        shared_barcode_i = metaData$barcode[metaData$barcode_origin %in% shared_barcode_origin_i & metaData$orig.ident == sample_i]
        shared_barcode_j = metaData$barcode[metaData$barcode_origin %in% shared_barcode_origin_j & metaData$orig.ident == sample_j]

        if(length(shared_barcode_i) == 0 || length(shared_barcode_j) == 0){next}

        seuratObj@meta.data$tempIdent = ""
        seuratObj@meta.data$tempIdent[rownames(seuratObj@meta.data) %in% shared_barcode_i] = "i"
        seuratObj@meta.data$tempIdent[rownames(seuratObj@meta.data) %in% shared_barcode_j] = "j"
        Idents(seuratObj) <- seuratObj$tempIdent
        DEGs_iVSj <- FindMarkers(seuratObj, ident.1 = "i", ident.2 = "j")

        DEGs_iVSj$gene <- rownames(DEGs_iVSj)

        write_tsv(DEGs_iVSj, file.path(dirname(dataPath), paste0('DEGs_sharedTCR_', sample_i, '_VS_', sample_j, '.tsv')))
        print("finish")
        print(paste0('DEGs_', sample_i, '_VS_', sample_j, '.tsv'))
    }
}
