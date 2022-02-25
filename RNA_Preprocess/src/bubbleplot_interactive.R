#' filename : bubbleplot_interactive.R
#' Date : 2020-06-24
#' contributor : Yanshuo Chu
#' function: bubbleplot_interactive

rm(list=ls())

library(readr)
library(rjson)
library(rlist)
library(SeuratData)
library(harmony)
library(Rmagic)
library(Seurat)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(viridis)

figD <- "/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/tempPlots"
seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/snn-harmony-umap.rds")

## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/SubC6/dist_0.2_resolution_0.7/umap.rds")
## seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/umap-harmony-origin.rds")

checkpoint.pattern.list <- list(
    "ADORA2"=grep("^ADORA2", rownames(seuratObj@assays$RNA@data), value = T),
    "BTLA"=grep("^BTLA$", rownames(seuratObj@assays$RNA@data), value = T),
    "VISTA"=grep("^(VIS(R|TA)|C10orf54)$", rownames(seuratObj@assays$RNA@data), value = T),
    "CD200"=grep("^CD200([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "TNFRSF"=grep("^TNFRSF(4|9|18|14|25|12)([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "OX40L|CD137L"=grep("^TNFSF(4|9)([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD27"=grep("^CD27([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "PDL-1"=grep("^CD274([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD276"=grep("^CD276([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD28"=grep("^CD28([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD40"=grep("^CD40([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD70"=grep("^CD70([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD80"=grep("^CD80([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD86"=grep("^CD86([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CEACAM1"=grep("^CEACAM1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CTLA4"=grep("^CTLA4([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "GAL3"=grep("^LGALS3([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "HAVCR|TIM-3"=grep("^HAVCR(1|2)([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "ICOS"=grep("^ICOS(LG$|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "IDO1"=grep("^IDO1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "IL2RB"=grep("^IL2RB$", rownames(seuratObj@assays$RNA@data), value = T),
    "KIR3DL1"=grep("^KIR3DL1$", rownames(seuratObj@assays$RNA@data), value = T),
    "LAG3"=grep("^LAG3([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "LAIRI"=grep("^LAIRI([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "PD1"=grep("^PDCD1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "PVR"=grep("^PVR([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "SLAM"=grep("^SLAMF1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "TIGIT"=grep("^TIGIT([0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "VTCN1"=grep("^(VTCN1([^0-9]|$)|B7[^0-9]|VCTN1)", rownames(seuratObj@assays$RNA@data), value = T)
)

gene.pattern.list <- list(
    "CD1C"=grep("^CD1C", rownames(seuratObj@assays$RNA@data), value = T),
    "CD3"=grep("^CD3([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD4"=grep("^CD4([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD8"=grep("^CD8([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD11"=grep("^ITGA", rownames(seuratObj@assays$RNA@data), value = T),
    "CD14"=grep("^CD14", rownames(seuratObj@assays$RNA@data), value = T),
    "CD16"=grep("^FCGR3", rownames(seuratObj@assays$RNA@data), value = T),
    "CD19"=grep("^CD19([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD20"=grep("^MS4A", rownames(seuratObj@assays$RNA@data), value = T),
    "CD38"=grep("^CD38([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD56"=grep("^NCAM1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD79"=grep("^CD79([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CD138"=grep("^SDC1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CXCR3"=grep("^CXCR3([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "CX3CR1"=grep("^CX3CR1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "S100A8"=grep("^S100A8([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "S100A9"=grep("^S100A9([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T),
    "TRGV"=grep("^TRGV", rownames(seuratObj@assays$RNA@data), value = T),
    "TRDV"=grep("^TRDV", rownames(seuratObj@assays$RNA@data), value = T),
    "TRAV"=grep("^TRAV", rownames(seuratObj@assays$RNA@data), value = T),
    "TRBV"=grep("^TRBV", rownames(seuratObj@assays$RNA@data), value = T),
    "HLA-D"=grep("^HLA-D", rownames(seuratObj@assays$RNA@data), value = T),
    "ITGAX"=grep("^ITGAX", rownames(seuratObj@assays$RNA@data), value = T),
    "PLA2G2A"=grep("^PLA2G2A", rownames(seuratObj@assays$RNA@data), value = T),
    "PDCD1"=grep("^PDCD1([^0-9]|$)", rownames(seuratObj@assays$RNA@data), value = T)
)

Tgene<-c('LEF1','CCR7','SELL','IL7R','TCF7','COTL1','CD4','CD40LG','S100A4',
        'FOXP3','IL2RA','IKZF2',#Treg
        'CTLA4','TNFRSF4','LAYN','BATF', 'LAG3','TIGIT','PDCD1','HAVCR2',#Inhibitory receptors
        "IL2","GZMA","GZMB",'GZMK','GZMH','GNLY','PRF1',"IFNG",'KLRF1','KLRD1','NCAM1','NKG7', #Cytokines_effectors
        'CD8A','CD8B','CD28',"TNFRSF14","ICOS","TNFRSF9", #CoStimulatory
        "EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX", #TFs
        'CD69','MAF','HLF','IL10RA','KLRG1','DDIT4','CLIC1','FCER1G','XCL2','CD3D','CST7','CCL5',
        "CXCR3","CXCR4","LTB","LYAR","CD44",'ENTPD1' #other
        )
# Th
Thgene<-c('STAT4','CXCR3','TBX21','IFNG', #Th1
        'STAT6','GATA3','CCR4','PTGDR2','IL4','IL5','IL13', #Th2
        'CCR6','KLRB1','IL17A','IL17F', #Th17
        'STAT5','FOXP3','IL2RA','IL7R','IL10','TGFB1',#Treg
        'STAT3','BCL6','CXCR5','IL21','CXCL13','CCR7' # Tfh
        )
#NK
NKgene <- c("GZMA","GZMB","GZMK","PRF1","GNLY", # Granules
                "CD2","SELL","NCAM1","FCGR3A","KLRC1","KLRD1",'KLRB1', # surface proteins
                "KIR2DL1","KIR3DL1","KIR2DL3","KIR3DL2","KIR3DL3", # inhibitory KIR
                "ACTB","ARPC3","ARPC4","CFL1","PFN1","CST7",#Negative regulators
                "XCL1","XCL2","PVRIG","TIGIT","CCL5")

RPgene<-c('RPLP2','RPS21','RPS12','RPS27A','RPS10','RPS23','RPS15A','RPS25','RPS25','RPL32')

gene.group.list <- list(
    "T_gamma|delta|alpha|beta" = c(
        gene.pattern.list[["TRGV"]],
        gene.pattern.list[["TRDV"]],
        gene.pattern.list[["TRAV"]],
        gene.pattern.list[["TRBV"]])
    ## "T" = c(Tgene),
    ## "Th" = c(Thgene),
    ## "NK" = c(NKgene),
    ## "RP" = c(RPgene),
    ## "B_Plasmablast_Plasma" = c(
    ##     gene.pattern.list[["CD19"]],
    ##     gene.pattern.list[["CD20"]],
    ##     gene.pattern.list[["CD38"]]
    ## ),
    ## "B_T_NK_Myeloid_Else" = c(
    ##     gene.pattern.list[["CD79"]],
    ##     gene.pattern.list[["CD3"]],
    ##     gene.pattern.list[["CD68"]],
    ##     gene.pattern.list[["LYZ"]],
    ##     gene.pattern.list[["CD14"]],
    ##     gene.pattern.list[["S100A8"]],
    ##     gene.pattern.list[["S100A9"]],
    ##     gene.pattern.list[["PLA2G2A"]]
    ## ),
    ## "DCs" = c(
    ##     gene.pattern.list[["HLA-D"]],
    ##     gene.pattern.list[["ITGAX"]],
    ##     gene.pattern.list[["CD1C"]]
    ## ),
    ## "Non-Lineage" = c(
    ##     gene.pattern.list[["CXCR3"]],
    ##     gene.pattern.list[["PDCD1"]],
    ##     gene.pattern.list[["CX3CR1"]]
    ## )
)


for(i in 1:length(gene.pattern.list)){
    pattern.name <- names(gene.pattern.list[i])
    gene.list <- gene.pattern.list[[i]]

    tempMarkerList <- gene.list

    print("marker list: ")
    print(tempMarkerList)
    if (length(tempMarkerList) < 1) {
        next
    }

    gp.list=list()
    for(i in 1: length(tempMarkerList)){
        gp.list = list.append(gp.list, annotate_figure(p = FeaturePlot(seuratObj, reduction = "umap", features = as.vector(tempMarkerList[i]), cols=c("lightgray", "blue", "black")), top = ggpubr::text_grob(label = pattern.name, face="bold", size = 20, color="red")))
    }

    combined.gp <- do.call(ggarrange, c(gp.list, ncol = 2, nrow = 2))
    pdf(paste0(figD, "/featureplot_", pattern.name, ".pdf"), width = 18, height = 12)
    print(combined.gp)
    dev.off()
}



for(i in 1:length(gene.group.list)){
    gene.list.name <- names(gene.group.list[i])
    gene.list <- gene.group.list[[i]]

    gene<-intersect(gene.list,rownames(seuratObj))
    p<-DotPlot(seuratObj, features = gene)
    data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
    #data[data$avg.exp>7,'avg.exp']<-7
    plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
      geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
      #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
      scale_fill_viridis()+
      scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
      theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_text(angle = 90,vjust = 1,hjust = 1) )
    ht <- 8/31 * length(unique(data$id))
    wh <- 4/10 * length(gene.list)
    if(ht < 5){ht = 5;}
    if(wh < 5){wh = 5;}
    ggsave(paste(figD, '/bubbleplot_', gene.list.name, '.pdf', sep=''), plotx,
           height=ht,
           width=wh)
}


seuratObj <- readRDS("/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony_RM_doublets/SubT/umap-harmony-origin.rds")

gene.list <- c()
for(item in checkpoint.pattern.list){gene.list <- c(item, gene.list)}
gene.list.name <- "SubT_checkpoints"

gene<-intersect(gene.list,rownames(seuratObj))
p<-DotPlot(seuratObj, features = gene)
data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
#data[data$avg.exp>7,'avg.exp']<-7
plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_viridis()+
  scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+             ## to tune the size of circles
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust = 1,hjust = 1) )
ht <- 8/31 * length(unique(data$id))
wh <- 4/10 * length(gene.list)
if(ht < 2){ht = 2;}
if(wh < 5){wh = 5;}
ggsave(paste(figD, '/bubbleplot_', gene.list.name, '.pdf', sep=''), plotx,
        height=ht,
        width=wh)
