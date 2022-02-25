rm /rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/p2rhp2.o /rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/p2rhp2.e -rf

module load python/3.7.3-anaconda
module load R/3.6.0

srcD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/src'
DataD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_CC'
##setup enviroment variables
ResD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_VI/harmony'
paramD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/params'
databaseD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/database'

runR="Rscript --no-save "


${runR} ${srcD}/load-cellranger.R -d $DataD -o $ResD/combined.data.rds -c ${paramD}/load-cellranger.json

### get doublet infor from scrublet ################################################
/usr/bin/bash /rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/rmdoublets.harmony.pipeline.sh

${runR} ${srcD}/visualize_sbt.R -d ${ResD}/combined.data.rds -o ${ResD}/doublet.pdf -c ${paramD}/qc-plotSize.json

${runR} ${srcD}/qc-doublets-viability.R -d ${ResD}/combined.data.rds -o ${ResD}/doublets-viability-plot.pdf -c ${paramD}/qc-plotSize.json
${runR} ${srcD}/cell-cycle.R -d ${ResD}/combined.data.rds -o ${ResD}/cell-cycle.tsv -s human 
${runR} ${srcD}/determinePC.R -d ${ResD}/combined.data.rds -o ${ResD}/elBowPlot.pdf
${runR} ${srcD}/filter-normalize.R -d $ResD/combined.data.rds -o $ResD/combined.filtered.data.rds -c ${paramD}/filter-normalize.json

########################## harmony
${runR} ${srcD}/addBatchInfo.R -d ${ResD}/combined.filtered.data.rds -o ${ResD}/combined.filtered.data.withbatch.rds
${runR} ${srcD}/run-harmony.R -d ${ResD}/combined.filtered.data.withbatch.rds -o ${ResD}/harmony.rds

###############################################################################
#      Run magic (after scale in harmany)
###############################################################################
${runR} ${srcD}/magicimpute.R -d ${ResD}/combined.filtered.data.rds -o ${ResD}/magic.rds

###############################################################################
#                               subset operation
###############################################################################

# /usr/bin/bash /rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/M.sub.pipeline.sh
# /usr/bin/bash /rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/T.sub.pipeline.sh

##########################################################################################################################################
# here to find best npc(elboplot determine pc) and dist (only view on umap) and res(affected to cluster number 20)
# rm -rf ${ResD}/umapfinder/
# mkdir -p ${ResD}/umapfinder/
# ${runR} ${srcD}/snn-harmony-best-umap-finder.R -d ${ResD}/combined.filtered.data.withbatch.rds -o ${ResD}/umapfinder/umap
# ${runR} ${srcD}/determinePC-harmony.R -d ${ResD}/harmony.rds -o ${ResD}/elBowPlot-harmony.pdf
##########################################################################################################################################

# ${runR} ${srcD}/umap-harmony.R -d ${ResD}/harmony.rds -o ${ResD}/snn-harmony-umap.rds -c ${paramD}/snn-harmony-umap.json
# ${runR} ${srcD}/addBatchInfo.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/snn-harmony-umap.rds

###############################################################################
#                              umap visualization                             #
###############################################################################

#${runR} ${srcD}/addBatchInfo.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/snn-harmony-umap.rds
# ${runR} ${srcD}/qc-doublets-viability-by-cluster.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/doublets-viability-by-cluster-plot.pdf -c ${paramD}/qc-plotSize.json

# ${runR} ${srcD}/visualize.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/umap.pdf -c ${paramD}/visualize_umap.json
# ${runR} ${srcD}/visualize_batch.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/umap_bach.pdf -c ${paramD}/visualize_umap_batch.json

# ${runR} ${srcD}/snn-marker.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/snn-markers.tsv -c ${paramD}/snn-marker.json
# ${runR} ${srcD}/snn-heatmap.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/markersHeatmap.pdf -c ${paramD}/snn-heatmap.json -m ${ResD}/snn-markers.tsv

###############################################################################
#                               Further analysis                              #
###############################################################################

# ${runR} ${srcD}/visualize_batch.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/umap_bach.pdf -c ${paramD}/visualize_umap_batch.json

###################################################################################

# python ${srcD}/statMarkers.py --markersTop ${ResD}/markers.top.tsv --markersDatabase ${databaseD}/immuneMarkers_minghaoandenyu.txt --out ${ResD}/markers.top_mhd_celltype.tsv
# python ${srcD}/statMarkers.py --markersTop ${ResD}/markers.top.tsv --markersDatabase ${databaseD}/immuneCellMarkerAllinBox_Yanshuo.txt --out ${ResD}/markers.top_ys_celltype.tsv
# python ${srcD}/statMarkers.py --markersTop ${ResD}/snn-markers.tsv --markersDatabase ${databaseD}/immuneMarkers_minghaoandenyu.txt --out ${ResD}/markers.mhd_celltype.tsv
# python ${srcD}/statMarkers.py --markersTop ${ResD}/snn-markers.tsv --markersDatabase ${databaseD}/immuneCellMarkerAllinBox_Yanshuo.txt --out ${ResD}/markers.ys_celltype.tsv


# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/robust -c ${paramD}/feature-plot-origin.json -m ${databaseD}/robustMarker.txt
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/immuneMarkersenyu -c ${paramD}/feature-plot-origin.json -m ${databaseD}/immuneMarkersenyu.txt
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/iCMAPurified -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/ppt -c ${paramD}/feature-plot-origin.json -m ${databaseD}/ppt.txt
# 
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/robust_split_control -c ${paramD}/feature-plot-origin.json -m ${databaseD}/robustMarker.txt --control
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/immuneMarkersenyu_split_control -c ${paramD}/feature-plot-origin.json -m ${databaseD}/immuneMarkersenyu.txt --control
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/iCMAPurified_split_control -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt --control
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/ppt_split_control -c ${paramD}/feature-plot-origin.json -m ${databaseD}/ppt.txt --control

# ${runR} ${srcD}/ridge-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/ridge-plot.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/immuneMarkersenyu.txt
# ${runR} ${srcD}/dot-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/dot-plot.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/immuneMarkersenyu.txt
# ${runR} ${srcD}/heatmap-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/heatmap-plot.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/immuneMarkersenyu.txt
# ${runR} ${srcD}/feature-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/feature-plot_yp.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt
# ${runR} ${srcD}/ridge-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/ridge-plot_yp.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt
# ${runR} ${srcD}/dot-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/dot-plot_yp.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt
# ${runR} ${srcD}/heatmap-plot.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/heatmap-plot_yp.pdf -c ${paramD}/feature-plot-origin.json -m ${databaseD}/iCMA_Y_purified.txt

# ${runR} ${srcD}/feature-plot-topmarker.R -d ${ResD}/snn-harmony-umap.rds -o ${ResD}/featureplot/topmarker -c ${paramD}/feature-plot-topmarker.json -m ${ResD}/markers.top_ys_celltype.tsv


##########################

# ${runR} "${srcD}"/sc3.R -d "${ResD}"/umap.rds --out-prefix "${ResD}"/sc3 -c "${srcD}"/param-sc3.json

##############################################################################################################################################################

###############################################################################
#                           result files management                           #
###############################################################################

# shopt -s extglob nullglob

# ls ~/data/tmp/outs/$(basename ${ResD})/*.@(pdf|tsv|txt) | xargs -n 1 rm -f
# ls ${ResD}/*.@(pdf|txt|tsv) | xargs -n 1 -I{} cp {} ~/data/tmp/outs/$(basename ${ResD})/

# rm -rf ~/data/tmp/outs/$(basename ${ResD})/featureplot
# cp -r ${ResD}/featureplot ~/data/tmp/outs/$(basename ${ResD})/


