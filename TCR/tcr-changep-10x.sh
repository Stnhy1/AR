module load singularity
module load cellranger/3.1.0

srcD='/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/scripts/src'
runR="Rscript --no-save "

DIST=0.15
NPROC=12

TCR_CC_DIR=/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/TCR_CC
LOCAL_BIN=/rsrch3/home/genomic_med/ychu2/.local/bin
IMGT_VDJ_REF_DIR=/rsrch3/home/genomic_med/ychu2/software/share/germlines/imgt/human/vdj
LOCAL_BLAST_IGDATA_DIR=/rsrch3/home/genomic_med/ychu2/.conda/pkgs/igblast-1.15.0-h18cd00f_0/share/igblast/bin

for SAMPLE_NAME in $(ls ${TCR_CC_DIR}); do
    TEMP_SAMPLE_DIR=/data/${SAMPLE_NAME}/outs
    READS=${TEMP_SAMPLE_DIR}/filtered_contig.fasta
    ANNOTATIONS=${TEMP_SAMPLE_DIR}/filtered_contig_annotations.csv
    OUT_DIR=${TEMP_SAMPLE_DIR}/changeo-10x-outs

    singularity exec -B ${TCR_CC_DIR}:/data ${HOME}/software/immcantation/immcantation-3.0.0.sif \
                changeo-10x -s $READS -a $ANNOTATIONS \
                -g human -t tr -x $DIST \
                -n $SAMPLE_NAME -o $OUT_DIR -f airr -p $NPROC

    CC_OUT_DIR=${TCR_CC_DIR}/${SAMPLE_NAME}/outs
    cd $CC_OUT_DIR
    ${LOCAL_BIN}/AssignGenes.py igblast -s filtered_contig.fasta -b ~/share/igblast \
                       --organism human --loci ig --format blast

    ${LOCAL_BIN}/MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta \
                  -r IMGT_Human_*.fasta --10x filtered_contig_annotations.csv --extended

    CHANGEO_OUT_DIR=${TCR_CC_DIR}/${SAMPLE_NAME}/outs/changeo-10x-outs
    cd ${CHANGEO_OUT_DIR}

    ${LOCAL_BIN}/DefineClones.py -d ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_db-pass.tsv --act set --model ham \
                        --norm len --dist 0.16
    ${LOCAL_BIN}/CreateGermlines.py -d ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_db-pass_clone-pass.tab \
                    -g dmask --cloned \
                    -r ${IMGT_VDJ_REF_DIR}/imgt_human_IGHV.fasta \
                    ${IMGT_VDJ_REF_DIR}/imgt_human_IGHD.fasta \
                    ${IMGT_VDJ_REF_DIR}/imgt_human_IGHJ.fasta

    ${LOCAL_BIN}/DefineClones.py -d ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T.tab --act set --model ham \
                --norm len --dist 0.16
    ${LOCAL_BIN}/CreateGermlines.py -d ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T_clone-pass.tab \
                -g dmask --cloned \
                -r ${IMGT_VDJ_REF_DIR}/imgt_human_IGHV.fasta \
                ${IMGT_VDJ_REF_DIR}/imgt_human_IGHD.fasta \
                ${IMGT_VDJ_REF_DIR}/imgt_human_IGHJ.fasta
    ${runR} ${srcD}/createTreeObj.r \
            --outputDir ${CHANGEO_OUT_DIR} \
            --germPassTable ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T_clone-pass_germ-pass.tab \
            --cloneOutputFile ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T_clone.rds \
            --dnapars_exec ${HOME}/software/phylip-3.697/exe/dnapars \
            --subDFOutputFile ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T_subDF.rds \
            --treeOutputFile ${CHANGEO_OUT_DIR}/${SAMPLE_NAME}_heavy_FUNCTIONAL-T_treeVIS.rds
done
