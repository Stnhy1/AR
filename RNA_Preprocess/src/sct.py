#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: scrublettrail.py
#        Author: Chu Yanshuo
#         Email: yanshuochu@qq.com
# =============================================================================
'''

import gzip
import scipy.io
import numpy as np
import scrublet as scr
import matplotlib.pyplot as plt


def load_barcodes(filename):
    """override the load_genes function in scrublet help funcs"""
    barcode_list=[]
    with gzip.open(filename, 'r') as f:
        for l in f:
            l = l.decode('utf8').strip()
            barcode_list.append(l)
    return barcode_list


def predict_doublet(input_dir, out_file):
    """get predict result of scrublet

    input_dir = '/rsrch3/scratch/genomic_med/ychu2/data/tmp/DP-RNA_10X_SC/analysis/outs/5GE_CC/65B/outs/filtered_feature_bc_matrix'

    export barcode for each sample. so to embed this result in to R's seurat package
    """

    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
    # genes = np.array(load_genes(input_dir + '/features.tsv.gz', delimiter='\t', column=1))
    barcodes = np.array(load_barcodes(input_dir + '/barcodes.tsv.gz'))

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                              min_cells=3,
                                                              min_gene_variability_pctl=85,
                                                              n_prin_comps=30)


    with open(input_dir + '/' + out_file, 'w') as outf:
        for barcode, isdt, scoredt in zip(barcodes, predicted_doublets, doublet_scores):
            outf.write("{2}\t{0}\t{1}\n".format(isdt, scoredt, barcode))
            pass
        pass
    pass

    predicted_doublets.tofile(input_dir + '/' + out_file + '_predicted_doublets.txt',sep='\t')
    doublet_scores.tofile(input_dir + '/' + out_file + '_doublet_scores.txt',sep='\t')
    print('savefile Done.')

    print('start Plot...')
    scrub.plot_histogram()
    plt.savefig(input_dir + '/' + out_file +'_scrub.plot_histogram.pdf')
    plt.close('all')
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(input_dir + '/' + out_file +'_scrub.plot_UMAP.pdf')
    plt.close('all')


def main():
    """parse options """
    import argparse
    parser = argparse.ArgumentParser(description='scrublet script')
    parser.add_argument('--input_dir', dest='input_dir', help='input put dir')
    parser.add_argument('--out_file', dest='out_file', help='out file name')
    args = parser.parse_args()

    predict_doublet(args.input_dir, args.out_file)

if __name__ == '__main__':
    main()

