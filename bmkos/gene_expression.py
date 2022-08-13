import os
import pysam
import numpy as np
import pandas as pd
from os import path
from scipy import io, sparse

def tags_bam(bam, out_bam, info, chrom=None):
    bam = pysam.AlignmentFile(bam, 'r')
    chrbam = (chrom is None) and bam or bam.fetch(chrom)
    out_bam = pysam.AlignmentFile(out_bam, 'wb', template=bam)
    for i in chrbam:
        aread = info.loc[i.qname, :]
        tags = (
            ('R1', aread.bc1, 'Z'),
            ('R2', aread.bc2, 'Z'),
            ('R3', aread.bc3, 'Z'),
            ('C1', aread.bc1_corr, 'Z'),
            ('C2', aread.bc2_corr, 'Z'),
            ('C3', aread.bc3_corr, 'Z'),
            ('UR', aread.umi, 'Z'),
            ('UB', aread.umi_corr, 'Z'),
            ('GN', aread.gene, 'Z'),
        )
        i.set_tags(tags)
        out_bam.write(i)
    bam.close()
    out_bam.close()


def filter_tags_bam(bam, tags_bam, info):
    bam = pysam.AlignmentFile(bam, 'rb')
    tags_bam = pysam.AlignmentFile(tags_bam, 'wb', template=bam)
    for i in bam:
        aread = info.loc[i.qname, :]
        if not aread.is_keep or aread.gene=='NA':
            continue
        tags = (
            ('R1', aread.bc1, 'Z'),
            ('R2', aread.bc2, 'Z'),
            ('R3', aread.bc3, 'Z'),
            ('C1', aread.bc1_corr, 'Z'),
            ('C2', aread.bc2_corr, 'Z'),
            ('C3', aread.bc3_corr, 'Z'),
            ('UR', aread.umi, 'Z'),
            ('UB', aread.umi_corr, 'Z'),
            ('GN', aread.gene, 'Z'),
        )
        i.set_tags(tags)
        tags_bam.write(i)
    bam.close()
    tags_bam.close()

'''
def get_expression(cumi, gtf, outdir):
    exp_raw = cumi.groupby(['cell', 'gene'])['umi_corr'].nunique().reset_index()
    count_cell = exp_raw.groupby(['cell'])['umi_corr'].sum()
    count_cell = count_cell.astype(np.uint16)
    cells = count_cell.sort_values().index[-1000:].to_list()
    exp_raw = pd.pivot_table(exp_raw, index=['gene'], columns=['cell'], values='umi_corr', fill_value=0)
    exp_filter = exp_raw[cells]
    rawdir = path.join(outdir, 'raw_feature_bc_matrix')
    filterdir = path.join(outdir, 'filtered_feature_bc_matrix')
    path.isdir(rawdir) or os.mkdir(rawdir)
    path.isdir(filterdir) or os.mkdir(filterdir)
    rawm = sparse.csc_matrix(exp_raw, dtype=np.uint16)
    filterm = sparse.csc_matrix(exp_filter, dtype=np.uint16)
    io.mmwrite(path.join(rawdir, 'matrix.mtx'),
               rawm, 'BMK full-length RNA expression')
    io.mmwrite(path.join(filterdir, 'matrix.mtx'),
               filterm, 'BMK full-length RNA expression')
    gene_filter = gtf[['attribute', 'gene_name']].set_index('attribute').loc[exp_filter.index, :].reset_index()
    gene_filter['expression'] = 'Gene Expression'
    gene_filter.to_csv(path.join(filterdir, 'features.tsv'), sep='\t', index=False, header=False)
    exp_filter.columns.to_frame().to_csv(path.join(filterdir, 'barcodes.tsv'), sep='\t', header=False, index=False)

    gene_raw = gtf[['attribute', 'gene_name']].set_index('attribute').loc[exp_raw.index, :].reset_index()
    gene_raw['expression'] = 'Gene Expression'
    gene_raw.to_csv(path.join(rawdir, 'features.tsv'), sep='\t', index=False, header=False)
    exp_raw.columns.to_frame().to_csv(path.join(rawdir, 'barcodes.tsv'), sep='\t', header=False, index=False)
    return exp_raw, exp_filter
'''


def get_expression(cumi, gtf, outdir):
    exp_raw = cumi.groupby(['cell', 'gene'])['umi_corr'].nunique().reset_index()
    count_cell = exp_raw.groupby(['cell'])['umi_corr'].sum()
    #count_cell = count_cell.astype(np.uint16)
    cells = count_cell.sort_values().index[-5000:].to_list()
    exp_raw['umi_corr'] = exp_raw.umi_corr.astype(np.uint16)
    exp_raw.to_csv(path.join(outdir, 'highrawcount.tsv'), sep='\t', header=True, index=False)
    exp_raw = pd.pivot_table(exp_raw, index=['gene'], columns=['cell'], values='umi_corr', fill_value=0)
    exp_filter = exp_raw[cells]
    rawdir = path.join(outdir, 'raw_feature_bc_matrix')
    filterdir = path.join(outdir, 'filtered_feature_bc_matrix')
    path.isdir(rawdir) or os.mkdir(rawdir)
    path.isdir(filterdir) or os.mkdir(filterdir)
    rawm = sparse.csc_matrix(exp_raw, dtype=np.uint16)
    filterm = sparse.csc_matrix(exp_filter, dtype=np.uint16)
    io.mmwrite(path.join(rawdir, 'matrix.mtx'),
               rawm, 'BMK full-length RNA expression')
    io.mmwrite(path.join(filterdir, 'matrix.mtx'),
               filterm, 'BMK full-length RNA expression')
    gene_filter = gtf[['attribute', 'gene_name']].set_index('attribute').loc[exp_filter.index, :].reset_index()
    gene_filter['expression'] = 'Gene Expression'
    gene_filter.to_csv(path.join(filterdir, 'features.tsv'), sep='\t', index=False, header=False)
    exp_filter.columns.to_frame().to_csv(path.join(filterdir, 'barcodes.tsv'), sep='\t', header=False, index=False)

    gene_raw = gtf[['attribute', 'gene_name']].set_index('attribute').loc[exp_raw.index, :].reset_index()
    gene_raw['expression'] = 'Gene Expression'
    gene_raw.to_csv(path.join(rawdir, 'features.tsv'), sep='\t', index=False, header=False)
    exp_raw.columns.to_frame().to_csv(path.join(rawdir, 'barcodes.tsv'), sep='\t', header=False, index=False)
    return exp_raw, exp_filter
