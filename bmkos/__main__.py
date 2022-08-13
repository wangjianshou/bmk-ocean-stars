import os, sys, pickle
import argparse
import gzip
import pysam
import parasail
import numpy as np
import pandas as pd

from os import path
from tqdm import tqdm
from scipy import sparse
from scipy.io import mmwrite
from io import BytesIO

from concurrent.futures import ProcessPoolExecutor, as_completed, wait

#sys.path.append(path.dirname(path.abspath(__file__)))
from bmkos.pipeline import pipeline, qc, process_unmap,generate_report
from bmkos.extract_barcode import load_whitelist, bc_align, add_bc_info, build_matrix
from bmkos.anno_gene import load_gtf, bam2bed, assign_gene
from bmkos.cluster_umi import get_umi, correct_umis
from bmkos.gene_expression import tags_bam, filter_tags_bam, get_expression
from bmkos.plotQC import plotLen, plotQscore

#read结构: read1-bc1-link1-bc2-ACGACTC-bc3-umi-polyT-CDS-ssp

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
            '--bam', '-b', 
            required=True,
            help='sorted bam file'
    )
    parser.add_argument(
        '--gtf', '-g',
        required=True,
        help='gtf file'
    )
    parser.add_argument(
        '--outdir', '-o',
        required=False, default='.',
        help='output direction'
    )
    parser.add_argument(
            '--threads', '-th',
            required=False, type=int, default=16,
            help=''
    )
    parser.add_argument(
        '--save-bam', '-sb',
        action='store_true',
        help='save tags bam or not'
    )
    parser.add_argument(
            '--product', '-pd',
            required=False, default='cell', choices=['cell', 'spatial'],
            help=''
    )
    parser.add_argument(
        '--debug', '-dg', action='store_true',
        help='debug mode output more: '
    )
    parser.add_argument(
            '--expect-cells', '-ec', default=3000, required=False,
            type=int, help='number of expected cells'
    )
    parser.add_argument(
            '--sampleid', '-id', default="SampleID", required=False,
             help='sample id in web summary'
    )
    parser.add_argument(
            '--link1', '-l1',
            required=False, default='GTCATCGCAGAGTACTACGT',
            help=''
    )
    parser.add_argument(
            '--read1', '-r1',
            required=False, default='CTACACGACGCTCTTCCGATCT',
            help=''
    )
    parser.add_argument(
            '--ssp', '-ssp',
            required=False, default='CCCAGCAATATCAGCACCAACAGAAA',
            help=''
    )
    parser.add_argument(
            '--read1_length', '-r1len',
            required=False, type=int, default=16,
            help=''
    )
    parser.add_argument(
            '--bc1', '-b1',
            required=False, default='sequence/barcode/bc1.txt',
            help=''
    )
    parser.add_argument('--bc2', '-b2',
            required=False, default='sequence/barcode/bc2.txt',
            help=''
    )
    parser.add_argument(
            '--bc3', '-b3',
            required=False, default='sequence/barcode/bc3.txt',
            help=''
    )
    parser.add_argument(
            '--umi_length', '-ul',
            required=False, type=int, default=12,
            help=''
    )
    parser.add_argument(
            '--polyT_length', '-T',
            required=False, type=int, default=10,
            help=''
    )
    parser.add_argument(
            "-go", "--gap_open",
            help="Gap open penalty",
            type=int, default=4
    )
    parser.add_argument(
            "-e", "--gap_extend",
            help="Gap extend penalty",
            type=int, default=2
    )
    parser.add_argument(
            "-m", "--match",
            type=int, default=5,
            help="Match score",
    )
    parser.add_argument(
            "-x", "--mismatch",
            help="Mismatch score",
            type=int, default=-1
    )
    parser.add_argument(
        "-n", "--acg_to_n_match",
        help="Score for A/C/G<-->N match",
        type=int, default=1,
    )
    parser.add_argument(
        "-s", "--t_to_n_match",
        help="Score for T<-->N match",
        type=int, default=1
    )
    parser.add_argument(
        "--max_read1_ed",
        type=int, default=6,
        help="Max edit distance with the read1 sequence",
    )
    parser.add_argument(
        "--max_link1_ed",
        type=int, default=4,
        help="Max edit distance with the link1 sequence",
    )
    parser.add_argument(
        "--min_align_score",
        type=int, default=200,
        help="Min align score between bm  and reads",
    )
    parser.add_argument(
        "--mstep", '-ms',
        type = int, default = 1000,
        help = "process Matrix each step",
    )
    args = parser.parse_args(args)
    return args
def each_exp(cumi, barcodes, genes, start, step=1000):
    x = cumi.loc[barcodes[start:start+step], :]
    raw_exp = pd.pivot_table(x, index=['gene'], columns=['cell'], values='umi_corr', fill_value=0)
    raw_exp = raw_exp.reindex(index = genes, ).fillna(0, downcast='int32')
    raw_exp = sparse.coo_matrix(raw_exp)
    return raw_exp

def mkdir_or_not(dirs):
    for i in dirs:
        path.isdir(i) or os.mkdir(i)
    return dirs




def main():
    args = parseArgs()
    rawdir = path.join(args.outdir, 'raw_feature_bc_matrix')
    filterdir = path.join(args.outdir, 'filtered_feature_bc_matrix')
    mkdir_or_not([args.outdir, rawdir, filterdir])
    with open(args.bc1) as f:
        bc1_list = [i.strip() for i in f]
    with open(args.bc2) as f:
        bc2_list = [i.strip() for i in f]
    with open(args.bc3) as f:
        bc3_list = [i.strip() for i in f]
    bc1_kmer_idx = load_whitelist(bc1_list, k=5)
    bc2_kmer_idx = load_whitelist(bc2_list, k=5)
    bc3_kmer_idx = load_whitelist(bc3_list, k=5)
    bc1_len = len(bc1_list[0])
    bc2_len = len(bc2_list[0])
    bc3_len = len(bc3_list[0])
    link1 = args.link1
    read1 = args.read1[-args.read1_length:]
    bm = (read1 + 'N' * bc1_len +
          link1 + 'N' * bc2_len + 'CGACTC' +
          'N' * (bc3_len + args.umi_length) +
          'T' * args.polyT_length)
    gtf = load_gtf(args.gtf)

    p = ProcessPoolExecutor(args.threads)
    baminfo = pysam.AlignmentFile(args.bam, 'r')
    #numReads = baminfo.mapped + baminfo.unmapped
    chroms = baminfo.references
    baminfo.close()

    tasks = [
                p.submit(pipeline, args.bam, i, gtf, link1, read1, args.ssp, bm, args.outdir,
                         bc1_list, bc2_list, bc3_list, bc1_kmer_idx,
                         bc2_kmer_idx, bc3_kmer_idx, is_save_bam = args.save_bam,
                         max_read1_ed=args.max_read1_ed, max_link1_ed=args.max_link1_ed,
                         min_align_score=args.min_align_score, match=args.match,
                         mismatch = args.mismatch, gap_open=args.gap_open, gap_extend=args.gap_extend,
                         acg_to_n_match=args.acg_to_n_match, t_to_n_match=args.t_to_n_match)
                for i in chroms
            ]
    pbar = tqdm(total=len(tasks), desc='process chorms:')
    for i in as_completed(tasks):
        pbar.update()

    matrix = build_matrix(args.match, args.mismatch, args.acg_to_n_match, args.t_to_n_match)
    bam = pysam.AlignmentFile(args.bam, 'r')
    unmap_info = bc_align(bam.fetch('*'), read1, link1, args.ssp, bm, matrix)
    unmap_info = process_unmap(unmap_info, bc1_list, bc2_list, bc3_list,
                  bc1_kmer_idx, bc2_kmer_idx, bc3_kmer_idx,
                  args.min_align_score, args.max_read1_ed, args.max_link1_ed)
    bam.close()

    wait(tasks)
    info = pd.concat([i.result()[0] for i in tasks], axis=0)
    info = pd.concat([info, unmap_info], axis=0)
    if args.debug:
        info.to_csv(path.join(args.outdir, 'reads_info.tsv'), index=True, header=True, sep='\t')
    qcd = qc(info)

    #with open(path.join(args.outdir, 'summary_qc.txt'), 'w') as f:
    #    for i in qcd.keys():
    #        f.write(i+'\t'+qcd[i]+'\n')

    cumi = [i.result()[1] for i in tasks]
    cumi = pd.concat(cumi, axis=0)
    barcodes = cumi.cell.unique()
    genes = gtf.attribute.unique()
    cumi.set_index('cell', inplace=True)
    if args.debug:
        cumi.to_csv(path.join(args.outdir, 'reads_cumi.tsv'), index=True, header=True, sep='\t')
    p = ProcessPoolExecutor(args.threads)
    #tasks = [p.submit(each_exp, cumi, barcodes, genes, i, 1000) for i in tqdm(range(0, barcodes.shape[0], 1000))]
    tasks = [p.submit(each_exp, cumi, barcodes, genes, i, args.mstep) for i in
             tqdm(range(0, barcodes.shape[0], args.mstep))]
    pbar = tqdm(total=len(tasks), desc='process expression matrix:')
    for i in as_completed(tasks):
        pbar.update()
    wait(tasks)

    features = gtf[['attribute', 'gene_name']].set_index('attribute').loc[genes, :].reset_index()
    features['exp'] = "Gene Expression"
    with gzip.open(path.join(rawdir, 'features.tsv.gz'), 'wb') as f:
        with BytesIO() as tmpf:
            features.to_csv(tmpf, header=False, index=False, sep='\t')
            f.write(tmpf.getvalue())

    with gzip.open(path.join(rawdir, 'barcodes.tsv.gz'), 'wt') as f:
        for i in barcodes:
            f.write(i+'\n')
    tmp = [i.result() for i in tasks]
    if args.debug:
        with open(path.join(args.outdir, 'cache_raw_expression.pkl'), 'wb') as f:
            pickle.dump(tmp, f)

    raw_exp = sparse.hstack(tmp)
    with gzip.open(path.join(rawdir, 'matrix.mtx.gz'), 'wb') as f:
        with BytesIO() as tmpf:
            mmwrite(tmpf, raw_exp, 'BMK full-length RNA expression')
            f.write(tmpf.getvalue())
    if args.product == 'cell':
        from bmkcc import cells
        matrix, filtered_matrix = cells.callcells(rawdir, args.outdir, args.expect_cells)
        cell_barcodes = filtered_matrix.bcs[:]
        counts_per_bc = matrix.get_counts_per_bc()
        raw_barcodes = matrix.bcs[:]
        cellfig =cells.barcode_umi_plot(counts_per_bc, raw_barcodes, cell_barcodes)
        qcd['cellfig'] = cellfig.to_json()
        cellqc = cells.cellQC(filtered_matrix)
        qcd.update(cellqc)

    qcd['ReadsLengthDistribution'] = plotLen(info.qlen).to_json()
    qcd['ReadsQscoreDistribution'] = plotQscore(info.qscore/info.qlen).to_json()
    qcd['Smapleid'] = args.sampleid
    with open(path.join(args.outdir, 'web_summary.html'), 'w') as f:
        f.write(generate_report(qcd, args.product))



if __name__=='__main__':
    main()



    










