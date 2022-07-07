import os, sys, argparse, pickle
import pysam
import parasail
import numpy as np
import pandas as pd

from os import path
from tqdm import tqdm
from scipy import sparse, io
from concurrent.futures import ProcessPoolExecutor, as_completed, wait

sys.path.append(path.dirname(path.abspath(__file__)))
from extract_barcode import load_whitelist, bc_align, add_bc_info, build_matrix
from anno_gene import load_gtf, bam2bed, assign_gene
from cluster_umi import get_umi, correct_umis
from gene_expression import tags_bam, filter_tags_bam, get_expression

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
        help='gtf file'
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
        type=int, default=8,
        help="Max edit distance with the read1 sequence",
    )
    parser.add_argument(
        "--max_link1_ed",
        type=int, default=8,
        help="Max edit distance with the link1 sequence",
    )
    parser.add_argument(
        "--min_align_score",
        type=int, default=200,
        help="Min align score between bm  and reads",
    )
    args = parser.parse_args(args)
    return args

def pipeline(bam, chrom, gtf, link1, read1, bm, outdir,
             bc1_list, bc2_list, bc3_list, bc1_kmer_idx,
             bc2_kmer_idx, bc3_kmer_idx, is_save_bam = False,
             max_read1_ed=6, max_link1_ed=6, min_align_score=200,
             match=5, mismatch = -1, gap_open=4, gap_extend=2, 
             acg_to_n_match=1, t_to_n_match=1):
    tmpdir = path.join(outdir, 'tmp')
    path.isdir(tmpdir) or os.mkdir(tmpdir)
    matrix = build_matrix(match, mismatch, acg_to_n_match, t_to_n_match)
    bam = pysam.AlignmentFile(bam, 'r')
    bed = bam2bed(bam.fetch(chrom))
    anno = assign_gene(bed, gtf)
    info = bc_align(bam.fetch(chrom), read1, link1, bm, matrix, gap_open, gap_extend)
    bam.close()
    is_keep = np.logical_and(info.score > min_align_score,
                             np.logical_or(info.editread1 <= max_read1_ed,
                                           info.editlink1 <= max_link1_ed
                                           )
                             )
    is_keep = np.logical_and(is_keep, info.bc1.str.count('[ATCG]') >= 5)
    is_keep = np.logical_and(is_keep, info.bc2.str.count('[ATCG]') >= 5)
    is_keep = np.logical_and(is_keep, info.bc3.str.count('[ATCG]') >= 5)
    info['is_keep'] = is_keep
    info, bc1_info = add_bc_info(info, bc1_list, bc1_kmer_idx, 'bc1')
    info, bc2_info = add_bc_info(info, bc2_list, bc2_kmer_idx, 'bc2')
    info, bc3_info = add_bc_info(info, bc3_list, bc3_kmer_idx, 'bc3')
    info[['gene', 'gene_name']] = anno[['gene', 'gene_name']]
    # 校正UMI，增加tag到bam文件中
    info['umi'] = info.apply(get_umi, axis=1, args=(matrix,))
    #cumi = info[np.logical_and(info.is_keep, info.gene != 'NA')][
    #    ['bc1', 'bc2', 'bc3', 'bc1_corr', 'bc2_corr', 'bc3_corr', 'umi', 'gene']]
    #cumi['cell'] = cumi['bc1_corr'] + cumi['bc2_corr'] + cumi['bc3_corr']
    cumi = info[np.logical_and(info.is_keep, info.gene != 'NA')][
           ['bc1_corr_idx', 'bc2_corr_idx', 'bc3_corr_idx', 'umi', 'gene']]
    bc_idx = ['bc1_corr_idx', 'bc2_corr_idx', 'bc3_corr_idx']
    cumi[bc_idx] = cumi[bc_idx].astype('str')
    cumi['cell'] = 'bc1_' + cumi.bc1_corr_idx + '-bc2_' + cumi.bc2_corr_idx + '-bc3_' + cumi.bc3_corr_idx

    cumi['gene_cell'] = cumi['cell'] + ':' + cumi['gene']
    counts_dict = dict(cumi.umi.value_counts())
    cumi["umi_corr"] = cumi.groupby(["gene_cell"])["umi"].transform(correct_umis)
    info['umi_corr'] = cumi.loc[:, 'umi_corr']
    info.umi_corr.fillna('N', inplace=True)
    info[['bc1', 'bc2', 'bc3']] = info[['bc1', 'bc2', 'bc3']].apply(lambda a: a.str.replace('-', ''), axis=0)
    if is_save_bam:
        tags_bam(bam, path.join(outdir, chrom+'.bmk_tags.sorted.bam'), info)
        filter_tags_bam(bam, path.join(outdir, chrom+'.bmk_filter_tags.sorted.bam'), info)
    # 计算表达矩阵
    exp_raw = cumi.groupby(['cell', 'gene'])['umi_corr'].nunique().reset_index()
    #exp_raw.to_csv(path.join(tmpdir, chrom+'.'+'cumi.tsv'), sep='\t', header=True, index=True)
    #exp_raw = pd.pivot_table(exp_raw, index=['gene'], columns=['cell'], values='umi_corr', fill_value=0)
    keep_cols = ['score','edit', 'editread1','editlink1', 'is_keep',
                 'bc1_match_ed','bc2_match_ed','bc3_match_ed','gene',
                 'gene_name','umi_corr']
    info = info[keep_cols]
    #info.to_csv(path.join(tmpdir, chrom+'.'+'info.tsv'), sep='\t', header=True, index=True)
    #exp_raw.to_csv(path.join(tmpdir, chrom+'.'+'raw.matrix.tsv'), sep='\t', header=True, index=True)
    return info, exp_raw

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
    with open(args.bc1) as f:
        bc1_list = sorted([i.strip() for i in f])
    with open(args.bc2) as f:
        bc2_list = sorted([i.strip() for i in f])
    with open(args.bc3) as f:
        bc3_list = sorted([i.strip() for i in f])
    bc1_kmer_idx = load_whitelist(bc1_list, k=5)
    bc2_kmer_idx = load_whitelist(bc2_list, k=5)
    bc3_kmer_idx = load_whitelist(bc3_list, k=5)
    bc1_len = len(bc1_list[0])
    bc2_len = len(bc2_list[0])
    bc3_len = len(bc3_list[0])
    link1 = args.link1
    read1 = args.read1[-args.read1_length:]
    bm = (read1 + 'N' * bc1_len +
          link1 + 'N' * bc2_len + 'ACGACTC' +
          'N' * (bc3_len + args.umi_length) +
          'T' * args.polyT_length)
    print(bm)
    gtf = load_gtf(args.gtf)

    p = ProcessPoolExecutor(args.threads)
    baminfo = pysam.AlignmentFile(args.bam, 'r')
    chroms = baminfo.references
    baminfo.close()

    tasks = [p.submit(pipeline, args.bam, i, gtf, link1, read1, bm, args.outdir,
             bc1_list, bc2_list, bc3_list, bc1_kmer_idx,
             bc2_kmer_idx, bc3_kmer_idx, is_save_bam = False,
             max_read1_ed=6, max_link1_ed=6, min_align_score=200,
             match=5, mismatch = -1, gap_open=args.gap_open, gap_extend=args.gap_extend,
             acg_to_n_match=1, t_to_n_match=1) for i in chroms]
    pbar = tqdm(total=len(tasks), desc='process chorms:')
    for i in as_completed(tasks):
        pbar.update()
    wait(tasks)
    info = [i.result()[0] for i in tasks]
    cumi = [i.result()[1] for i in tasks]
    cumi = pd.concat(cumi, axis=0)
    barcodes = cumi.cell.unique()
    genes = gtf.attribute.unique()
    cumi.set_index('cell', inplace=True)
    p = ProcessPoolExecutor(args.threads)
    tasks = [p.submit(each_exp, cumi, barcodes, genes, i, 1000) for i in tqdm(range(0, barcodes.shape[0], 1000))]
    pbar = tqdm(total=len(tasks), desc='process expression matrix:')
    for i in as_completed(tasks):
        pbar.update()
    wait(tasks)

    rawdir = path.join(args.outdir, 'raw_feature_bc_matrix')
    filterdir = path.join(args.outdir, 'filtered_feature_bc_matrix')
    mkdir_or_not([rawdir, filterdir])
    with open(path.join(rawdir, 'features.tsv'), 'w') as f:
        for i in genes:
            f.write(i+'\n')
    with open(path.join(rawdir, 'barcodes.tsv'), 'w') as f:
        for i in barcodes:
            f.write(i+'\n')
    tmp = [i.result() for i in tasks]
    with open(path.join(args.outdir, 'cache_raw_expression.pkl'), 'wb') as f:
        pickle.dump(tmp, f)

    raw_exp = sparse.hstack(tmp)
    io.mmwrite(path.join(rawdir, 'matrix.mtx'), raw_exp, 'BMK full-length RNA expression')


if __name__=='__main__':
    main()



    










