import pysam
import numpy as np

from bmkos.extract_barcode import load_whitelist, bc_align, add_bc_info, build_matrix
from bmkos.anno_gene import load_gtf, bam2bed, assign_gene
from bmkos.cluster_umi import get_umi, correct_umis
from bmkos.gene_expression import tags_bam, filter_tags_bam, get_expression


def qc(bam, info):
    qcd = dict()
    baminfo = pysam.AlignmentFile(bam, 'r')
    numMapped = baminfo.mapped
    numUnmapped = baminfo.unmapped
    baminfo.close()
    numReads = numMapped + numUnmapped
    qcd['Reads with Valid Barcode'] = '{:.2f}%'.format(info.is_keep.sum()/numReads*100)
    qcd['Reads Mapped to Genome'] = '{:.2f}%'.format(numMapped/numReads*100)
    qcd['Reads Mapped to Transcriptome'] = '{:.2f}%'.format((info.gene!='NA').sum()/numReads*100)
    qcd['Reads Mapped to Transcriptome And with Valid Barcode'] =  \
            '{:.2f}%'.format(np.logical_and(info.gene!='NA', info.is_keep).sum()/numReads*100)
    qcd['Full Length Rate'] = '{:.2f}%'.format(info.is_full_length.sum()/info.shape[0]*100)

    return qcd


def pipeline(bam, chrom, gtf, link1, read1, ssp, bm, outdir,
             bc1_list, bc2_list, bc3_list, bc1_kmer_idx,
             bc2_kmer_idx, bc3_kmer_idx, is_save_bam = False,
             max_read1_ed=6, max_link1_ed=6, min_align_score=200,
             match=5, mismatch = -1, gap_open=4, gap_extend=2, 
             acg_to_n_match=1, t_to_n_match=1):
    matrix = build_matrix(match, mismatch, acg_to_n_match, t_to_n_match)
    bam = pysam.AlignmentFile(bam, 'r')
    bed = bam2bed(bam.fetch(chrom))
    anno = assign_gene(bed, gtf)
    info = bc_align(bam.fetch(chrom), read1, link1, ssp, bm, matrix, gap_open, gap_extend)
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
                 'gene_name','umi_corr', 'ssp_score', 'is_full_length']
    info = info[keep_cols]
    #info.to_csv(path.join(tmpdir, chrom+'.'+'info.tsv'), sep='\t', header=True, index=True)
    #exp_raw.to_csv(path.join(tmpdir, chrom+'.'+'raw.matrix.tsv'), sep='\t', header=True, index=True)
    return info, exp_raw
