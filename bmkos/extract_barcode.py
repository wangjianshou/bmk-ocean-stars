import parasail
import numpy as np
import pandas as pd
import editdistance as ed
import pysam


def build_matrix(match, mismatch, acg_to_n_match, t_to_n_match):
    matrix = parasail.matrix_create("ACGTN", match, mismatch)

    ############################
    # SCORING MATRIX POSITIONS #
    ############################
    #     A   C   G   T   N
    # A   0   1   2   3   4   5
    # C   6   7   8   9   10  11
    # G   12  13  14  15  16  17
    # T   18  19  20  21  22  23
    # N   24  25  26  27  28  29

    # Update scoring matrix so that N matches A/C/G/N
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = t_to_n_match
    matrix.pointer[0].matrix[27] = t_to_n_match
    return matrix




def split_seq_into_kmers(seq, k=5):
    #assert len(seq) >= k, "Pick a value for k that is less than len(barcode)"
    if len(seq) < k:
        return [seq, ]
    return [seq[i : i + k] for i in range(0, len(seq) - k + 1)]


def load_whitelist(wl, k=5):
    kmer_to_bc_index = {}
    for index, bc in enumerate(wl):
        bc_kmers = split_seq_into_kmers(bc, k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = {index}
            else:
                kmer_to_bc_index[bc_kmer].add(index)
    return kmer_to_bc_index


def filter_whitelist_by_kmers(wl, kmers, kmer_to_bc_index):
    # collect sets of indices that each kmer points to
    id_sets = [
        kmer_to_bc_index[kmer] for kmer in kmers if kmer in kmer_to_bc_index.keys()
    ]
    # retain all barcodes that have at least one kmer match with the query barcode
    all_filt_indices = list(set().union(*id_sets))
    filt_wl = [wl[i] for i in all_filt_indices]
    return filt_wl, all_filt_indices


def get_cor_bc(uncor_bc, bc_list, bc_kmer_idx, k=5):
    uncor_kmers = split_seq_into_kmers(uncor_bc, k)
    filt_bc_list,filt_bc_idx = filter_whitelist_by_kmers(bc_list, uncor_kmers, bc_kmer_idx)
    if not filt_bc_list:
        return -1, 'N', 100, 0
    if len(filt_bc_list)==1:
        return filt_bc_idx[0]+1, filt_bc_list[0], ed.eval(filt_bc_list[0], uncor_bc), 100
    ed_list = [ed.eval(uncor_bc, i) for i in filt_bc_list]
    ax = np.argpartition(ed_list, 1)
    #bc_match_idx = ed_list[ax[0]] < ed_list[ax[1]] and ax[0]  or ax[1]
    bc_match_ed = ed_list[ax[0]]
    bc_match = filt_bc_list[ax[0]]
    next_match_diff = ed_list[ax[1]]-ed_list[ax[0]]
    return filt_bc_idx[ax[0]]+1, bc_match, bc_match_ed, next_match_diff


def add_bc_info(df, bc_list, bc_kmer_idx, name):
    df = df[:]
    match_ed = name+'_match_ed'
    match_ed_diff = name + '_match_ed_diff'
    bc_corr = name + '_corr'
    bc_corr_idx = bc_corr + '_idx'
    x = df[name].str.replace('-', '').apply(get_cor_bc, True, (bc_list, bc_kmer_idx))
    bc_info = pd.DataFrame(x.to_list(), columns=[bc_corr_idx, bc_corr, match_ed, match_ed_diff], index=x.index)
    bc_info = bc_info.astype({bc_corr_idx:np.int32, match_ed:np.int32, match_ed_diff:np.int32})
    df = pd.concat([df, bc_info], axis=1)
    #keepreads = keepreads[np.logical_and(keepreads[match_ed] <= 6, keepreads[match_ed_diff] >= 1)]
    #keepreads.index = range(keepreads.shape[0])
    is_keep = np.logical_and(df[match_ed] <= 6, df[match_ed_diff] >= 1)
    df['is_keep'] = np.logical_and(is_keep, df.is_keep)
    return df, bc_info



#barcode的组合序列bm在比对过程中会出现 有N被clip的read设置为极端情况，
# 目前测试的数据这种情况比较少，因此会直接去除，也就是barcode和umi会被设置为N
# 后续会实现更优方案：设置被clip的N的阈值，低于3个N，可以正常提取
'''
def bc_align(sam, read1, link1, bm, matrix, gap_open=2, gap_extend=4):
    #sam = pysam.AlignmentFile(sam, 'r')
    tmp = []
    for i in sam:
        seq = i.get_forward_sequence()[:250]
        rseq = seq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]
        frd = parasail.sw_trace(seq, bm, gap_open, gap_extend, matrix)
        rev = parasail.sw_trace(rseq, bm, gap_open, gap_extend, matrix)
        align = (frd.score >= rev.score) and frd or rev
        nidx = [i for i in range(len(align.traceback.ref)) if align.traceback.ref[i]=='N']
        if len(nidx) == 17+19+19+12:
            r1 = align.traceback.query[:nidx[0]]
            l1 = align.traceback.query[nidx[16]+1:nidx[17]]
            bc1 = align.traceback.query[nidx[0]:nidx[16]+1]
            bc2 = align.traceback.query[nidx[17]:nidx[35]+1]
            bc3 = align.traceback.query[nidx[36]:nidx[54]+1]
            rawumi = align.traceback.query[nidx[55]:nidx[66]+1]
        else:  # 有N被clip的read设置为极端情况，后续会考虑其他方案
            r1 = l1 = 'N' * 20
            bc1, bc2, bc3, rawumi = 'N', 'N', 'N', 'N'
        editread1 = ed.eval(read1, r1)
        editlink1 = ed.eval(link1, l1)
        edit = ed.eval(align.traceback.query, align.traceback.ref)
        #tmp.append([i.query_name,  bc1, bc2, bc3, rawumi, i.is_mapped, align.score, align.traceback.comp,
        #           align.traceback.ref, align.traceback.query,
        #           align.query, frd.score, rev.score,
        #           edit, editread1, editlink1,])
        tmp.append([i.query_name,  bc1, bc2, bc3, rawumi, i.is_mapped, align.score,
                    align.query, edit, editread1, editlink1])

    #keep_cols = ['id', 'bc1', 'bc2', 'bc3', 'rawumi', 'is_mapped', 'score', 'comp',
    #             'ref', 'query', 'rawquery', 'frdscore', 'revscore', 'edit',
    #             'editread1', 'editlink1']
    keep_cols = ['id', 'bc1', 'bc2', 'bc3', 'rawumi', 'is_mapped', 'score',
                 'rawquery', 'edit', 'editread1', 'editlink1']
    info = pd.DataFrame(tmp, columns=keep_cols)
    info.set_index('id', inplace=True)
    #sam.close()
    return info
'''

def bc_align(sam, read1, link1, ssp, bm, matrix, gap_open=2, gap_extend=4):
    #sam = pysam.AlignmentFile(sam, 'r')
    tmp = []
    for i in sam:
        readseq = i.get_forward_sequence()
        revreadseq = readseq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]
        
        seq = readseq[:250]
        rseq = readseq[-250:]
        rseq = rseq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]
        frd = parasail.sw_trace(seq, bm, gap_open, gap_extend, matrix)
        rev = parasail.sw_trace(rseq, bm, gap_open, gap_extend, matrix)
        
        frd_ssp_align = parasail.sw_trace(readseq, ssp, gap_open, gap_extend, matrix)
        rev_ssp_align = parasail.sw_trace(revreadseq, ssp, gap_open, gap_extend, matrix)
        #align = (frd.score >= rev.score) and frd or rev
        if frd.score >= rev.score:
            align = frd
            direction = '+'
        else:
            align = rev
            ssp_direction = '-'
        if frd_ssp_align.score >= rev_ssp_align.score:
            ssp_align = frd_ssp_align
            ssp_direction = '+'
        else:
            ssp_align = rev_ssp_align
            ssp_direction = '-'
        nidx = [i for i in range(len(align.traceback.ref)) if align.traceback.ref[i]=='N']
        if len(nidx) == 17+19+19+12:
            r1 = align.traceback.query[:nidx[0]]
            l1 = align.traceback.query[nidx[16]+1:nidx[17]]
            bc1 = align.traceback.query[nidx[0]:nidx[16]+1]
            bc2 = align.traceback.query[nidx[17]:nidx[35]+1]
            bc3 = align.traceback.query[nidx[36]:nidx[54]+1]
            rawumi = align.traceback.query[nidx[55]:nidx[66]+1]
        else:  # 有N被clip的read设置为极端情况，后续会考虑其他方案
            r1 = l1 = 'N' * 20
            bc1, bc2, bc3, rawumi = 'N', 'N', 'N', 'N'
        editread1 = ed.eval(read1, r1)
        editlink1 = ed.eval(link1, l1)
        edit = ed.eval(align.traceback.query, align.traceback.ref)
        #tmp.append([i.query_name,  bc1, bc2, bc3, rawumi, i.is_mapped, align.score, align.traceback.comp,
        #            align.traceback.ref, align.traceback.query,
        #            align.query, frd.score, rev.score,
        #            edit, editread1, editlink1, direction])
        tmp.append([i.query_name,  bc1, bc2, bc3, rawumi, i.is_mapped, align.score,
                    align.query, edit, editread1, editlink1, ssp_align.score,
                    sum(i.query_qualities), i.query_length])
    #keep_cols = ['id', 'bc1', 'bc2', 'bc3', 'rawumi', 'is_mapped', 'score', 'comp',
    #             'ref', 'query', 'rawquery', 'frdscore', 'revscore', 'edit',
    #             'editread1', 'editlink1', 'direction']
    keep_cols = ['id', 'bc1', 'bc2', 'bc3', 'rawumi', 'is_mapped', 'score',
                 'rawquery', 'edit', 'editread1', 'editlink1', 'ssp_score',
                 'qscore', 'qlen']
    info = pd.DataFrame(tmp, columns=keep_cols)
    info.set_index('id', inplace=True)
    info['is_full_length'] =  np.logical_and(info.score>200, info.ssp_score>100)
    #sam.close()
    return info

