from itertools import combinations
import editdistance as ed
import parasail

#UR: 原始UMI序列
#UB: 校正UMI序列
#R1: 原始barcode 1序列
#R2: 原始barcode 2序列
#R3: 原始barcode 3序列
#C1: 校正barcode 1序列
#C2: 校正barcode 2序列
#C3: 校正barcode 3序列
#GN: 基因
#GT: 转录本
#**N 代表木有，木有，木有呀**

# 存在umi_corr为'TTTTTTTTTTTT', 'GGGGGGGGGGGG'的情况，原因大概是UMI序列没有接上去，
# 暂时不做特别的处理，这可能会导致UMI的错误合并，而导致低估基因的表达量，
# 在后续的测试中会持续优化处理方法



#umis是包含了所有原始UMI系列的list或tuple等；
#counts是{UMI:number}的字典，
#threshold是编辑距离阈值
#函数的作用是遍历任意两个UMI序列，计算编辑距离，小于threshold的，并且
#个数差别大，个数少的被合并到个数多的中，
#或者说加入一条有向边，个数多的UMI为起点
# 当两个UMI的个数都是1时，会产生双向边
# 返回字典：邻接表，键是起点，值是终点
def get_adj_list_directional(umis, counts, threshold=3):
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if ed.eval(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


#umis是UMI序列的列表或tuple
#graph是有向图的邻接表，以字典表示
#counts是字典，每个UMI的个数
# 返回components列表，每一个元素是一个集合，集合包含了与一个node连接的所有node
def get_connected_components_adjacency(umis, graph, counts):
    found = set()
    components = list()
    #每个节点按照count数由高到低排序
    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            # component = self.search(node, graph)
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)
    return components


# node是一个节点
# adj_list是邻接表
# 返回node节点能连接的所有节点的集合
def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))
    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


# clusters: 是components
# adj_list: 是邻接表
# counts: UMI个数
# 返回值groups: 是一个列表，每个元素也是列表
def group_directional(clusters, adj_list, counts):
    observed = set()
    groups = []
    for cluster in clusters:
        if len(cluster) == 1:
            groups.append(list(cluster))
            observed.update(cluster)
        else:
            cluster = sorted(cluster, key=lambda x: counts[x], reverse=True)
            # need to remove any node which has already been observed
            temp_cluster = []
            for node in cluster:
                if node not in observed:
                    temp_cluster.append(node)
                    observed.add(node)
            groups.append(temp_cluster)

    return groups


def cluster(counts_dict, threshold=2):
    adj_list = get_adj_list_directional(counts_dict.keys(), counts_dict, threshold)
    clusters = get_connected_components_adjacency(
        counts_dict.keys(), adj_list, counts_dict
    )
    final_umis = [list(x) for x in group_directional(clusters, adj_list, counts_dict)]
    return final_umis


# cluster_list: 也就是groups，每一个group以第一个UMI序列作为校正后的UMI序列
# 返回值my_map: 字典, 键是UMI原始序列，值是校正后的UMI序列
def create_map_to_correct_umi(cluster_list):
    my_map = {y: x[0] for x in cluster_list for y in x}
    return my_map


def correct_umis(umis):
    counts_dict = dict(umis.value_counts())
    umi_map = create_map_to_correct_umi(cluster(counts_dict))
    return umis.replace(umi_map)


def get_umi(aread, matrix):
    if aread.is_keep:
        seq, bc = aread['rawquery'], aread['bc3_corr'] + 'NNNNNNNNNNNNTTTTTTTTTTTT'
        align = parasail.sw_trace(seq, bc, 4, 2, matrix)
        nidx = [i for i in range(len(align.traceback.ref)) if align.traceback.ref[i]=='N']
        umi = align.traceback.query[nidx[0]:nidx[-1]+1]
    else:
        umi = aread.rawumi
    umi = umi.replace('-', '')
    return umi

