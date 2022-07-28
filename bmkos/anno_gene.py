import pandas as pd
import bioframe as bf
import pysam


def load_gtf(gtf):
    cols = [
        "chrom", "source", "feature", "start",
        "end", "score", "strand", "frame", "attribute",
    ]
    df = pd.read_csv(gtf, sep="\t", comment="#", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "GTF file not loading as a valid dataframe!"
        # Restrict annotations to the gene level and isolate gene name from attributes
        df = df[df["feature"] == "gene"]
        df['gene_name'] = df.attribute.str.extract(r'gene_name \"(.*?)\";')
        df['attribute'] = df.attribute.str.extract(r'gene_id \"(.*?)\";')

    return df



def bam2bed(bam):
    cols = [ "chrom", "start", "end",
            "name", "score", "strand",]
    #bam = pysam.AlignmentFile(bam, 'r')
    record = [(i.reference_name, i.reference_start,
               i.reference_end, i.query_name, i.mapq,
               i.is_reverse and '-' or '+')
               for i in bam if i.is_mapped]
    #bam.close()
    return pd.DataFrame(record, columns=cols)




def assign_gene(bed, gtf):
    df = bf.overlap(bed, gtf, how="left", suffixes=("_bed", "_gtf"),
                    return_overlap=True, return_index=True)
    
    df = df[["index_bed", "name_bed", "chrom_bed", "score_bed",
    "strand_gtf", "attribute_gtf", "overlap_start",
    "overlap_end", 'gene_name_gtf']].fillna(0)
    
    df = df.rename(
        columns={
            "name_bed": "read",
            "chrom_bed": "chrom",
            "score_bed": "score",
            "strand_gtf": "strand",
            "attribute_gtf": "gene",
            "gene_name_gtf": "gene_name",
        }
    )
    
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]
    df["status"] = "Unknown"
    df.loc[df["score"] < 30, "status"] = "Unassigned_mapq"
    df.loc[df["score"] < 30, "gene"] = "NA"
    is_ambiguous = df[df["status"] == "Unknown"].duplicated(
        subset=["index_bed", "overlap_bp"], keep=False
    )
    ambiguous_idx = is_ambiguous.index[is_ambiguous]
    df.loc[ambiguous_idx, "status"] = "Unassigned_ambiguous"
    df.loc[ambiguous_idx, "gene"] = "NA"
    df = df.drop_duplicates(subset=["index_bed", "overlap_bp", "status"])
    df.loc[df["gene"] == 0, "status"] = "Unassigned_no_features"
    df.loc[df["gene"] == 0, "gene"] = "NA"
    max_ovlp_idx = df.groupby(["index_bed"])["overlap_bp"].idxmax().sort_values().values
    df = df.loc[max_ovlp_idx, :]
    unknown_idx = [i for i in max_ovlp_idx if df.loc[i, "status"] == "Unknown"]
    df.loc[unknown_idx, "status"] = "Assigned"
    df = df[["read", "status", "score", "gene", "gene_name", "index_bed"]]
    df.reset_index(drop=True, inplace=True)
    df = df.drop(["index_bed"], axis=1)
    df.set_index('read', inplace=True)
    df.loc[df[df.gene == 'NA'].index, 'gene_name'] = 'NA'
    return df

