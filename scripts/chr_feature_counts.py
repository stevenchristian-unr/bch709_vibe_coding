#!/usr/bin/env python3
import gzip
import os
from collections import defaultdict

GFF_PATH = "data/saccharomyces_cerevisiae.gff.gz"
CHROM_SIZES_PATH = "data/chrom.sizes"
OUT_COUNTS = "results/chr_feature_counts.tsv"
OUT_DROPPED = "results/dropped_seqids.txt"


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def load_chrom_sizes(path):
    sizes = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes


def main():
    chrom_sizes = load_chrom_sizes(CHROM_SIZES_PATH)
    chrom_set = set(chrom_sizes.keys())

    n_gene = defaultdict(int)
    n_trna = defaultdict(int)
    n_snorna = defaultdict(int)
    exon_intervals = defaultdict(set)  # (start, end, strand) unique per chrom
    dropped_seqids = set()
    excluded_lines = 0

    with open_maybe_gzip(GFF_PATH) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts

            if seqid not in chrom_set:
                dropped_seqids.add(seqid)
                excluded_lines += 1
                continue

            if ftype == "gene":
                n_gene[seqid] += 1
            elif ftype == "exon":
                exon_intervals[seqid].add((start, end, strand))
            elif ftype == "tRNA":
                n_trna[seqid] += 1
            elif ftype == "snoRNA":
                n_snorna[seqid] += 1

    # Ensure output directory exists
    os.makedirs(os.path.dirname(OUT_COUNTS), exist_ok=True)

    # Write dropped seqids
    with open(OUT_DROPPED, "w") as f:
        for s in sorted(dropped_seqids):
            f.write(f"{s}\n")

    # Build output table
    rows = []
    for chrom in sorted(chrom_sizes.keys()):
        chrom_len = chrom_sizes[chrom]
        denom = chrom_len / 1e6 if chrom_len else 0
        gene_ct = n_gene.get(chrom, 0)
        exon_ct = len(exon_intervals.get(chrom, set()))
        trna_ct = n_trna.get(chrom, 0)
        snorna_ct = n_snorna.get(chrom, 0)
        gene_per_mb = round(gene_ct / denom, 4) if denom else 0.0
        exon_per_mb = round(exon_ct / denom, 4) if denom else 0.0
        trna_per_mb = round(trna_ct / denom, 4) if denom else 0.0
        snorna_per_mb = round(snorna_ct / denom, 4) if denom else 0.0
        rows.append(
            {
                "chrom": chrom,
                "chrom_length_bp": chrom_len,
                "n_gene": gene_ct,
                "n_exon_unique": exon_ct,
                "n_tRNA": trna_ct,
                "n_snoRNA": snorna_ct,
                "gene_per_Mb": gene_per_mb,
                "exon_unique_per_Mb": exon_per_mb,
                "tRNA_per_Mb": trna_per_mb,
                "snoRNA_per_Mb": snorna_per_mb,
            }
        )

    rows.sort(key=lambda r: r["gene_per_Mb"], reverse=True)

    # Write TSV with header
    header = [
        "chrom",
        "chrom_length_bp",
        "n_gene",
        "n_exon_unique",
        "n_tRNA",
        "n_snoRNA",
        "gene_per_Mb",
        "exon_unique_per_Mb",
        "tRNA_per_Mb",
        "snoRNA_per_Mb",
    ]
    with open(OUT_COUNTS, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write(
                f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
                f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}\n"
            )

    # Print summary and top 5 rows
    print(f"dropped_seqids: {len(dropped_seqids)}")
    print(f"excluded_feature_lines: {excluded_lines}")
    print("\t".join(header))
    for r in rows[:5]:
        print(
            f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
            f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}\t"
            f"{r['gene_per_Mb']}\t{r['exon_unique_per_Mb']}\t"
            f"{r['tRNA_per_Mb']}\t{r['snoRNA_per_Mb']}"
        )


if __name__ == "__main__":
    main()
