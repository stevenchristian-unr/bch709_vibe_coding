#!/usr/bin/env python3
# scripts/feature_density_by_chr.py

import argparse
import gzip
import sys
from collections import defaultdict
from math import isnan
import csv


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_chrom_sizes(path):
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


def pearson_corr(x, y):
    # basic Pearson correlation
    n = len(x)
    if n == 0:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = sum((a - mx) ** 2 for a in x)
    deny = sum((b - my) ** 2 for b in y)
    if denx == 0 or deny == 0:
        return float("nan")
    return num / (denx**0.5 * deny**0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--gff", required=True, help="Path to saccharomyces_cerevisiae.gff.gz"
    )
    ap.add_argument("--sizes", required=True, help="Path to chrom.sizes")
    ap.add_argument("--out", default="feature_counts_by_chr.csv", help="Output CSV")
    args = ap.parse_args()

    sizes = read_chrom_sizes(args.sizes)

    # counts[chr][feature] = count
    counts = defaultdict(lambda: defaultdict(int))
    # track max end seen per chr for integrity check
    max_end = defaultdict(int)

    feature_types = {"gene", "exon", "tRNA", "snoRNA"}

    with open_maybe_gzip(args.gff) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype in feature_types:
                counts[seqid][ftype] += 1
            try:
                end_i = int(end)
                if end_i > max_end[seqid]:
                    max_end[seqid] = end_i
            except ValueError:
                pass

    # build output table
    rows = []
    for chr_name in sorted(sizes.keys()):
        chr_size = sizes.get(chr_name)
        gene_ct = counts[chr_name].get("gene", 0)
        exon_ct = counts[chr_name].get("exon", 0)
        trna_ct = counts[chr_name].get("tRNA", 0)
        snorna_ct = counts[chr_name].get("snoRNA", 0)

        # densities per Mb
        mb = chr_size / 1_000_000 if chr_size else 0
        rows.append(
            {
                "chromosome": chr_name,
                "size_bp": chr_size,
                "gene_count": gene_ct,
                "exon_count": exon_ct,
                "tRNA_count": trna_ct,
                "snoRNA_count": snorna_ct,
                "gene_per_mb": gene_ct / mb if mb else 0.0,
                "exon_per_mb": exon_ct / mb if mb else 0.0,
                "tRNA_per_mb": trna_ct / mb if mb else 0.0,
                "snoRNA_per_mb": snorna_ct / mb if mb else 0.0,
                "max_end_in_gff": max_end.get(chr_name, 0),
                "size_ok": (max_end.get(chr_name, 0) <= chr_size)
                if chr_size
                else False,
            }
        )

    # correlation between chromosome size and feature counts
    sizes_list = [r["size_bp"] for r in rows]
    gene_list = [r["gene_count"] for r in rows]
    exon_list = [r["exon_count"] for r in rows]
    trna_list = [r["tRNA_count"] for r in rows]
    snorna_list = [r["snoRNA_count"] for r in rows]

    corr_gene = pearson_corr(sizes_list, gene_list)
    corr_exon = pearson_corr(sizes_list, exon_list)
    corr_trna = pearson_corr(sizes_list, trna_list)
    corr_snorna = pearson_corr(sizes_list, snorna_list)

    with open(args.out, "w", newline="") as out_f:
        w = csv.DictWriter(out_f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # summary to stdout
    print("Wrote:", args.out)
    print("Correlation (size vs feature count):")
    print(f"  gene:  {corr_gene:.4f}")
    print(f"  exon:  {corr_exon:.4f}")
    print(f"  tRNA:  {corr_trna:.4f}")
    print(f"  snoRNA:{corr_snorna:.4f}")
    print(
        "Integrity check: any size_ok == False indicates GFF end exceeds chrom.sizes."
    )


if __name__ == "__main__":
    main()
