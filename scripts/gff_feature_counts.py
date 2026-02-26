#!/usr/bin/env python3
# scripts/gff_feature_counts.py

import gzip
from collections import defaultdict

GFF_PATH = "data/saccharomyces_cerevisiae.gff.gz"

FEATURES = {"gene", "exon", "tRNA", "snoRNA"}


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def main():
    counts = defaultdict(lambda: defaultdict(int))

    with open_maybe_gzip(GFF_PATH) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype in FEATURES:
                counts[seqid][ftype] += 1

    # print header
    print("chromosome\tgene\texon\ttRNA\tsnoRNA")
    for chr_name in sorted(counts.keys()):
        print(
            f"{chr_name}\t"
            f"{counts[chr_name].get('gene', 0)}\t"
            f"{counts[chr_name].get('exon', 0)}\t"
            f"{counts[chr_name].get('tRNA', 0)}\t"
            f"{counts[chr_name].get('snoRNA', 0)}"
        )


if __name__ == "__main__":
    main()
