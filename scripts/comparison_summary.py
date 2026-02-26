#!/usr/bin/env python3
import os

OUT_TSV = "results/comparison_summary.tsv"

HEADER = ["Aspect", "Stage 1 (Vague)", "Stage 2 (Format)", "Stage 3 (Detailed)"]
ROWS = [
    ["Chromosome scope", "Everything", "chrom.sizes only", "chrom.sizes + zero-fill"],
    ["Exon definition", "Duplicate-counted", "Unique interval", "Unique interval"],
    ["QC artifact", "None", "dropped_seqids.txt", "Count + line count + file"],
    ["Density", "None", "None", "4 per_Mb columns"],
    ["Sorting", "None", "None", "gene_per_Mb descending"],
    ["Reusability", "Low", "Medium", "High (publication-ready)"],
]


def main():
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w") as f:
        f.write("\t".join(HEADER) + "\n")
        for row in ROWS:
            f.write("\t".join(row) + "\n")

    # Print to console
    print("\t".join(HEADER))
    for row in ROWS:
        print("\t".join(row))
    print(f"Saved: {OUT_TSV}")


if __name__ == "__main__":
    main()

