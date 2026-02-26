from pathlib import Path
import gzip
from collections import defaultdict

DATA_GFF = Path("data/saccharomyces_cerevisiae.gff.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV = Path("results/chr_feature_counts.tsv")
OUT_DROP = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# 1) Chrom sizes
chroms = []
chrom_len = {}
with open(DATA_SIZES, "r") as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        chroms.append(parts[0])
        chrom_len[parts[0]] = int(parts[1])
chrom_set = set(chroms)

# 2) Counters
n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_tRNA = defaultdict(int)
n_snoRNA = defaultdict(int)
dropped_seqids = set()
dropped_lines = 0

with gzip.open(DATA_GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = fields

        if seqid not in chrom_set:
            dropped_seqids.add(seqid)
            dropped_lines += 1
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype == "tRNA":
            n_tRNA[seqid] += 1
        elif ftype == "snoRNA":
            n_snoRNA[seqid] += 1

# 3) Build result table (include all chroms from chrom.sizes; fill 0 where no features)
rows = []
for chrom in chroms:
    L = float(chrom_len[chrom])
    g = n_gene.get(chrom, 0)
    ex = len(exon_intervals.get(chrom, set()))
    tr = n_tRNA.get(chrom, 0)
    sn = n_snoRNA.get(chrom, 0)
    Mb = L / 1e6 if L > 0 else 1

    rows.append(
        {
            "chrom": chrom,
            "chrom_length_bp": int(L),
            "n_gene": g,
            "n_exon_unique": ex,
            "n_tRNA": tr,
            "n_snoRNA": sn,
            "gene_per_Mb": round(g / Mb, 4),
            "exon_unique_per_Mb": round(ex / Mb, 4),
            "tRNA_per_Mb": round(tr / Mb, 4),
            "snoRNA_per_Mb": round(sn / Mb, 4),
        }
    )

rows.sort(key=lambda r: r["gene_per_Mb"], reverse=True)
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
with open(OUT_TSV, "w") as f:
    f.write("\t".join(header) + "\n")
    for r in rows:
        f.write(
            f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
            f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}\t"
            f"{r['gene_per_Mb']}\t{r['exon_unique_per_Mb']}\t"
            f"{r['tRNA_per_Mb']}\t{r['snoRNA_per_Mb']}\n"
        )
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print(f"Saved: {OUT_TSV}")
print(f"Saved: {OUT_DROP}")
print(f"Dropped seqids: {len(dropped_seqids)}")
print(f"Dropped feature lines: {dropped_lines}")
print("\t".join(header))
for r in rows[:5]:
    print(
        f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
        f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}\t"
        f"{r['gene_per_Mb']}\t{r['exon_unique_per_Mb']}\t"
        f"{r['tRNA_per_Mb']}\t{r['snoRNA_per_Mb']}"
    )
