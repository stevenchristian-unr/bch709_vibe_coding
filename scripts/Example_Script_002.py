from pathlib import Path
import gzip
from collections import defaultdict

DATA_GFF = Path("data/saccharomyces_cerevisiae.gff.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV = Path("results/chr_feature_counts.tsv")
OUT_DROP = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# Load chrom sizes
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

n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_tRNA = defaultdict(int)
n_snoRNA = defaultdict(int)
dropped_seqids = set()

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
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype == "tRNA":
            n_tRNA[seqid] += 1
        elif ftype == "snoRNA":
            n_snoRNA[seqid] += 1

rows = []
for chrom in chroms:
    rows.append(
        {
            "chrom": chrom,
            "chrom_length_bp": int(chrom_len[chrom]),
            "n_gene": n_gene.get(chrom, 0),
            "n_exon_unique": len(exon_intervals.get(chrom, set())),
            "n_tRNA": n_tRNA.get(chrom, 0),
            "n_snoRNA": n_snoRNA.get(chrom, 0),
        }
    )

header = [
    "chrom",
    "chrom_length_bp",
    "n_gene",
    "n_exon_unique",
    "n_tRNA",
    "n_snoRNA",
]
with open(OUT_TSV, "w") as f:
    f.write("\t".join(header) + "\n")
    for r in rows:
        f.write(
            f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
            f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}\n"
        )
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print("Saved:", OUT_TSV)
print("Dropped seqids:", len(dropped_seqids))
print("\t".join(header))
for r in rows[:5]:
    print(
        f"{r['chrom']}\t{r['chrom_length_bp']}\t{r['n_gene']}\t"
        f"{r['n_exon_unique']}\t{r['n_tRNA']}\t{r['n_snoRNA']}"
    )
