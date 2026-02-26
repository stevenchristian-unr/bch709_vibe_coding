import gzip, re
from collections import defaultdict

counts = defaultdict(lambda: defaultdict(int))

with gzip.open("data/saccharomyces_cerevisiae.gff.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        ftype = fields[2]
        if ftype in ("gene", "exon", "tRNA", "snoRNA"):
            counts[chrom][ftype] += 1

for chrom in sorted(counts):
    print(chrom, dict(counts[chrom]))
