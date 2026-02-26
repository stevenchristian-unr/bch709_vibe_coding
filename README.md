# bch709_vibe_coding

Yeast genomics data processing and analysis scripts for BCH 709.

**Scope**
- Parse `GFF3` annotations and `chrom.sizes`
- Count features by chromosome (genes, exons, tRNAs, snoRNAs)
- Compute per‑Mb densities and export tabular results
- Summarize stress‑response expression variability (Gasch 2000 dataset)

**Quick Start**
1. Activate your conda environment for this repo.
2. Run the stage‑2 feature counts:
   ```sh
   python3 scripts/chr_feature_counts.py
   ```

Outputs will appear in `results/`.

**Inputs**
- `data/saccharomyces_cerevisiae.gff.gz`: GFF3 annotations
- `data/chrom.sizes`: chromosome lengths (TSV: `chrom`, `length_bp`)
- `data/gasch2000.txt`: expression matrix used by the stress‑variability script

**Outputs**
- `results/chr_feature_counts.tsv`
  - Columns: `chrom`, `chrom_length_bp`, `n_gene`, `n_exon_unique`, `n_tRNA`, `n_snoRNA`,
    `gene_per_Mb`, `exon_unique_per_Mb`, `tRNA_per_Mb`, `snoRNA_per_Mb`
  - Sorted by `gene_per_Mb` descending
- `results/dropped_seqids.txt`
  - GFF seqids not found in `chrom.sizes` (alphabetically sorted)

**Scripts**
- `scripts/chr_feature_counts.py`
  - Stage 2 implementation used for the final counts and densities
  - Drops seqids not in `chrom.sizes` and reports excluded feature lines
- `scripts/Example_Script_002.py`
  - Feature counts (no pandas dependency)
- `scripts/Example_Script_003.py`
  - Feature counts + per‑Mb densities (no pandas dependency)
- `scripts/variable_stress_genes.R`
  - Computes most variable genes across stress conditions from `data/gasch2000.txt`

**Manual / Reference**
- Only chromosomes present in `data/chrom.sizes` are included in counts and densities.
- Exon counts are de‑duplicated by unique `(start, end, strand)` to avoid isoform overcounting.
- Counts for missing features on a chromosome are reported as `0`.
- Density is `count / (chrom_length_bp / 1e6)` and rounded to 4 decimals.
- GFF seqids not found in `chrom.sizes` are collected in `results/dropped_seqids.txt`.

**Notes**
- If a script fails with missing packages, use the no‑pandas versions in `scripts/Example_Script_002.py`
  or `scripts/Example_Script_003.py`.
