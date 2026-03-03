# bch709_vibe_coding

Yeast genomics data processing and analysis scripts for BCH 709.

**Scope**
- Parse `GFF3` annotations and `chrom.sizes`
- Count features by chromosome (genes, exons, tRNAs, snoRNAs)
- Compute per‑Mb densities and export tabular results
- Summarize stress‑response expression variability (Gasch 2000 dataset)

**Installation**
1. Create and activate the `bch709_vibe_coding` conda environment.
2. Ensure the inputs listed below are present in `data/`.
3. If you want the full class context, see the BCH709 Vibe Coding page:
   `https://bch709.plantgenomicslab.org/vibe_coding/index.html` citeturn0view0

**How To Use**
1. Run the stage‑2 feature counts:
   ```sh
   python3 scripts/chr_feature_counts.py
   ```
2. (Optional) Run the comparison summary table:
   ```sh
   python3 scripts/comparison_summary.py
   ```
3. (Optional, R) Run the stress‑response variability script:
   ```sh
   Rscript scripts/variable_stress_genes.R
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


**Homework 1 — Yeast mRNA GC Content**

Prompt (Homework 1)
Write a script in Python to:
- Take the yeast mRNA FASTA (`mrna.fa.gz`) and yeast genome FASTA reference (`sacCer3.fa.gz`)
- Compute per-transcript metrics and output `results/mrna_metrics.tsv` with:
  - `accession`, `length`, `gc_content` (4 decimals), sorted by `gc_content` descending
- Generate `results/gc_content_distribution.png`:
  - Histogram of GC content (0–1) with a density curve overlay
  - Mean and median as vertical dashed lines
  - Caption showing `n`, `mean`, `median`, `sd`
  - 1600×900 px, dpi 200

Grading Criteria:
- Prompt quality (40%): input format, accession parsing, gzip handling, output specs (filenames, columns, decimals, graph size)
- Code correctness (40%): correct gzip parse, accurate GC computation, both files generated
- Result interpretation (20%): explain 2 possible biological reasons for the GC content distribution pattern in yeast mRNA

Interpretation (GC content distribution)
Two plausible biological reasons yeast mRNA GC content shows a characteristic distribution:

Gene expression and translation optimization (codon usage bias):
   Highly expressed genes often evolve codon usage patterns that match abundant tRNAs for efficient translation. Because codons differ in GC composition, selection on codon usage can shift GC content for subsets of transcripts, shaping the overall distribution.

Regional/functional constraints on sequence composition:
   Yeast genes vary in features like UTR length, regulatory motifs, and amino acid composition constraints (which indirectly constrain codons). In addition, mutational processes and DNA repair/replication biases can differ across genomic regions, and transcript sets sampled from those regions can inherit distinct GC tendencies—contributing to multimodality or skew in transcript GC content.
