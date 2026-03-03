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


# Homework 2 — Clustering Top 200 Yeast Stress-Response Genes  
**Course:** BCH709 Bioinformatics  

---

## Project Overview

This assignment investigates whether the **top 200 most variable yeast stress-response genes** form **distinct expression clusters**, and what biological patterns characterize those clusters.

Using hierarchical clustering on normalized expression data, we determine whether stress-responsive genes group into reproducible expression modules.

---

## Input Data

**File:**
```
results/yeast_stress_cv_top200.tsv
```

**Format assumptions:**
- 200 genes (rows)
- A `gene_id` column
- Remaining columns contain **log2 expression ratios**
- The first 30 stress-condition columns are used for clustering
- Columns are kept in their **original order**
- Any metadata columns must be excluded from calculations

---

## Analysis Workflow

The R script performs the following steps:

### 1. Column Selection
- Extract the first **30 stress-condition columns**
- Do **not** cluster columns in the heatmap

---

### 2. Row-wise Z-score Normalization

For each gene:

\[
z = \frac{x - \text{row\_mean}}{\text{row\_sd}}
\]

Where:
- `row_mean` = mean expression across the 30 conditions
- `row_sd` = standard deviation across the 30 conditions
- If `row_sd = 0`, the row is set to 0 (documented in code)

---

### 3. Hierarchical Clustering

- **Distance metric:** Euclidean  
- **Linkage method:** ward.D2  
- Cluster **genes only (rows)**  
- Columns remain in original order  

---

### 4. Cluster Assignment

- Cut dendrogram at **k = 4**
- Assign cluster labels **1–4** using `cutree()`

---

## Expected Outputs

All outputs are written to the `results/` directory.

---

### 1️⃣ Clustered Heatmap

**File:**
```
results/cv_top200_cluster_heatmap.pdf
```

**Specifications:**

| Property | Value |
|-----------|--------|
| Data | Log2 ratios → row-wise Z-score |
| Rows | Hierarchical clustering (Euclidean, ward.D2) |
| Columns | First 30 stress conditions, original order |
| Annotation | k = 4 cluster color bar |
| Size | 8 × 12 inches |
| Background | White |
| Device | `pdf(width = 8, height = 12)` |

---

### 2️⃣ Cluster Assignment Table

**File:**
```
results/cluster_assignment.tsv
```

| Column | Description |
|----------|-------------|
| gene_id | Yeast gene identifier (e.g., YAL001C) |
| cluster | Cluster label (1–4) |

**Sorting:**
- Cluster ascending  
- gene_id ascending (tie-breaker)

---

## How to Run

1. Ensure the input file exists:
   ```
   results/yeast_stress_cv_top200.tsv
   ```

2. Run the script:
   ```bash
   Rscript hw2_cluster_top200.R
   ```

3. Confirm outputs:
   ```
   results/cv_top200_cluster_heatmap.pdf
   results/cluster_assignment.tsv
   ```

---

## Interpretation (Required for Submission)

Provide **exactly 4 sentences**, one per cluster, describing the biological meaning of each cluster based on observed expression patterns.

Example structure:

- **Cluster 1:** Genes broadly induced across multiple stress conditions, consistent with general environmental stress response (ESR) activation.  
- **Cluster 2:** Genes consistently repressed under stress, likely associated with growth-related or ribosomal functions.  
- **Cluster 3:** Genes specifically induced under heat or oxidative stress conditions, suggesting pathway-specific activation.  
- **Cluster 4:** Genes with mixed or transient patterns, possibly representing regulatory or condition-specific processes.  

If GO enrichment is not performed, clearly state that interpretation is based on expression patterns alone.

---

## Grading Criteria

| Criterion | Weight | Description |
|------------|---------|------------|
| Prompt Quality | 40% | Correct Z-score formula, clustering method (Euclidean + ward.D2), k=4 cutree, column handling, output specifications |
| Code Correctness | 40% | Accurate normalization, clustering, cutree, and both files generated |
| Result Interpretation | 20% | Clear, biologically reasonable description of all four clusters |

---

## Project Structure

```
.
├── hw2_cluster_top200.R
├── README.md
└── results/
    ├── yeast_stress_cv_top200.tsv
    ├── cv_top200_cluster_heatmap.pdf
    └── cluster_assignment.tsv
```

---

## Notes

- Ensure only numeric expression columns are used for calculations.
- Use deterministic clustering (no randomness).
- Explicitly control PDF device size.
- Verify normalization was **row-wise**, not column-wise.
- Document any assumptions in comments inside your R script.

---
