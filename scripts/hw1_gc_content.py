#!/usr/bin/env python3
"""
PROMPT USED
Write a Python script that reads the gzipped yeast mRNA FASTA file at data/mrna.fa.gz, parses each FASTA record, and computes for each sequence: accession (token after '>' up to first space), length, and GC content (fraction of G or C bases, case-insensitive). Output results to results/mrna_metrics.tsv with columns accession, length, gc_content and gc_content formatted to 4 decimals, sorted by gc_content descending. Then create a histogram of GC content from 0 to 1 with a density curve overlay. The plot must be 1600x900 px at 200 dpi (figsize 8x4.5 inches), include mean and median vertical dashed lines, and a caption showing n, mean, median, and sd. Save to results/gc_content_distribution.png.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = Path("data/mrna.fa.gz")
RESULTS_DIR = Path("results")
TSV_PATH = RESULTS_DIR / "mrna_metrics.tsv"
PLOT_PATH = RESULTS_DIR / "gc_content_distribution.png"


def parse_fasta_gz(path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (accession, sequence) from a gzipped FASTA file."""
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        accession = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    yield accession, "".join(seq_chunks)
                header = line[1:]
                accession = header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if accession is not None:
            yield accession, "".join(seq_chunks)


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    seq_upper = seq.upper()
    gc = seq_upper.count("G") + seq_upper.count("C")
    return gc / len(seq_upper)


def gaussian_kde_1d(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Simple Gaussian KDE using Silverman's rule of thumb for bandwidth."""
    n = values.size
    if n == 0:
        return np.zeros_like(grid)
    std = values.std(ddof=1) if n > 1 else 0.0
    if std == 0.0:
        # All values identical; return a narrow spike
        density = np.zeros_like(grid)
        idx = np.argmin(np.abs(grid - values[0]))
        density[idx] = 1.0
        return density
    bandwidth = 1.06 * std * (n ** (-1 / 5))
    if bandwidth == 0.0:
        bandwidth = 1e-6
    # Vectorized KDE
    diff = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * diff ** 2).sum(axis=1) / (n * bandwidth * math.sqrt(2 * math.pi))
    return density


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing input file: {DATA_PATH}. Download it with: curl -L -o data/mrna.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    gc_values = []
    for acc, seq in parse_fasta_gz(DATA_PATH):
        length = len(seq)
        gc = gc_content(seq)
        rows.append((acc, length, gc))
        gc_values.append(gc)

    rows.sort(key=lambda r: r[2], reverse=True)

    with TSV_PATH.open("w", encoding="utf-8") as out:
        out.write("accession\tlength\tgc_content\n")
        for acc, length, gc in rows:
            out.write(f"{acc}\t{length}\t{gc:.4f}\n")

    gc_array = np.array(gc_values, dtype=float)
    n = gc_array.size
    mean = gc_array.mean() if n else float("nan")
    median = np.median(gc_array) if n else float("nan")
    sd = gc_array.std(ddof=1) if n > 1 else float("nan")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)
    bins = np.linspace(0, 1, 51)
    ax.hist(gc_array, bins=bins, density=True, color="#5B8FF9", alpha=0.7, edgecolor="#1F3A93")

    grid = np.linspace(0, 1, 400)
    density = gaussian_kde_1d(gc_array, grid)
    ax.plot(grid, density, color="#E86452", linewidth=2)

    ax.axvline(mean, color="#333333", linestyle="--", linewidth=1.5, label="Mean")
    ax.axvline(median, color="#333333", linestyle="--", linewidth=1.5, label="Median")

    ax.set_xlim(0, 1)
    ax.set_xlabel("GC content")
    ax.set_ylabel("Density")
    ax.set_title("Yeast mRNA GC Content Distribution")

    caption = f"n={n}, mean={mean:.4f}, median={median:.4f}, sd={sd:.4f}"
    fig.text(0.5, 0.01, caption, ha="center", fontsize=9)

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(PLOT_PATH)
    plt.close(fig)


if __name__ == "__main__":
    main()

# Interpretation (Homework 1):
# 1) Yeast mRNAs may show multiple GC-content subpopulations due to gene classes with different codon usage
#    biases, which reflect expression level and translational efficiency (highly expressed genes often have
#    stronger codon bias and can skew GC content).
# 2) Local genomic context and evolutionary constraints (e.g., GC-biased gene conversion or selection on
#    mRNA secondary structure/stability) can shift GC content for subsets of transcripts, producing
#    distinct peaks or shoulders in the distribution.
