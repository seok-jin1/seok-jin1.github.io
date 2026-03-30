---
layout: post
title: "GWAS Summary Statistics for Precision Medicine: From Manhattan Plot to Candidate Genes"
date: 2026-03-29
permalink: /blog/gwas-summary-statistics-precision-medicine/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - genomics
  - precision-medicine
  - python
  - tutorial
description: "A practical guide to reading GWAS summary statistics, visualizing association signals, and turning loci into biologically meaningful hypotheses."
---

Precision medicine projects often start with a deceptively simple file: a table of variants and association statistics. You may not have access to the raw genotype data, the imputation pipeline, or the original cohort. But if you have **GWAS summary statistics**, you can already do a surprising amount of useful work: inspect signal quality, identify sentinel loci, prioritize genes, compare traits, and connect disease risk to expression, splicing, and cellular state.

That is why summary statistics matter. They are usually the first bridge between **population genetics** and **mechanism**.

This tutorial focuses on the practical side. We will not run a full GWAS from raw genotype files. Instead, we will assume that you already have a summary-statistics file and want to answer the questions that matter downstream:

1. Is the file usable?
2. Where are the strongest loci?
3. Which signals are likely independent?
4. How do I go from a locus to a precision-medicine hypothesis?

By the end, you will have a compact Python workflow for loading and cleaning summary statistics, making Manhattan and QQ plots, extracting sentinel variants, and preparing the file for downstream tools such as PLINK 2 clumping, colocalization, and single-cell integration.

{% include figure.liquid loading="eager" path="assets/img/blog/gwas-summary/figure1-workflow.svg" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 1. A practical workflow from raw GWAS summary statistics to precision-medicine hypotheses. The key transition is not association alone, but moving from signal to mechanism, cellular context, and clinical interpretation.
</div>

---

## 1. What Summary Statistics Usually Contain

At minimum, a usable GWAS summary-statistics file should contain:

| Column         | Meaning                               |
| -------------- | ------------------------------------- |
| `CHR`          | chromosome                            |
| `POS`          | genomic position                      |
| `SNP` or `ID`  | variant identifier, often rsID        |
| `EA`           | effect allele                         |
| `OA`           | other allele                          |
| `BETA` or `OR` | effect size                           |
| `SE`           | standard error                        |
| `P`            | association p-value                   |
| `EAF`          | effect-allele frequency, if available |
| `N`            | sample size, if available             |

In practice, files are messy. Some studies use `A1/A2`, some use `ALT/REF`, some report `OR` instead of `BETA`, and some use hg19 coordinates while others use hg38. Before you interpret anything, you need to standardize the schema.

For public studies, the [GWAS Catalog summary statistics portal](https://www.ebi.ac.uk/gwas/summary-statistics) is a good starting point. If you run association analysis yourself, [PLINK 2](https://www.cog-genomics.org/plink/2.0/assoc) is one of the most widely used tools, and it can also reformat association outputs into GWAS-SSF for the GWAS Catalog.

---

## 2. Setup

We will use a small Python stack:

```bash
conda create -n gwas-summary python=3.11 -y
conda activate gwas-summary
pip install pandas numpy matplotlib scipy
```

If you also want LD-based clumping later, install PLINK 2 separately and make sure `plink2` is on your `PATH`.

---

## 3. Load and Standardize the File

Suppose your input file is `trait.sumstats.tsv.gz`.

```python
from pathlib import Path
import pandas as pd

sumstats_path = Path("trait.sumstats.tsv.gz")
df = pd.read_csv(sumstats_path, sep="\t")

print(df.head())
print(df.columns.tolist())
print(df.shape)
```

Now standardize common column aliases into one schema.

```python
COLUMN_MAP = {
    "chr": "CHR",
    "chrom": "CHR",
    "chromosome": "CHR",
    "bp": "POS",
    "position": "POS",
    "rsid": "SNP",
    "variant_id": "SNP",
    "id": "SNP",
    "a1": "EA",
    "alt": "EA",
    "effect_allele": "EA",
    "a2": "OA",
    "ref": "OA",
    "other_allele": "OA",
    "beta": "BETA",
    "effect": "BETA",
    "or": "OR",
    "se": "SE",
    "p": "P",
    "pval": "P",
    "p_value": "P",
    "eaf": "EAF",
    "effect_allele_frequency": "EAF",
    "n": "N",
}

df = df.rename(columns={c: COLUMN_MAP.get(c.lower(), c) for c in df.columns})

required = ["CHR", "POS", "EA", "OA", "P"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

print(df[required].head())
```

If your file reports odds ratios instead of betas, convert to log-odds.

```python
import numpy as np

if "BETA" not in df.columns and "OR" in df.columns:
    df["BETA"] = np.log(df["OR"])
```

---

## 4. Basic QC Before Interpretation

Summary-statistics QC is not glamorous, but it prevents avoidable mistakes.

### 4.1 Remove obviously invalid rows

```python
df = df.copy()

df["P"] = pd.to_numeric(df["P"], errors="coerce")
df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
df["CHR"] = df["CHR"].astype(str).str.replace("chr", "", regex=False)

df = df.dropna(subset=["CHR", "POS", "EA", "OA", "P"])
df = df[(df["P"] > 0) & (df["P"] <= 1)]
df = df[df["EA"].str.len() == 1]
df = df[df["OA"].str.len() == 1]

print(df.shape)
```

### 4.2 Remove duplicated variants

```python
key_cols = ["CHR", "POS", "EA", "OA"]
df = df.drop_duplicates(subset=key_cols)
```

### 4.3 Flag strand-ambiguous SNPs

`A/T` and `C/G` SNPs are strand-ambiguous. They are not always wrong, but they become risky when you merge datasets or compare studies across builds and allele conventions.

```python
AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
df["is_ambiguous"] = list(zip(df["EA"], df["OA"]))
df["is_ambiguous"] = df["is_ambiguous"].isin(AMBIGUOUS)

print(df["is_ambiguous"].mean())
```

If you are only plotting one file, you can keep them. If you are harmonizing across traits, it is often safer to exclude them unless allele frequencies clearly rescue orientation.

### 4.4 Compute a Z-score if possible

```python
if {"BETA", "SE"}.issubset(df.columns):
    df["Z"] = df["BETA"] / df["SE"]
```

---

## 5. Manhattan Plot

The Manhattan plot answers a simple question: where in the genome are the strongest association signals?

```python
import matplotlib.pyplot as plt
import numpy as np

plot_df = df.copy()
plot_df["CHR_NUM"] = pd.to_numeric(plot_df["CHR"], errors="coerce")
plot_df = plot_df.dropna(subset=["CHR_NUM"]).sort_values(["CHR_NUM", "POS"])
plot_df["minus_log10_p"] = -np.log10(plot_df["P"])

chrom_offsets = {}
offset = 0
xticks = []
xticklabels = []

for chrom, subdf in plot_df.groupby("CHR_NUM", sort=True):
    chrom_offsets[chrom] = offset
    midpoint = offset + (subdf["POS"].max() - subdf["POS"].min()) / 2
    xticks.append(midpoint)
    xticklabels.append(str(int(chrom)))
    offset += subdf["POS"].max()

plot_df["genome_pos"] = plot_df.apply(
    lambda row: row["POS"] + chrom_offsets[row["CHR_NUM"]],
    axis=1,
)

fig, ax = plt.subplots(figsize=(14, 5))
colors = ["#4c78a8", "#f58518"]

for i, (chrom, subdf) in enumerate(plot_df.groupby("CHR_NUM", sort=True)):
    ax.scatter(
        subdf["genome_pos"],
        subdf["minus_log10_p"],
        s=6,
        color=colors[i % 2],
        alpha=0.8,
        rasterized=True,
    )

ax.axhline(-np.log10(5e-8), color="crimson", linestyle="--", linewidth=1)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
ax.set_xlabel("Chromosome")
ax.set_ylabel("-log10(P)")
ax.set_title("Manhattan plot")
plt.tight_layout()
plt.show()
```

The dashed line at `5e-8` is the conventional genome-wide significance threshold for common-variant GWAS. It is a useful convention, not a law of nature. In sequencing studies, rare-variant analyses, or ancestry-specific studies, the appropriate threshold may differ.

{% include figure.liquid loading="eager" path="assets/img/blog/gwas-summary/figure3-manhattan-plot.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. Manhattan plot of GWAS association results. Each dot represents a tested variant; the y-axis shows the strength of association. The red dashed line marks genome-wide significance (5×10⁻⁸). Peaks indicate genomic loci harboring disease-associated variants.
</div>

---

## 6. QQ Plot

The QQ plot helps you decide whether the observed test statistics deviate from the null more than expected. If the entire curve inflates upward, you may be looking at population structure, cryptic relatedness, poor QC, or batch effects. If only the far tail deviates, that is more consistent with genuine signals.

```python
qq_df = df[["P"]].dropna().copy()
qq_df = qq_df[(qq_df["P"] > 0) & (qq_df["P"] <= 1)].sort_values("P")

observed = -np.log10(qq_df["P"].values)
expected = -np.log10(np.arange(1, len(observed) + 1) / (len(observed) + 1))

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(expected, observed, s=8, alpha=0.7, color="#4c78a8")
max_val = max(expected.max(), observed.max())
ax.plot([0, max_val], [0, max_val], linestyle="--", color="gray")
ax.set_xlabel("Expected -log10(P)")
ax.set_ylabel("Observed -log10(P)")
ax.set_title("QQ plot")
plt.tight_layout()
plt.show()
```

If you have access to lambda GC or LDSC intercept estimates from the original study, interpret the QQ plot together with those numbers. The plot alone cannot distinguish polygenicity from confounding.

{% include figure.liquid loading="eager" path="assets/img/blog/gwas-summary/figure4-qq-plot.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. QQ plot comparing observed versus expected p-value distributions. The diagonal line represents the null expectation. Deviation in the upper tail indicates genuine association signals. The genomic inflation factor (λ_GC) quantifies systematic inflation.
</div>

---

## 7. Extract Sentinel Variants

Before more sophisticated fine-mapping, it is useful to extract a set of lead variants. Here is a simple distance-based approach that greedily picks the lowest-p variant and suppresses neighbors within `+/-500 kb`.

```python
window_bp = 500_000
lead_rows = []

for chrom, subdf in df.sort_values("P").groupby("CHR", sort=False):
    used_positions = []
    for _, row in subdf.sort_values("P").iterrows():
        if all(abs(row["POS"] - pos) > window_bp for pos in used_positions):
            lead_rows.append(row)
            used_positions.append(row["POS"])

lead_df = pd.DataFrame(lead_rows).sort_values(["CHR", "POS"])
print(lead_df[["CHR", "POS", "SNP", "P"]].head(20))
```

This is useful for exploration, but it is not LD-aware. In real analyses you usually want **LD-based clumping** or fine-mapping.

With PLINK 2, clumping is the standard next step:

```bash
plink2 \
  --pfile reference_panel \
  --clump trait.sumstats.tsv \
  --clump-p1 5e-8 \
  --clump-p2 1e-4 \
  --clump-r2 0.1 \
  --clump-kb 500 \
  --out trait
```

According to the [PLINK 2 report postprocessing documentation](https://www.cog-genomics.org/plink/2.0/postproc), `--clump` groups association results into LD-based clumps so that nearby correlated hits are not over-counted as independent biology.

{% include figure.liquid loading="eager" path="assets/img/blog/gwas-summary/figure5-locus-zoom.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Locus zoom plot of the top GWAS signal. Variants are colored by linkage disequilibrium (r²) with the lead SNP (purple diamond). This view reveals the extent of the association signal and helps distinguish independent signals from correlated hits.
</div>

---

## 8. From Locus to Candidate Genes

This is where many analyses go wrong. The nearest gene is not always the causal gene. A risk variant can act through:

- a coding change in the nearest gene
- a splice-altering effect in a nearby gene
- a distal enhancer that regulates a gene hundreds of kilobases away
- a cell-type-specific chromatin loop

That is why a precision-medicine workflow should move in layers:

1. **Association layer**: sentinel SNP or credible set
2. **Regulatory layer**: eQTL, sQTL, chromatin accessibility, enhancer links
3. **Cellular layer**: which cell types and cell states show enrichment
4. **Clinical layer**: is the gene druggable, prognostic, or predictive of treatment response

If your blog or research direction already includes AlphaGenome, splicing models, single-cell atlases, and cellular dynamics, this is where GWAS becomes useful. GWAS does not replace those methods. It gives you a ranked list of loci that those methods can interpret mechanistically.

---

## 9. Why This Matters for Precision Medicine

Summary statistics are not just for population geneticists. They become clinically useful when they help you answer questions like:

- Which loci differentiate responders from non-responders?
- Which immune cell states carry disease-associated regulatory burden?
- Which variants converge on a druggable pathway?
- Which loci should be prioritized for functional validation?

For example, a disease GWAS hit that colocalizes with an eQTL in exhausted CD8 T cells is immediately more actionable than a purely positional hit with no expression evidence. Likewise, a locus that points to a splicing mechanism can connect directly to variant interpretation models such as SpliceAI and to therapeutic ideas such as antisense oligonucleotides.

In other words, **GWAS summary statistics become precision medicine when they are integrated with mechanism and context**.

{% include figure.liquid loading="eager" path="assets/img/blog/gwas-summary/figure2-integration-layers.svg" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 2. A layered interpretation model for GWAS. The most common mistake is stopping at the nearest gene; the useful path is locus to regulatory mechanism to cell-state context to clinical actionability.
</div>

---

## 10. Common Mistakes

### Treating the nearest gene as the causal gene

This is the most common over-interpretation in beginner analyses.

### Ignoring ancestry

LD structure, allele frequency, and transferability differ across ancestries. A locus that is easy to fine-map in one ancestry may be diffuse in another.

### Forgetting allele harmonization

Cross-study comparisons break quickly if effect alleles are not aligned.

### Over-reading p-values without effect sizes

Large cohorts can make tiny effects look overwhelming. Always inspect effect size, standard error, and biological plausibility together.

### Confusing locus discovery with mechanism

A Manhattan peak is the start of a story, not the end.

---

## Key Takeaways

- GWAS summary statistics are often the fastest entry point into genetic risk interpretation.
- Before any biological interpretation, standardize schema and perform basic QC.
- Manhattan and QQ plots tell you whether the file behaves like a real GWAS signal or a technical artifact.
- Distance-based lead SNP extraction is fine for exploration, but LD-aware clumping or fine-mapping is better for serious analysis.
- The real value for precision medicine comes from integrating GWAS with eQTL, splicing, chromatin, and single-cell context.

---

## References and Resources

- GWAS Catalog summary statistics portal: https://www.ebi.ac.uk/gwas/summary-statistics
- GWAS Catalog training page on summary statistics: https://www.ebi.ac.uk/training/online/courses/gwas-catalogue-exploring-snp-trait-associations/summary-statistics/
- PLINK 2 association analysis documentation: https://www.cog-genomics.org/plink/2.0/assoc
- PLINK 2 report postprocessing and clumping: https://www.cog-genomics.org/plink/2.0/postproc
- PLINK 2 output and file-format index: https://www.cog-genomics.org/plink/2.0/index
