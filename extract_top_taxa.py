#!/usr/bin/env python3
"""
Extract top differentially abundant taxa from BIRDMAn results for KEGG pathway analysis.

This script:
1. Parses BIRDMAn output to identify significant features
2. Filters based on effect size and credible intervals
3. Exports results for downstream KEGG pathway enrichment analysis
"""

import pandas as pd
import numpy as np
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


def load_birdman_results(results_dir):
    """
    Load BIRDMAn results from exported directory.

    Parameters:
    -----------
    results_dir : str
        Path to exported BIRDMAn results

    Returns:
    --------
    pd.DataFrame
        BIRDMAn results with effect sizes and credible intervals
    """
    results_path = Path(results_dir)

    # BIRDMAn exports different files - look for the main results
    potential_files = [
        results_path / "differentials.tsv",
        results_path / "feature_table.tsv",
        results_path / "summary.tsv",
    ]

    for file_path in potential_files:
        if file_path.exists():
            print(f"Loading results from: {file_path}")
            df = pd.read_csv(file_path, sep="\t")
            return df

    # If no standard file found, try to find any TSV
    tsv_files = list(results_path.glob("*.tsv"))
    if tsv_files:
        print(f"Loading results from: {tsv_files[0]}")
        return pd.read_csv(tsv_files[0], sep="\t")

    raise FileNotFoundError(f"No BIRDMAn results found in {results_dir}")


def identify_significant_features(
    df,
    effect_size_col="log_fold_change",
    lower_ci_col="ci_lower",
    upper_ci_col="ci_upper",
    effect_size_threshold=2,
    ci_exclude_zero=True,
):
    """
    Identify statistically significant and biologically meaningful features.

    Parameters:
    -----------
    df : pd.DataFrame
        BIRDMAn results dataframe
    effect_size_col : str
        Column name for effect size
    lower_ci_col : str
        Column name for lower credible interval
    upper_ci_col : str
        Column name for upper credible interval
    effect_size_threshold : float
        Minimum absolute effect size
    ci_exclude_zero : bool
        Require credible interval to exclude zero

    Returns:
    --------
    tuple
        (significant_features_df, enriched_df, depleted_df)
    """
    # Make a copy to avoid modifying original
    results = df.copy()

    # Calculate absolute effect size
    results["abs_effect_size"] = results[effect_size_col].abs()

    # Filter by effect size
    significant = results[results["abs_effect_size"] >= effect_size_threshold].copy()

    # Filter by credible interval (if CI doesn't include 0)
    if (
        ci_exclude_zero
        and lower_ci_col in results.columns
        and upper_ci_col in results.columns
    ):
        significant = significant[
            ~((significant[lower_ci_col] < 0) & (significant[upper_ci_col] > 0))
        ]

    # Separate enriched (positive) and depleted (negative)
    enriched = significant[significant[effect_size_col] > 0].copy()
    depleted = significant[significant[effect_size_col] < 0].copy()

    # Sort by absolute effect size
    enriched = enriched.sort_values("abs_effect_size", ascending=False)
    depleted = depleted.sort_values("abs_effect_size", ascending=False)

    print(f"\nSignificant Features Summary:")
    print(f"  Total significant: {len(significant)}")
    print(f"  Enriched (positive): {len(enriched)}")
    print(f"  Depleted (negative): {len(depleted)}")

    return significant, enriched, depleted


def export_for_kegg_analysis(
    enriched_df, depleted_df, output_dir, feature_id_col="feature_id"
):
    """
    Export feature lists for KEGG pathway analysis.

    Parameters:
    -----------
    enriched_df : pd.DataFrame
        Enriched features
    depleted_df : pd.DataFrame
        Depleted features
    output_dir : str
        Output directory
    feature_id_col : str
        Column name for feature IDs
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Export enriched features
    enriched_ids = enriched_df[feature_id_col].tolist()
    with open(output_path / "enriched_features.txt", "w") as f:
        f.write("\n".join(enriched_ids))

    # Export depleted features
    depleted_ids = depleted_df[feature_id_col].tolist()
    with open(output_path / "depleted_features.txt", "w") as f:
        f.write("\n".join(depleted_ids))

    # Export all significant features
    all_significant = pd.concat([enriched_df, depleted_df])
    all_significant.to_csv(
        output_path / "significant_features.tsv", sep="\t", index=False
    )

    # Export summary
    summary = {
        "total_significant": len(all_significant),
        "enriched": len(enriched_df),
        "depleted": len(depleted_df),
        "enriched_features": enriched_ids[:20],  # Top 20
        "depleted_features": depleted_ids[:20],  # Top 20
    }

    with open(output_path / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n✓ Results exported to: {output_path}")
    print(f"  - enriched_features.txt: {len(enriched_ids)} features")
    print(f"  - depleted_features.txt: {len(depleted_ids)} features")
    print(f"  - significant_features.tsv: Complete results")


def create_visualizations(
    enriched_df,
    depleted_df,
    output_dir,
    effect_size_col="log_fold_change",
    feature_id_col="feature_id",
):
    """
    Create visualizations of differential abundance results.
    """
    output_path = Path(output_dir)

    # Combine data
    all_data = pd.concat(
        [
            enriched_df.assign(direction="Enriched"),
            depleted_df.assign(direction="Depleted"),
        ]
    )

    # 1. Forest plot of top features
    fig, ax = plt.subplots(figsize=(10, 12))

    # Get top 20 features by absolute effect size
    top_features = all_data.nlargest(20, "abs_effect_size")

    # Create forest plot
    colors = ["red" if d == "Depleted" else "blue" for d in top_features["direction"]]
    y_pos = np.arange(len(top_features))

    ax.barh(y_pos, top_features[effect_size_col], color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_features[feature_id_col])
    ax.set_xlabel("Log Fold Change (Effect Size)")
    ax.set_title("Top 20 Differentially Abundant Features")
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.5)

    # Add legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="blue", alpha=0.7, label="Enriched with SCFA"),
        Patch(facecolor="red", alpha=0.7, label="Depleted with SCFA"),
    ]
    ax.legend(handles=legend_elements)

    plt.tight_layout()
    plt.savefig(output_path / "forest_plot.png", dpi=300, bbox_inches="tight")
    print(f"  - forest_plot.png: Top features visualization")

    # 2. Distribution of effect sizes
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(
        enriched_df[effect_size_col], bins=30, alpha=0.7, label="Enriched", color="blue"
    )
    ax.hist(
        depleted_df[effect_size_col], bins=30, alpha=0.7, label="Depleted", color="red"
    )
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Effect Sizes")
    ax.legend()
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig(
        output_path / "effect_size_distribution.png", dpi=300, bbox_inches="tight"
    )
    print(f"  - effect_size_distribution.png: Effect size distribution")

    plt.close("all")


def generate_kegg_analysis_guide(output_dir):
    """
    Generate a guide for KEGG pathway analysis.
    """
    output_path = Path(output_dir)

    guide = """
# KEGG Pathway Analysis Guide

## Step 1: Functional Prediction

Use PICRUSt2 to predict functional potential from your 16S data:

```bash
# If you have representative sequences
picrust2_pipeline.py \\
  -s rep-seqs.fna \\
  -i feature-table.biom \\
  -o picrust2_output \\
  -p 32

# Or use QIIME2 plugin
qiime picrust2 full-pipeline \\
  --i-table feature-table.qza \\
  --i-seq rep-seqs.qza \\
  --output-dir picrust2_results \\
  --p-threads 32
```

## Step 2: Identify KEGG Pathways

PICRUSt2 will generate:
- `KO_metagenome_out/pred_metagenome_unstrat.tsv`: KEGG Orthology predictions
- `pathways_out/path_abun_unstrat.tsv`: Pathway abundances

## Step 3: Extract Features of Interest

Use the feature IDs from `enriched_features.txt` and `depleted_features.txt` 
to extract their predicted KEGG functions.

## Step 4: Over-Representation Analysis

### Expected Lipid/SCFA-Related Pathways:

1. **Butanoate metabolism (map00650)**
   - Butyrate production pathway
   - Key for SCFA hypothesis

2. **Propanoate metabolism (map00640)**
   - Propionate production
   - Linked to gut health

3. **Fatty acid biosynthesis (map00061)**
   - Short-chain fatty acid synthesis

4. **Fatty acid degradation (map00071)**
   - Beta-oxidation of fatty acids

5. **Lipopolysaccharide biosynthesis (map00540)**
   - Cell wall components

### Statistical Test

Use hypergeometric test or Fisher's exact test to determine if 
differentially abundant taxa are enriched in these pathways compared 
to background.

## Step 5: Visualization

Create pathway enrichment plots showing:
- Number of differentially abundant taxa per pathway
- p-values for enrichment
- Gene counts

## Tools for Pathway Analysis

### Option 1: KEGG Mapper
- Web interface: https://www.genome.jp/kegg/mapper/
- Upload KO list for pathway mapping

### Option 2: Python (GOATOOLS/scipy)
```python
from scipy.stats import hypergeom
import pandas as pd

# Calculate enrichment p-value
M = total_genes  # Total genes in database
n = genes_in_pathway  # Genes in specific pathway
N = diff_abundant_genes  # Your differentially abundant genes
k = overlap  # Overlap between your genes and pathway

pval = hypergeom.sf(k-1, M, n, N)
```

### Option 3: R (clusterProfiler)
```r
library(clusterProfiler)

enrichKEGG(gene = your_gene_list,
          organism = 'ko',
          pvalueCutoff = 0.05)
```

## Expected Results (if hypothesis is true)

✓ Enrichment of taxa in SCFA production pathways (butanoate, propanoate)
✓ Association with lipid metabolism pathways
✓ Higher representation than expected by chance (p < 0.05)

## Interpreting Results

- **Positive enrichment**: More differentially abundant taxa in pathway than expected
- **p-value < 0.05**: Statistically significant enrichment
- **Fold enrichment > 2**: Strong biological signal

## References

1. Douglas GM, et al. (2020). PICRUSt2 for prediction of metagenome functions. 
   Nature Biotechnology, 38(6), 685-688.

2. Kanehisa M, Goto S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. 
   Nucleic Acids Research, 28(1), 27-30.
"""

    with open(output_path / "KEGG_ANALYSIS_GUIDE.md", "w") as f:
        f.write(guide)

    print(f"  - KEGG_ANALYSIS_GUIDE.md: Step-by-step guide")


def main():
    parser = argparse.ArgumentParser(
        description="Extract top differentially abundant taxa for KEGG analysis"
    )
    parser.add_argument(
        "--birdman-results",
        required=True,
        help="Path to exported BIRDMAn results directory",
    )
    parser.add_argument(
        "--output-dir", default="./kegg_analysis", help="Output directory for results"
    )
    parser.add_argument(
        "--effect-size-threshold",
        type=float,
        default=2.0,
        help="Minimum absolute effect size (default: 2.0)",
    )
    parser.add_argument(
        "--exclude-zero",
        action="store_true",
        help="Require credible interval to exclude zero",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Extracting Differentially Abundant Taxa for KEGG Analysis")
    print("=" * 60)

    # Load results
    print("\n1. Loading BIRDMAn results...")
    results_df = load_birdman_results(args.birdman_results)
    print(f"   Total features: {len(results_df)}")

    # Identify significant features
    print("\n2. Identifying significant features...")
    significant, enriched, depleted = identify_significant_features(
        results_df,
        effect_size_threshold=args.effect_size_threshold,
        ci_exclude_zero=args.exclude_zero,
    )

    # Export results
    print("\n3. Exporting results...")
    export_for_kegg_analysis(enriched, depleted, args.output_dir)

    # Create visualizations
    print("\n4. Creating visualizations...")
    create_visualizations(enriched, depleted, args.output_dir)

    # Generate KEGG analysis guide
    print("\n5. Generating KEGG analysis guide...")
    generate_kegg_analysis_guide(args.output_dir)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)
    print(f"\nResults saved to: {args.output_dir}")
    print("\nNext steps:")
    print("1. Run PICRUSt2 for functional prediction")
    print("2. Follow KEGG_ANALYSIS_GUIDE.md for pathway analysis")
    print(f"3. Use enriched_features.txt and depleted_features.txt for ORA")


if __name__ == "__main__":
    main()
