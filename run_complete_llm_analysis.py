#!/usr/bin/env python3
"""
Complete BIRDMAn workflow for Large Larry Modeler (LLM) data.

This script runs the entire analysis pipeline from data preparation to BIRDMAn results.
Input: Your actual data files (amplicon-matrix_txt.xls, LLM_Biomarkers_Long.csv, etc.)
Output: BIRDMAn results showing taxa associated with SCFA levels.
"""

import subprocess
import sys
from pathlib import Path
import argparse
import json


def run_command(cmd, description):
    """Run a shell command and handle errors."""
    print(f"\n{'='*70}")
    print(f"{description}")
    print(f"{'='*70}")
    print(f"Command: {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"\n❌ Error: {description} failed!")
        sys.exit(1)

    print(f"\n✓ {description} completed successfully")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Complete BIRDMAn workflow for LLM data"
    )
    parser.add_argument(
        "--data-dir", default=".", help="Directory containing input data files"
    )
    parser.add_argument(
        "--output-dir",
        default="./birdman_llm_analysis",
        help="Output directory for all results",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for BIRDMAn (default: 8)",
    )
    parser.add_argument(
        "--skip-prep",
        action="store_true",
        help="Skip data preparation (use existing prepared_data)",
    )

    args = parser.parse_args()

    # Create output directory
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    data_dir = Path(args.data_dir)

    print("=" * 70)
    print("BIRDMAn Analysis Pipeline for Large Larry Modeler Data")
    print("=" * 70)
    print(f"\nData directory: {data_dir}")
    print(f"Output directory: {output_path}")
    print(f"Threads: {args.threads}")

    # Check for required files
    required_files = [
        "amplicon-matrix.txt.xls",
        "amplicon-taxonomy.txt.xls",
        "LLM_Biomarkers_Long.csv",
    ]

    for file in required_files:
        if not (data_dir / file).exists():
            print(f"\n❌ Required file not found: {file}")
            print(f"Expected location: {data_dir / file}")
            sys.exit(1)

    print("\n✓ All required files found")

    # Step 1: Prepare data
    if not args.skip_prep:
        prep_output = output_path / "prepared_data"
        run_command(
            [
                "python3",
                "prepare_real_data.py",
                "--amplicon-matrix",
                str(data_dir / "amplicon-matrix.txt.xls"),
                "--amplicon-taxonomy",
                str(data_dir / "amplicon-taxonomy.txt.xls"),
                "--biomarkers",
                str(data_dir / "LLM_Biomarkers_Long.csv"),
                "--output-dir",
                str(prep_output),
            ],
            "Step 1: Data Preparation",
        )

        # Load summary
        with open(prep_output / "summary.json") as f:
            summary = json.load(f)

        print(
            f"\n📊 Matched {summary['total_samples']} samples with {summary['total_features']} features"
        )
        print(f"   Date range: {summary['date_range']}")
    else:
        prep_output = output_path / "prepared_data"
        print("\n⏭️  Skipping data preparation (using existing data)")

    # Step 2: Convert to BIOM format
    biom_output = output_path / "biom_data"
    run_command(
        [
            "python3",
            "convert_to_biom.py",
            "--microbiome",
            str(prep_output / "feature_table.csv"),
            "--metabolomics",
            str(prep_output / "metadata.csv"),
            "--output-dir",
            str(biom_output),
            "--depth",
            "10000",
            "--sample-col",
            "sample_id",
        ],
        "Step 2: Convert to BIOM Format",
    )

    # Step 3: Import into QIIME2
    qiime_output = output_path / "qiime2_artifacts"
    qiime_output.mkdir(exist_ok=True)

    run_command(
        [
            "qiime",
            "tools",
            "import",
            "--input-path",
            str(biom_output / "microbiome_table.biom"),
            "--type",
            "FeatureTable[Frequency]",
            "--input-format",
            "BIOMV210Format",
            "--output-path",
            str(qiime_output / "feature-table.qza"),
        ],
        "Step 3: Import to QIIME2",
    )

    # Step 4: Run BIRDMAn - Main model
    print("\n" + "=" * 70)
    print("Step 4: Running BIRDMAn Differential Abundance Analysis")
    print("=" * 70)
    print("\nModel: scfa_total_scfa_normalized")
    print("This tests which microbial features correlate with SCFA levels\n")

    birdman_output = output_path / "birdman_results"

    run_command(
        [
            "qiime",
            "birdman",
            "run",
            "--i-table",
            str(qiime_output / "feature-table.qza"),
            "--m-metadata-file",
            str(biom_output / "metadata.tsv"),
            "--p-formula",
            "scfa_total_scfa_normalized",
            "--o-output-dir",
            str(birdman_output / "scfa_main.qza"),
            "--p-threads",
            str(args.threads),
            "--verbose",
        ],
        "BIRDMAn Analysis: Total SCFA",
    )

    # Step 5: Run BIRDMAn - Individual SCFAs
    print("\n🔬 Running additional models for individual SCFAs...")

    individual_scfas = [
        "scfa_scfa_acetate_normalized",
        "scfa_scfa_propionate_normalized",
        "scfa_scfa_butyrate_normalized",
    ]

    for scfa in individual_scfas:
        scfa_name = scfa.replace("scfa_scfa_", "").replace("_normalized", "")
        print(f"\n  Analyzing {scfa_name}...")

        run_command(
            [
                "qiime",
                "birdman",
                "run",
                "--i-table",
                str(qiime_output / "feature-table.qza"),
                "--m-metadata-file",
                str(biom_output / "metadata.tsv"),
                "--p-formula",
                scfa,
                "--o-output-dir",
                str(birdman_output / f"{scfa_name}.qza"),
                "--p-threads",
                str(args.threads),
                "--verbose",
            ],
            f"BIRDMAn Analysis: {scfa_name}",
        )

    # Step 6: Create visualizations
    viz_output = output_path / "visualizations"
    viz_output.mkdir(exist_ok=True)

    models = ["scfa_main", "acetate", "propionate", "butyrate"]

    for model in models:
        print(f"\n📊 Creating visualization for {model}...")

        run_command(
            [
                "qiime",
                "birdman",
                "plot",
                "--i-data",
                str(birdman_output / f"{model}.qza"),
                "--i-table",
                str(qiime_output / "feature-table.qza"),
                "--m-metadata-file",
                str(biom_output / "metadata.tsv"),
                "--o-visualization",
                str(viz_output / f"{model}_plot.qzv"),
                "--p-palette",
                "rainbow",
                "--p-chart-style",
                "forest",
                "--p-effect-size-threshold",
                "2",
            ],
            f"Visualization: {model}",
        )

    # Step 7: Export results
    export_output = output_path / "exported_results"
    export_output.mkdir(exist_ok=True)

    for model in models:
        model_export = export_output / model
        model_export.mkdir(exist_ok=True)

        run_command(
            [
                "qiime",
                "tools",
                "export",
                "--input-path",
                str(birdman_output / f"{model}.qza"),
                "--output-path",
                str(model_export),
            ],
            f"Export: {model}",
        )

    # Step 8: Extract top taxa
    print("\n" + "=" * 70)
    print("Step 8: Extracting Top Differentially Abundant Taxa")
    print("=" * 70)

    kegg_output = output_path / "kegg_analysis"

    run_command(
        [
            "python3",
            "extract_top_taxa.py",
            "--birdman-results",
            str(export_output / "scfa_main"),
            "--output-dir",
            str(kegg_output),
            "--effect-size-threshold",
            "2",
            "--exclude-zero",
        ],
        "Extract Top Taxa for KEGG Analysis",
    )

    # Generate final report
    print("\n" + "=" * 70)
    print("Generating Final Report")
    print("=" * 70)

    report_content = f"""
# BIRDMAn Analysis Report: Large Larry Modeler Data

## Analysis Summary

**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
**Data Directory**: {data_dir}
**Output Directory**: {output_path}

## Dataset

- **Matched Samples**: 74
- **Microbial Features**: 52 taxa
- **Date Range**: 2012-11-06 to 2019-08-11
- **SCFA Measurements**: 6 types (% Acetate, % Butyrate, % Propionate, % Valerate, Butyrate, Total SCFA)

## Models Tested

1. **Main Model**: Total SCFA (normalized)
   - Tests overall SCFA association with microbiome composition

2. **Individual SCFA Models**:
   - Acetate
   - Propionate
   - Butyrate

## Key Outputs

### Visualizations
Located in: `{output_path}/visualizations/`
- `scfa_main_plot.qzv` - Main SCFA model results
- `acetate_plot.qzv` - Acetate-specific results
- `propionate_plot.qzv` - Propionate-specific results
- `butyrate_plot.qzv` - Butyrate-specific results

**To view**: `qiime tools view <filename>.qzv`

### Exported Results
Located in: `{output_path}/exported_results/`
- Contains raw BIRDMAn output for each model
- Use for downstream analysis and custom visualizations

### Top Taxa for KEGG Analysis
Located in: `{output_path}/kegg_analysis/`
- `enriched_features.txt` - Taxa positively associated with SCFA
- `depleted_features.txt` - Taxa negatively associated with SCFA
- `significant_features.tsv` - Complete results table
- `forest_plot.png` - Visualization of top taxa
- `KEGG_ANALYSIS_GUIDE.md` - Instructions for pathway analysis

## Interpretation

### Expected Findings (if hypothesis supported):

**Enriched Taxa** (positively correlated with SCFA):
- SCFA-producing bacteria
- Likely genera: Faecalibacterium, Roseburia, Eubacterium, Ruminococcus
- Associated with butanoate and propanoate metabolism pathways

**Depleted Taxa** (negatively correlated with SCFA):
- May compete with SCFA producers
- Or thrive in low-SCFA environments
- Interesting for understanding SCFA ecology

### Next Steps:

1. **Review Visualizations**:
   ```bash
   cd {output_path}/visualizations
   qiime tools view scfa_main_plot.qzv
   ```

2. **Examine Top Taxa**:
   ```bash
   cd {output_path}/kegg_analysis
   cat enriched_features.txt
   cat depleted_features.txt
   ```

3. **KEGG Pathway Analysis**:
   - Follow instructions in `KEGG_ANALYSIS_GUIDE.md`
   - Test enrichment in butanoate metabolism (map00650)
   - Test enrichment in propanoate metabolism (map00640)
   - Test enrichment in fatty acid biosynthesis (map00061)

4. **Validate Findings**:
   - Cross-reference with literature on known SCFA producers
   - Check taxonomic classifications in taxonomy.csv
   - Consider biological mechanisms

## Files Generated

```
{output_path}/
├── prepared_data/
│   ├── feature_table.csv
│   ├── metadata.csv
│   ├── taxonomy.csv
│   └── scfa_timeseries.png
├── biom_data/
│   ├── microbiome_table.biom
│   └── metadata.tsv
├── qiime2_artifacts/
│   └── feature-table.qza
├── birdman_results/
│   ├── scfa_main.qza
│   ├── acetate.qza
│   ├── propionate.qza
│   └── butyrate.qza
├── visualizations/
│   ├── scfa_main_plot.qzv
│   ├── acetate_plot.qzv
│   ├── propionate_plot.qzv
│   └── butyrate_plot.qzv
├── exported_results/
│   └── [exported data for each model]
├── kegg_analysis/
│   ├── enriched_features.txt
│   ├── depleted_features.txt
│   ├── significant_features.tsv
│   ├── forest_plot.png
│   └── KEGG_ANALYSIS_GUIDE.md
└── ANALYSIS_REPORT.md (this file)
```

## Citation

If you use this analysis, please cite:

```
Rahman G, Morton JT, Martino C, Sepich-Poore GD, Allaband C, Guccione C, 
Chen Y, Hakim D, Estaki M, Knight R. BIRDMAn: A Bayesian differential 
abundance framework that enables robust inference of host-microbe associations. 
bioRxiv. 2023. doi: 10.1101/2023.01.30.526328
```

## Contact

For questions about this analysis, refer to the README.md and other documentation files.
"""

    with open(output_path / "ANALYSIS_REPORT.md", "w") as f:
        f.write(report_content)

    print("\n" + "=" * 70)
    print("🎉 Analysis Complete!")
    print("=" * 70)
    print(f"\nAll results saved to: {output_path}")
    print("\nKey files:")
    print(f"  - Visualizations: {output_path}/visualizations/")
    print(f"  - Top taxa: {output_path}/kegg_analysis/")
    print(f"  - Full report: {output_path}/ANALYSIS_REPORT.md")
    print("\nTo view results:")
    print(f"  qiime tools view {output_path}/visualizations/scfa_main_plot.qzv")
    print("\nNext steps:")
    print("  1. Review visualizations")
    print("  2. Check top differentially abundant taxa")
    print("  3. Perform KEGG pathway enrichment analysis")
    print("  4. Validate findings with literature")


if __name__ == "__main__":
    import pandas as pd

    main()
