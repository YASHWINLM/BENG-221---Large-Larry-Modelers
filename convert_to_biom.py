#!/usr/bin/env python3
"""
Convert Excel metabolomics and microbiome data to BIOM format for BIRDMAn analysis.

This script:
1. Reads Excel files with metabolomics (SCFA) and microbiome data
2. Converts relative abundance to counts (if needed)
3. Creates BIOM tables compatible with QIIME2/BIRDMAn
4. Generates metadata file with SCFA and covariate information
"""

import pandas as pd
import numpy as np
import biom
from biom.table import Table
import argparse
import json
from pathlib import Path


def relative_to_counts(relative_abundance, depth=10000):
    """
    Convert relative abundance values to count data.

    This replicates common R functionality for converting proportions to counts:
    counts = round(relative_abundance * depth)

    Parameters:
    -----------
    relative_abundance : array-like
        Relative abundance values (proportions summing to 1)
    depth : int
        Sequencing depth to use for conversion (default: 10000)

    Returns:
    --------
    numpy.ndarray
        Integer counts
    """
    # Ensure data is numeric
    rel_abund = np.array(relative_abundance, dtype=float)

    # Normalize to sum to 1 if not already
    if not np.isclose(rel_abund.sum(), 1.0, atol=0.01):
        print(
            f"Warning: Values don't sum to 1 (sum={rel_abund.sum():.4f}). Normalizing..."
        )
        rel_abund = rel_abund / rel_abund.sum()

    # Convert to counts
    counts = np.round(rel_abund * depth).astype(int)

    # Adjust for rounding errors to maintain exact depth
    diff = depth - counts.sum()
    if diff != 0:
        # Add/subtract from the largest count
        max_idx = np.argmax(counts)
        counts[max_idx] += diff

    return counts


def read_microbiome_data(filepath, sample_col="sample_id", depth=10000):
    """
    Read microbiome data from Excel or CSV file.

    Parameters:
    -----------
    filepath : str
        Path to Excel or CSV file with microbiome data
    sample_col : str
        Column name containing sample IDs
    depth : int
        Sequencing depth for count conversion

    Returns:
    --------
    tuple
        (feature_table_df, sample_ids, feature_ids)
    """
    # Auto-detect file format
    if filepath.endswith(".csv"):
        df = pd.read_csv(filepath, index_col=0)  # Sample IDs are in first column/index
    else:
        df = pd.read_excel(filepath)

    # Check if sample_col exists as a column
    if sample_col in df.columns:
        # Extract sample IDs from column
        sample_ids = df[sample_col].values
        # Get feature columns (all columns except sample_col)
        feature_cols = [col for col in df.columns if col != sample_col]
        # Extract feature data
        feature_data = df[feature_cols]
    else:
        # Sample IDs are in the index
        sample_ids = df.index.values
        # All columns are features
        feature_cols = df.columns.tolist()
        feature_data = df

    # Check if data is relative abundance or already counts
    is_relative = all(
        (feature_data.sum(axis=1) > 0.9) & (feature_data.sum(axis=1) <= 1.1)
    )

    if is_relative:
        print("Data appears to be relative abundance. Converting to counts...")
        count_data = np.array(
            [
                relative_to_counts(row.values, depth=depth)
                for _, row in feature_data.iterrows()
            ]
        )
        feature_table = pd.DataFrame(count_data, columns=feature_cols, index=sample_ids)
    else:
        print("Data appears to be counts already.")
        feature_table = feature_data.copy()
        feature_table.index = sample_ids

    return (
        feature_table.T,
        sample_ids,
        feature_cols,
    )  # Transpose: features as rows, samples as columns


def read_metabolomics_data(filepath, sample_col="sample_id", scfa_cols=None):
    """
    Read metabolomics data from Excel or CSV file.

    Parameters:
    -----------
    filepath : str
        Path to Excel or CSV file with metabolomics data
    sample_col : str
        Column name containing sample IDs
    scfa_cols : list
        List of SCFA column names (if None, auto-detect)

    Returns:
    --------
    tuple
        (metabolite_table_df, sample_ids, metabolite_ids)
    """
    # Auto-detect file format
    if filepath.endswith(".csv"):
        df = pd.read_csv(filepath)
    else:
        df = pd.read_excel(filepath)

    # Check if sample_col exists
    if sample_col in df.columns:
        # Extract sample IDs from column
        sample_ids = df[sample_col].values
    elif "Unnamed: 0" in df.columns:
        # Sample IDs might be in first unnamed column
        sample_ids = df["Unnamed: 0"].values
        sample_col = "Unnamed: 0"
    else:
        # Use index as sample IDs
        sample_ids = df.index.values
        # Create a temporary column
        df["__sample_id__"] = sample_ids
        sample_col = "__sample_id__"

    # Auto-detect SCFA columns if not provided
    if scfa_cols is None:
        # Common SCFA names
        common_scfas = [
            "acetate",
            "propionate",
            "butyrate",
            "valerate",
            "isobutyrate",
            "isovalerate",
            "total_scfa",
            "scfa",
        ]
        scfa_cols = [
            col
            for col in df.columns
            if any(scfa.lower() in col.lower() for scfa in common_scfas)
        ]

    if not scfa_cols:
        # If no SCFA columns found, use all numeric columns except sample_col
        scfa_cols = [
            col
            for col in df.columns
            if col != sample_col and pd.api.types.is_numeric_dtype(df[col])
        ]

    print(f"Using metabolite columns: {scfa_cols}")

    # Extract metabolite data
    metabolite_data = df[scfa_cols]
    metabolite_data.index = sample_ids

    return (
        metabolite_data.T,
        sample_ids,
        scfa_cols,
    )  # Transpose: metabolites as rows, samples as columns


def create_biom_table(
    feature_table_df, observation_metadata=None, sample_metadata=None
):
    """
    Create a BIOM table from a pandas DataFrame.

    Parameters:
    -----------
    feature_table_df : pd.DataFrame
        Feature table with features as rows and samples as columns
    observation_metadata : list of dict
        Metadata for each feature/observation
    sample_metadata : list of dict
        Metadata for each sample

    Returns:
    --------
    biom.Table
        BIOM table object
    """
    # Convert to numpy array
    data = feature_table_df.values.astype(float)

    # Get observation (feature) and sample IDs
    observation_ids = feature_table_df.index.tolist()
    sample_ids = feature_table_df.columns.tolist()

    # Create BIOM table
    table = Table(
        data,
        observation_ids,
        sample_ids,
        observation_metadata=observation_metadata,
        sample_metadata=sample_metadata,
    )

    return table


def prepare_metadata_file(
    microbiome_samples,
    metabolomics_df,
    covariates_df,
    output_path,
    scfa_col="total_scfa",
):
    """
    Create a QIIME2-compatible metadata file.

    Parameters:
    -----------
    microbiome_samples : list
        List of sample IDs from microbiome data
    metabolomics_df : pd.DataFrame
        Metabolomics data with SCFA values
    covariates_df : pd.DataFrame
        Covariate data (age, BMI, cholesterol, etc.)
    output_path : str
        Path to save metadata TSV file
    scfa_col : str
        Column name for total SCFA or the SCFA of interest
    """
    # Create metadata dataframe
    metadata = pd.DataFrame({"sample_id": microbiome_samples})
    metadata.set_index("sample_id", inplace=True)

    # Add SCFA data
    if scfa_col in metabolomics_df.columns:
        metadata["scfa_level"] = metabolomics_df[scfa_col]
    else:
        # If no total_scfa, calculate it
        scfa_cols = [
            col
            for col in metabolomics_df.columns
            if any(s in col.lower() for s in ["acetate", "propionate", "butyrate"])
        ]
        if scfa_cols:
            metadata["scfa_level"] = metabolomics_df[scfa_cols].sum(axis=1)

    # Normalize SCFA levels (z-score normalization)
    if "scfa_level" in metadata.columns:
        metadata["scfa_level_normalized"] = (
            metadata["scfa_level"] - metadata["scfa_level"].mean()
        ) / metadata["scfa_level"].std()

    # Add covariates if available
    if covariates_df is not None:
        # Common covariate columns
        covariate_cols = [
            "age",
            "time",
            "body_weight",
            "bmi",
            "cholesterol",
            "triglycerides",
            "sex",
            "host_subject_id",
        ]

        for col in covariate_cols:
            if col in covariates_df.columns:
                metadata[col] = covariates_df[col]

    # Write metadata file with QIIME2 format
    with open(output_path, "w") as f:
        # Write header comment
        f.write("# QIIME2 metadata file for BIRDMAn analysis\n")
        f.write("#\n")
        # Write column names with first column as 'sample-id'
        f.write("sample-id\t" + "\t".join(metadata.columns) + "\n")

        # Write data
        for sample_id, row in metadata.iterrows():
            f.write(f"{sample_id}\t" + "\t".join(map(str, row.values)) + "\n")

    print(f"Metadata file saved to: {output_path}")
    return metadata


def main():
    parser = argparse.ArgumentParser(
        description="Convert Excel data to BIOM format for BIRDMAn analysis"
    )
    parser.add_argument(
        "--microbiome", required=True, help="Path to microbiome Excel file"
    )
    parser.add_argument(
        "--metabolomics", required=True, help="Path to metabolomics/SCFA Excel file"
    )
    parser.add_argument("--covariates", help="Path to covariates Excel file (optional)")
    parser.add_argument(
        "--output-dir", default="./biom_output", help="Output directory for BIOM files"
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=10000,
        help="Sequencing depth for count conversion (default: 10000)",
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
        help="Column name for sample IDs (default: sample_id)",
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Converting data to BIOM format for BIRDMAn")
    print("=" * 60)

    # Read microbiome data
    print("\n1. Reading microbiome data...")
    micro_table, micro_samples, micro_features = read_microbiome_data(
        args.microbiome, sample_col=args.sample_col, depth=args.depth
    )
    print(f"   - Samples: {len(micro_samples)}")
    print(f"   - Features: {len(micro_features)}")

    # Read metabolomics data
    print("\n2. Reading metabolomics data...")
    metab_table, metab_samples, metab_features = read_metabolomics_data(
        args.metabolomics, sample_col=args.sample_col
    )
    print(f"   - Samples: {len(metab_samples)}")
    print(f"   - Metabolites: {len(metab_features)}")

    # Read covariates if provided
    covariates_df = None
    if args.covariates:
        print("\n3. Reading covariate data...")
        if args.covariates.endswith(".csv"):
            covariates_df = pd.read_csv(args.covariates)
        else:
            covariates_df = pd.read_excel(args.covariates)
        print(f"   - Samples: {len(covariates_df)}")

    # Create BIOM tables
    print("\n4. Creating BIOM tables...")

    # Microbiome BIOM table
    micro_biom = create_biom_table(micro_table)
    micro_biom_path = output_dir / "microbiome_table.biom"
    with biom.util.biom_open(str(micro_biom_path), "w") as f:
        micro_biom.to_hdf5(f, "BIRDMAn microbiome analysis")
    print(f"   - Microbiome BIOM: {micro_biom_path}")

    # Metabolomics BIOM table
    metab_biom = create_biom_table(metab_table)
    metab_biom_path = output_dir / "metabolomics_table.biom"
    with biom.util.biom_open(str(metab_biom_path), "w") as f:
        metab_biom.to_hdf5(f, "BIRDMAn metabolomics analysis")
    print(f"   - Metabolomics BIOM: {metab_biom_path}")

    # Create metadata file
    print("\n5. Creating metadata file...")
    if args.metabolomics.endswith(".csv"):
        metabolomics_df = pd.read_csv(args.metabolomics)
    else:
        metabolomics_df = pd.read_excel(args.metabolomics)

    # Handle sample IDs in index vs column
    if "sample_id" not in metabolomics_df.columns:
        if metabolomics_df.index.name == "sample_id" or "sample_id" in str(
            metabolomics_df.index.name
        ):
            metabolomics_df = metabolomics_df.reset_index()
        else:
            metabolomics_df["sample_id"] = metabolomics_df.index
            metabolomics_df = metabolomics_df.reset_index(drop=True)

    metabolomics_df.set_index("sample_id", inplace=True)

    metadata_path = output_dir / "metadata.tsv"
    prepare_metadata_file(micro_samples, metabolomics_df, covariates_df, metadata_path)

    # Save summary
    summary = {
        "microbiome": {
            "samples": len(micro_samples),
            "features": len(micro_features),
            "biom_path": str(micro_biom_path),
        },
        "metabolomics": {
            "samples": len(metab_samples),
            "features": len(metab_features),
            "biom_path": str(metab_biom_path),
        },
        "metadata_path": str(metadata_path),
    }

    summary_path = output_dir / "conversion_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("Conversion complete!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print(f"Summary: {summary_path}")
    print("\nNext steps:")
    print("1. Import BIOM tables into QIIME2:")
    print(f"   qiime tools import \\")
    print(f"     --input-path {micro_biom_path} \\")
    print(f"     --type 'FeatureTable[Frequency]' \\")
    print(f"     --output-path feature-table.qza")
    print("\n2. Run BIRDMAn analysis:")
    print(f"   qiime birdman run \\")
    print(f"     --i-table feature-table.qza \\")
    print(f"     --m-metadata-file {metadata_path} \\")
    print(f"     --p-formula 'scfa_level_normalized+time+body_weight+cholesterol' \\")
    print(f"     --o-output-dir birdman_results.qza \\")
    print(f"     --p-threads 32")


if __name__ == "__main__":
    main()
