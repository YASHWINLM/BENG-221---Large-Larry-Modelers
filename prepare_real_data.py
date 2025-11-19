#!/usr/bin/env python3
"""
Prepare Large Larry Modeler data for BIRDMAn analysis.

This script processes YOUR actual data files:
- amplicon-matrix_txt.xls: 52 taxa x 156 samples (2012-2019)
- LLM_Biomarkers_Long.csv: SCFA and biomarker data (184 dates with SCFA)
- Matches samples by date to create paired microbiome-metabolomics datasets
"""

import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re


def parse_sample_date(sample_id):
    """
    Parse date from sample ID like 'LLM.2012.11.06' -> '2012-11-06'
    """
    parts = sample_id.replace('LLM.', '').split('.')
    if len(parts) == 3:
        year, month, day = parts
        return f"{year}-{month.zfill(2)}-{day.zfill(2)}"
    return None


def load_amplicon_data(amplicon_file, taxonomy_file):
    """
    Load amplicon microbiome data.
    """
    print("\n1. Loading amplicon data...")
    
    # Load feature table
    feature_table = pd.read_csv(amplicon_file, sep='\t', index_col=0)
    print(f"   Amplicon shape: {feature_table.shape}")
    print(f"   Features (taxa): {feature_table.shape[0]}")
    print(f"   Samples: {feature_table.shape[1]}")
    
    # Load taxonomy
    taxonomy = pd.read_csv(taxonomy_file, sep='\t', index_col=0)
    print(f"   Taxonomy loaded: {taxonomy.shape[0]} taxa")
    
    # Parse dates from sample IDs
    sample_dates = {col: parse_sample_date(col) for col in feature_table.columns}
    date_range = [d for d in sample_dates.values() if d]
    print(f"   Date range: {min(date_range)} to {max(date_range)}")
    
    return feature_table, taxonomy, sample_dates


def load_scfa_data(biomarker_file):
    """
    Load SCFA data from LLM_Biomarkers_Long.csv
    """
    print("\n2. Loading SCFA biomarker data...")
    
    df = pd.read_csv(biomarker_file)
    
    # Extract SCFA measurements
    scfa_keywords = ['SCFA', 'Butyrate', 'Propionate', 'Acetate', 'Valerate']
    scfa_data = df[df['Biomarker'].str.contains('|'.join(scfa_keywords), case=False, na=False)].copy()
    
    # Convert Value column to numeric
    scfa_data['Value'] = pd.to_numeric(scfa_data['Value'], errors='coerce')
    
    print(f"   Total biomarker records: {len(df)}")
    print(f"   SCFA records: {len(scfa_data)}")
    print(f"   SCFA biomarkers: {scfa_data['Biomarker'].unique()}")
    
    # Pivot to wide format
    scfa_wide = scfa_data.pivot_table(
        index='Date',
        columns='Biomarker',
        values='Value',
        aggfunc='mean'
    ).reset_index()
    
    print(f"   Unique dates with SCFA: {len(scfa_wide)}")
    
    # Calculate total SCFA if individual components available
    scfa_cols = [col for col in scfa_wide.columns 
                if any(s in col for s in ['Acetate', 'Propionate', 'Butyrate']) 
                and '%' not in col]
    
    if scfa_cols and 'Total SCFA' not in scfa_wide.columns:
        scfa_wide['Total SCFA'] = scfa_wide[scfa_cols].sum(axis=1, min_count=1)
        print(f"   Calculated Total SCFA from: {scfa_cols}")
    
    return scfa_wide


def load_additional_biomarkers(biomarker_file, selected_biomarkers=None):
    """
    Load additional biomarkers (cholesterol, body weight, etc.)
    """
    print("\n3. Loading additional biomarkers...")
    
    df = pd.read_csv(biomarker_file)
    
    if selected_biomarkers is None:
        # Common biomarkers for covariates
        selected_biomarkers = [
            'Sum Total Cholesterol',
            'LDL',
            'HDL',
            'Triglycerides',
            'Body Weight',
            'BMI',
            'Fasting Glucose'
        ]
    
    # Filter for selected biomarkers
    bio_data = df[df['Biomarker'].isin(selected_biomarkers)].copy()
    
    # Convert Value column to numeric
    bio_data['Value'] = pd.to_numeric(bio_data['Value'], errors='coerce')
    
    print(f"   Selected biomarkers: {selected_biomarkers}")
    print(f"   Found biomarkers: {bio_data['Biomarker'].unique().tolist()}")
    
    # Pivot to wide format
    bio_wide = bio_data.pivot_table(
        index='Date',
        columns='Biomarker',
        values='Value',
        aggfunc='mean'
    ).reset_index()
    
    print(f"   Dates with biomarkers: {len(bio_wide)}")
    
    return bio_wide


def match_samples(feature_table, sample_dates, scfa_data, bio_data):
    """
    Match microbiome samples with SCFA and biomarker data by date.
    """
    print("\n4. Matching samples by date...")
    
    # Get dates from each dataset
    microbiome_dates = {date: sample_id for sample_id, date in sample_dates.items() if date}
    scfa_dates_set = set(scfa_data['Date'].values)
    
    # Find overlapping dates
    overlap_dates = set(microbiome_dates.keys()) & scfa_dates_set
    
    print(f"   Microbiome samples: {len(microbiome_dates)}")
    print(f"   SCFA dates: {len(scfa_dates_set)}")
    print(f"   Overlapping dates: {len(overlap_dates)}")
    
    if len(overlap_dates) == 0:
        print("\n   ⚠️  No direct date matches found!")
        print("   This might be because:")
        print("   - Date formats don't match")
        print("   - Samples were taken on different days")
        print("   - Need to match within a time window")
        return None, None, None
    
    # Create matched dataset
    matched_samples = []
    matched_metadata = []
    
    for date in sorted(overlap_dates):
        sample_id = microbiome_dates[date]
        
        # Get SCFA data for this date
        scfa_row = scfa_data[scfa_data['Date'] == date].iloc[0]
        
        # Get biomarker data if available
        bio_row = bio_data[bio_data['Date'] == date]
        if len(bio_row) > 0:
            bio_row = bio_row.iloc[0]
        else:
            bio_row = None
        
        matched_samples.append(sample_id)
        
        # Create metadata row
        meta_row = {'sample_id': sample_id, 'date': date}
        
        # Add SCFA data
        for col in scfa_data.columns:
            if col != 'Date':
                meta_row[f'scfa_{col.lower().replace(" ", "_").replace("%_", "")}'] = scfa_row[col]
        
        # Add biomarker data if available
        if bio_row is not None:
            for col in bio_data.columns:
                if col != 'Date':
                    meta_row[col.lower().replace(" ", "_")] = bio_row[col]
        
        matched_metadata.append(meta_row)
    
    # Create metadata dataframe
    metadata_df = pd.DataFrame(matched_metadata)
    
    # Filter feature table to matched samples
    matched_feature_table = feature_table[matched_samples].T
    
    print(f"\n   Matched samples: {len(matched_samples)}")
    print(f"   Date range: {min(overlap_dates)} to {max(overlap_dates)}")
    
    return matched_feature_table, metadata_df, list(overlap_dates)


def prepare_for_birdman(feature_table, metadata_df):
    """
    Prepare data for BIRDMAn analysis.
    """
    print("\n5. Preparing for BIRDMAn...")
    
    # Normalize SCFA values (z-score)
    scfa_cols = [col for col in metadata_df.columns if 'scfa' in col.lower()]
    for col in scfa_cols:
        if metadata_df[col].notna().sum() > 0:
            mean_val = metadata_df[col].mean()
            std_val = metadata_df[col].std()
            if std_val > 0:
                metadata_df[f'{col}_normalized'] = (metadata_df[col] - mean_val) / std_val
                print(f"   Normalized {col}")
    
    # Check data types
    print("\n   Checking covariate data types...")
    numeric_covariates = []
    for col in metadata_df.columns:
        if col not in ['sample_id', 'date']:
            if pd.api.types.is_numeric_dtype(metadata_df[col]):
                numeric_covariates.append(col)
                non_na = metadata_df[col].notna().sum()
                print(f"     {col}: {non_na}/{len(metadata_df)} values ({100*non_na/len(metadata_df):.1f}%)")
    
    return feature_table, metadata_df


def create_visualizations(matched_dates, scfa_data, metadata_df, output_dir):
    """
    Create visualizations of the matched data.
    """
    print("\n6. Creating visualizations...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: SCFA over time
    fig, ax = plt.subplots(figsize=(14, 6))
    
    scfa_cols = [col for col in metadata_df.columns if 'scfa' in col.lower() and 'normalized' not in col.lower()]
    for col in scfa_cols:
        if metadata_df[col].notna().sum() > 2:
            dates = pd.to_datetime(metadata_df['date'])
            ax.plot(dates, metadata_df[col], marker='o', label=col, alpha=0.7)
    
    ax.set_xlabel('Date')
    ax.set_ylabel('SCFA Concentration')
    ax.set_title('SCFA Levels Over Time (Matched Samples)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path / 'scfa_timeseries.png', dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}/scfa_timeseries.png")
    
    # Plot 2: Sample availability timeline
    fig, ax = plt.subplots(figsize=(14, 4))
    
    dates = pd.to_datetime(matched_dates)
    ax.scatter(dates, [1]*len(dates), s=100, alpha=0.6, color='green')
    ax.set_xlabel('Date')
    ax.set_yticks([])
    ax.set_title(f'Matched Microbiome-SCFA Samples Timeline (n={len(matched_dates)})')
    ax.grid(True, alpha=0.3, axis='x')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path / 'sample_timeline.png', dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}/sample_timeline.png")
    
    plt.close('all')


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Large Larry Modeler data for BIRDMAn analysis"
    )
    parser.add_argument(
        "--amplicon-matrix",
        default="amplicon-matrix_txt.xls",
        help="Path to amplicon matrix file"
    )
    parser.add_argument(
        "--amplicon-taxonomy",
        default="amplicon-taxonomy_txt.xls",
        help="Path to amplicon taxonomy file"
    )
    parser.add_argument(
        "--biomarkers",
        default="LLM_Biomarkers_Long.csv",
        help="Path to biomarkers file"
    )
    parser.add_argument(
        "--output-dir",
        default="./prepared_data_llm",
        help="Output directory"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Large Larry Modeler Data Preparation for BIRDMAn")
    print("=" * 70)
    
    # Load data
    feature_table, taxonomy, sample_dates = load_amplicon_data(
        args.amplicon_matrix,
        args.amplicon_taxonomy
    )
    
    scfa_data = load_scfa_data(args.biomarkers)
    
    bio_data = load_additional_biomarkers(args.biomarkers)
    
    # Match samples
    matched_feature_table, metadata_df, matched_dates = match_samples(
        feature_table, sample_dates, scfa_data, bio_data
    )
    
    if matched_feature_table is None:
        print("\n❌ No matched samples found. Cannot proceed with analysis.")
        return
    
    # Prepare for BIRDMAn
    matched_feature_table, metadata_df = prepare_for_birdman(
        matched_feature_table, metadata_df
    )
    
    # Save outputs
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n7. Saving outputs...")
    
    # Save feature table
    matched_feature_table.to_csv(output_path / "feature_table.csv")
    print(f"   Saved: {output_path}/feature_table.csv")
    
    # Save metadata
    metadata_df.to_csv(output_path / "metadata.csv", index=False)
    print(f"   Saved: {output_path}/metadata.csv")
    
    # Save taxonomy
    # Filter taxonomy to features present in matched samples
    features_present = matched_feature_table.columns
    taxonomy_filtered = taxonomy[taxonomy.index.isin(features_present)]
    taxonomy_filtered.to_csv(output_path / "taxonomy.csv")
    print(f"   Saved: {output_path}/taxonomy.csv")
    
    # Create visualizations
    create_visualizations(matched_dates, scfa_data, metadata_df, args.output_dir)
    
    # Save summary
    summary = {
        "total_features": len(matched_feature_table.columns),
        "total_samples": len(matched_feature_table),
        "date_range": f"{min(matched_dates)} to {max(matched_dates)}",
        "scfa_columns": [col for col in metadata_df.columns if 'scfa' in col.lower()],
        "covariate_columns": [col for col in metadata_df.columns 
                            if col not in ['sample_id', 'date'] and 'scfa' not in col.lower()],
        "matched_dates": sorted(matched_dates)
    }
    
    import json
    with open(output_path / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "=" * 70)
    print("Data preparation complete!")
    print("=" * 70)
    print(f"\nOutput directory: {args.output_dir}")
    print(f"\nMatched samples: {len(matched_feature_table)}")
    print(f"Features: {len(matched_feature_table.columns)}")
    print(f"\nNext steps:")
    print(f"1. Review metadata.csv to check SCFA and covariate data")
    print(f"2. Run convert_to_biom.py to create BIOM tables")
    print(f"3. Run BIRDMAn analysis")


if __name__ == "__main__":
    main()
