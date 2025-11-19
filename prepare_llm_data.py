#!/usr/bin/env python3
"""
Data Preparation for BENG-221 Large Larry Modelers BIRDMAn Analysis

This script specifically handles the LLM_Biomarkers_Long.csv file and 
prepares it for BIRDMAn analysis by:
1. Filtering values that meet "good range" criteria
2. Harmonizing biomarkers across timepoints
3. Identifying timepoints with sufficient data coverage
4. Matching microbiome and metabolomics samples
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse


def load_biomarker_data(filepath):
    """Load the LLM_Biomarkers_Long.csv file."""
    df = pd.read_csv(filepath)
    print(f"Loaded biomarker data: {df.shape[0]} rows, {df.shape[1]} columns")
    print(f"Columns: {df.columns.tolist()}")
    return df


def filter_good_range(df, range_column='good_range'):
    """
    Filter biomarkers to only include values in the "good range".
    
    Parameters:
    -----------
    df : pd.DataFrame
        Biomarker dataframe
    range_column : str
        Column indicating if value is in good range
    """
    if range_column in df.columns:
        good_data = df[df[range_column] == True].copy()
        print(f"\nFiltering by good range:")
        print(f"  Original rows: {len(df)}")
        print(f"  Good range rows: {len(good_data)}")
        print(f"  Removed: {len(df) - len(good_data)} ({100*(len(df)-len(good_data))/len(df):.1f}%)")
        return good_data
    else:
        print(f"Warning: '{range_column}' column not found. Returning all data.")
        return df


def calculate_biomarker_coverage(df, biomarker_col='biomarker', 
                                 timepoint_col='day', value_col='value'):
    """
    Calculate data coverage for each biomarker across timepoints.
    """
    coverage = df.groupby(biomarker_col)[timepoint_col].agg([
        ('n_timepoints', 'nunique'),
        ('n_measurements', 'count')
    ]).reset_index()
    
    # Calculate completion rate
    total_timepoints = df[timepoint_col].nunique()
    coverage['completion_rate'] = coverage['n_timepoints'] / total_timepoints * 100
    
    coverage = coverage.sort_values('completion_rate', ascending=False)
    
    print("\nBiomarker Coverage Summary:")
    print("=" * 70)
    print(f"{'Biomarker':<30} {'Timepoints':<12} {'Measurements':<12} {'Coverage'}")
    print("-" * 70)
    for _, row in coverage.iterrows():
        print(f"{row[biomarker_col]:<30} {row['n_timepoints']:<12} "
              f"{row['n_measurements']:<12} {row['completion_rate']:.1f}%")
    print("=" * 70)
    
    return coverage


def identify_scfa_timepoints(df, scfa_biomarkers=None, 
                            biomarker_col='biomarker',
                            timepoint_col='day',
                            min_coverage=0.6):
    """
    Identify timepoints with sufficient SCFA data coverage.
    """
    if scfa_biomarkers is None:
        # Common SCFA names
        scfa_biomarkers = [
            'acetate', 'propionate', 'butyrate', 'valerate',
            'isobutyrate', 'isovalerate', 'total_scfa', 'scfa'
        ]
    
    # Filter for SCFA biomarkers
    scfa_mask = df[biomarker_col].str.lower().str.contains('|'.join(scfa_biomarkers))
    scfa_data = df[scfa_mask].copy()
    
    print(f"\nSCFA Biomarkers found: {scfa_data[biomarker_col].unique().tolist()}")
    
    # Count SCFA measurements per timepoint
    timepoint_coverage = scfa_data.groupby(timepoint_col).agg({
        biomarker_col: 'nunique',
        'value': 'count'
    }).reset_index()
    
    timepoint_coverage.columns = [timepoint_col, 'n_scfa_types', 'n_measurements']
    
    # Identify good timepoints
    max_scfa_types = timepoint_coverage['n_scfa_types'].max()
    good_timepoints = timepoint_coverage[
        timepoint_coverage['n_scfa_types'] >= max_scfa_types * min_coverage
    ][timepoint_col].tolist()
    
    print(f"\nTimepoints with sufficient SCFA coverage (>{min_coverage*100}%):")
    print(f"  Good timepoints: {len(good_timepoints)}")
    print(f"  Days: {sorted(good_timepoints)}")
    
    return good_timepoints, timepoint_coverage


def harmonize_biomarkers(df, timepoints, biomarker_col='biomarker',
                        timepoint_col='day', value_col='value',
                        subject_col='subject_id'):
    """
    Harmonize biomarkers for selected timepoints.
    Handle missing values through interpolation or forward-fill.
    """
    # Filter to selected timepoints
    harmonized = df[df[timepoint_col].isin(timepoints)].copy()
    
    # Pivot to wide format
    wide = harmonized.pivot_table(
        index=[subject_col, timepoint_col],
        columns=biomarker_col,
        values=value_col,
        aggfunc='mean'  # Average if multiple measurements
    ).reset_index()
    
    print(f"\nHarmonized data shape: {wide.shape}")
    print(f"Columns: {wide.columns.tolist()}")
    
    # Check for missing values
    missing_pct = (wide.isnull().sum() / len(wide) * 100).sort_values(ascending=False)
    print("\nMissing data by biomarker:")
    for biomarker, pct in missing_pct.items():
        if pct > 0:
            print(f"  {biomarker}: {pct:.1f}%")
    
    return wide


def match_microbiome_metabolome(microbiome_samples, metabolome_samples):
    """
    Find overlapping samples between microbiome and metabolome data.
    """
    micro_set = set(microbiome_samples)
    metab_set = set(metabolome_samples)
    
    overlap = micro_set.intersection(metab_set)
    micro_only = micro_set - metab_set
    metab_only = metab_set - micro_set
    
    print("\nSample Matching:")
    print(f"  Microbiome samples: {len(micro_set)}")
    print(f"  Metabolome samples: {len(metab_set)}")
    print(f"  Overlapping samples: {len(overlap)}")
    print(f"  Microbiome only: {len(micro_only)}")
    print(f"  Metabolome only: {len(metab_only)}")
    
    if len(overlap) == 0:
        print("\n⚠️  WARNING: No overlapping samples found!")
        print("  Check that sample IDs match between datasets")
    
    return list(overlap)


def calculate_total_scfa(df, scfa_cols=None):
    """
    Calculate total SCFA from individual SCFA measurements.
    """
    if scfa_cols is None:
        # Auto-detect SCFA columns
        scfa_keywords = ['acetate', 'propionate', 'butyrate', 'valerate']
        scfa_cols = [col for col in df.columns 
                    if any(kw in col.lower() for kw in scfa_keywords)
                    and 'total' not in col.lower()]
    
    if scfa_cols:
        df['total_scfa'] = df[scfa_cols].sum(axis=1)
        print(f"\nCalculated total_scfa from: {scfa_cols}")
        print(f"  Mean: {df['total_scfa'].mean():.2f}")
        print(f"  Std: {df['total_scfa'].std():.2f}")
        print(f"  Range: [{df['total_scfa'].min():.2f}, {df['total_scfa'].max():.2f}]")
    else:
        print("\nWarning: No SCFA columns found for total calculation")
    
    return df


def visualize_data_coverage(coverage_df, timepoint_coverage, output_dir):
    """
    Create visualizations of data coverage.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: Biomarker completion rates
    fig, ax = plt.subplots(figsize=(12, 8))
    
    coverage_df_sorted = coverage_df.sort_values('completion_rate')
    colors = ['green' if x >= 60 else 'orange' if x >= 40 else 'red' 
             for x in coverage_df_sorted['completion_rate']]
    
    ax.barh(range(len(coverage_df_sorted)), 
           coverage_df_sorted['completion_rate'],
           color=colors, alpha=0.7)
    ax.set_yticks(range(len(coverage_df_sorted)))
    ax.set_yticklabels(coverage_df_sorted['biomarker'])
    ax.set_xlabel('Completion Rate (%)')
    ax.set_title('Biomarker Data Coverage')
    ax.axvline(x=60, color='black', linestyle='--', linewidth=0.5, label='60% threshold')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path / "biomarker_coverage.png", dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}/biomarker_coverage.png")
    
    # Plot 2: SCFA coverage by timepoint
    if timepoint_coverage is not None:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        ax.plot(timepoint_coverage['day'], 
               timepoint_coverage['n_scfa_types'], 
               marker='o', linewidth=2, markersize=8)
        ax.set_xlabel('Day')
        ax.set_ylabel('Number of SCFA Types')
        ax.set_title('SCFA Data Coverage Across Timepoints')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_path / "scfa_timepoint_coverage.png", dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}/scfa_timepoint_coverage.png")
    
    plt.close('all')


def main():
    parser = argparse.ArgumentParser(
        description="Prepare LLM biomarker data for BIRDMAn analysis"
    )
    parser.add_argument(
        "--biomarkers",
        required=True,
        help="Path to LLM_Biomarkers_Long.csv"
    )
    parser.add_argument(
        "--output-dir",
        default="./prepared_data",
        help="Output directory"
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.6,
        help="Minimum biomarker coverage (0-1) for timepoint selection"
    )
    parser.add_argument(
        "--filter-good-range",
        action="store_true",
        help="Filter values by good_range column"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Data Preparation for BIRDMAn Analysis")
    print("=" * 70)
    
    # Load data
    print("\n1. Loading biomarker data...")
    df = load_biomarker_data(args.biomarkers)
    
    # Filter by good range if requested
    if args.filter_good_range:
        print("\n2. Filtering by good range...")
        df = filter_good_range(df)
    
    # Calculate coverage
    print("\n3. Calculating biomarker coverage...")
    coverage = calculate_biomarker_coverage(df)
    
    # Identify good SCFA timepoints
    print("\n4. Identifying timepoints with good SCFA coverage...")
    good_timepoints, tp_coverage = identify_scfa_timepoints(
        df, 
        min_coverage=args.min_coverage
    )
    
    # Harmonize data
    if good_timepoints:
        print("\n5. Harmonizing biomarker data...")
        harmonized = harmonize_biomarkers(df, good_timepoints)
        
        # Calculate total SCFA
        harmonized = calculate_total_scfa(harmonized)
        
        # Save harmonized data
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        harmonized_file = output_path / "harmonized_biomarkers.csv"
        harmonized.to_csv(harmonized_file, index=False)
        print(f"\n  Saved: {harmonized_file}")
    else:
        print("\n⚠️  No timepoints with sufficient SCFA coverage found")
        harmonized = None
    
    # Create visualizations
    print("\n6. Creating visualizations...")
    visualize_data_coverage(coverage, tp_coverage, args.output_dir)
    
    # Save summary
    summary = {
        "total_biomarkers": len(coverage),
        "good_range_filter": args.filter_good_range,
        "min_coverage_threshold": args.min_coverage,
        "good_timepoints": good_timepoints if good_timepoints else [],
        "n_good_timepoints": len(good_timepoints) if good_timepoints else 0
    }
    
    import json
    with open(output_path / "preparation_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "=" * 70)
    print("Data preparation complete!")
    print("=" * 70)
    print(f"\nOutput directory: {args.output_dir}")
    print("\nNext steps:")
    print("1. Match harmonized_biomarkers.csv with microbiome samples")
    print("2. Run convert_to_biom.py to create BIOM tables")
    print("3. Run BIRDMAn analysis")


if __name__ == "__main__":
    main()
