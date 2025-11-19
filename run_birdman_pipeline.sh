#!/bin/bash
# BIRDMAn Analysis Workflow for SCFA-Microbiome Correlation Study
# 
# Hypothesis: Patterns in gut ecology correlate with fluctuations in SCFA 
# abundance across the dataset, irrespective of time, age, body weight, or cholesterol.
#
# This script automates the complete BIRDMAn analysis pipeline

set -e  # Exit on error

# Configuration
THREADS=32
OUTPUT_BASE="birdman_analysis"
EFFECT_SIZE_THRESHOLD=2

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_step() {
    echo -e "\n${GREEN}===================================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${GREEN}===================================================${NC}\n"
}

print_error() {
    echo -e "${RED}ERROR: $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

# Check if required files exist
check_required_files() {
    print_step "Step 1: Checking required files"
    
    if [ ! -f "$MICROBIOME_DATA" ]; then
        print_error "Microbiome data file not found: $MICROBIOME_DATA"
        exit 1
    fi
    
    if [ ! -f "$METABOLOMICS_DATA" ]; then
        print_error "Metabolomics data file not found: $METABOLOMICS_DATA"
        exit 1
    fi
    
    echo "✓ All required files found"
}

# Convert Excel data to BIOM format
convert_data() {
    print_step "Step 2: Converting Excel data to BIOM format"
    
    python convert_to_biom.py \
        --microbiome "$MICROBIOME_DATA" \
        --metabolomics "$METABOLOMICS_DATA" \
        ${COVARIATES_DATA:+--covariates "$COVARIATES_DATA"} \
        --output-dir "$OUTPUT_BASE/biom_data" \
        --depth 10000 \
        --sample-col sample_id
    
    if [ $? -eq 0 ]; then
        echo "✓ Data conversion completed successfully"
    else
        print_error "Data conversion failed"
        exit 1
    fi
}

# Import BIOM tables into QIIME2
import_qiime2() {
    print_step "Step 3: Importing data into QIIME2 format"
    
    # Import microbiome feature table
    qiime tools import \
        --input-path "$OUTPUT_BASE/biom_data/microbiome_table.biom" \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path "$OUTPUT_BASE/feature-table.qza"
    
    echo "✓ Microbiome feature table imported"
    
    # Import metabolomics feature table (optional, for separate analysis)
    qiime tools import \
        --input-path "$OUTPUT_BASE/biom_data/metabolomics_table.biom" \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path "$OUTPUT_BASE/metabolomics-table.qza"
    
    echo "✓ Metabolomics feature table imported"
}

# Run BIRDMAn differential abundance analysis
run_birdman() {
    print_step "Step 4: Running BIRDMAn differential abundance analysis"
    
    # Model: SCFA level as primary variable of interest + covariates
    # Formula: scfa_level_normalized + time + body_weight + cholesterol
    # This will identify which microbial features are associated with SCFA levels
    
    qiime birdman run \
        --i-table "$OUTPUT_BASE/feature-table.qza" \
        --m-metadata-file "$OUTPUT_BASE/biom_data/metadata.tsv" \
        --p-formula "scfa_level_normalized+time+body_weight+cholesterol" \
        --o-output-dir "$OUTPUT_BASE/birdman_results.qza" \
        --p-threads $THREADS \
        --verbose
    
    if [ $? -eq 0 ]; then
        echo "✓ BIRDMAn analysis completed"
    else
        print_error "BIRDMAn analysis failed"
        exit 1
    fi
}

# Run alternative models for comparison
run_alternative_models() {
    print_step "Step 5: Running alternative models for comparison"
    
    # Model 1: SCFA only (main effect)
    print_warning "Running Model 1: SCFA level only"
    qiime birdman run \
        --i-table "$OUTPUT_BASE/feature-table.qza" \
        --m-metadata-file "$OUTPUT_BASE/biom_data/metadata.tsv" \
        --p-formula "scfa_level_normalized" \
        --o-output-dir "$OUTPUT_BASE/birdman_results_scfa_only.qza" \
        --p-threads $THREADS \
        --verbose
    
    # Model 2: SCFA + time interaction (if you want to explore time effects)
    print_warning "Running Model 2: SCFA + time"
    qiime birdman run \
        --i-table "$OUTPUT_BASE/feature-table.qza" \
        --m-metadata-file "$OUTPUT_BASE/biom_data/metadata.tsv" \
        --p-formula "scfa_level_normalized+time" \
        --o-output-dir "$OUTPUT_BASE/birdman_results_scfa_time.qza" \
        --p-threads $THREADS \
        --verbose
    
    echo "✓ Alternative models completed"
}

# Visualize BIRDMAn results
visualize_results() {
    print_step "Step 6: Creating visualizations"
    
    # Main model visualization
    qiime birdman plot \
        --i-data "$OUTPUT_BASE/birdman_results.qza" \
        --i-table "$OUTPUT_BASE/feature-table.qza" \
        --m-metadata-file "$OUTPUT_BASE/biom_data/metadata.tsv" \
        --o-visualization "$OUTPUT_BASE/birdman_plot.qzv" \
        --p-palette rainbow \
        --p-chart-style forest \
        --p-effect-size-threshold $EFFECT_SIZE_THRESHOLD \
        --verbose
    
    echo "✓ Main visualization created: $OUTPUT_BASE/birdman_plot.qzv"
    
    # SCFA-only model visualization
    qiime birdman plot \
        --i-data "$OUTPUT_BASE/birdman_results_scfa_only.qza" \
        --i-table "$OUTPUT_BASE/feature-table.qza" \
        --m-metadata-file "$OUTPUT_BASE/biom_data/metadata.tsv" \
        --o-visualization "$OUTPUT_BASE/birdman_plot_scfa_only.qzv" \
        --p-palette rainbow \
        --p-chart-style forest \
        --p-effect-size-threshold $EFFECT_SIZE_THRESHOLD \
        --verbose
    
    echo "✓ SCFA-only visualization created"
}

# Export results for downstream analysis
export_results() {
    print_step "Step 7: Exporting results for downstream analysis"
    
    # Export differentially abundant features
    qiime tools export \
        --input-path "$OUTPUT_BASE/birdman_results.qza" \
        --output-path "$OUTPUT_BASE/exported_results"
    
    echo "✓ Results exported to: $OUTPUT_BASE/exported_results"
    
    # Create summary report
    cat > "$OUTPUT_BASE/analysis_summary.txt" << EOF
BIRDMAn Analysis Summary
========================

Hypothesis:
Patterns in gut ecology correlate with fluctuations in SCFA abundance 
across the dataset, irrespective of time, age, body weight, or cholesterol.

Analysis Date: $(date)

Input Data:
- Microbiome: $MICROBIOME_DATA
- Metabolomics: $METABOLOMICS_DATA
- Covariates: ${COVARIATES_DATA:-"Not provided"}

Model Formula:
scfa_level_normalized + time + body_weight + cholesterol

Parameters:
- Threads: $THREADS
- Effect Size Threshold: $EFFECT_SIZE_THRESHOLD

Output Files:
- Feature table: $OUTPUT_BASE/feature-table.qza
- BIRDMAn results: $OUTPUT_BASE/birdman_results.qza
- Visualization: $OUTPUT_BASE/birdman_plot.qzv
- Exported results: $OUTPUT_BASE/exported_results/

Next Steps:
1. View visualization:
   qiime tools view $OUTPUT_BASE/birdman_plot.qzv

2. Identify top differentially abundant taxa for KEGG pathway analysis

3. Conduct over-representation analysis using KEGG to investigate if 
   differentially abundant taxa are associated with lipid sequestration pathways

EOF
    
    echo "✓ Summary report created: $OUTPUT_BASE/analysis_summary.txt"
}

# Main workflow
main() {
    echo -e "${GREEN}"
    echo "=========================================================="
    echo "  BIRDMAn Analysis Pipeline for SCFA-Microbiome Study"
    echo "=========================================================="
    echo -e "${NC}"
    
    # Parse command line arguments
    if [ "$#" -lt 2 ]; then
        echo "Usage: $0 <microbiome_excel> <metabolomics_excel> [covariates_excel]"
        echo ""
        echo "Example:"
        echo "  $0 microbiome_data.xlsx scfa_data.xlsx covariates.xlsx"
        exit 1
    fi
    
    MICROBIOME_DATA="$1"
    METABOLOMICS_DATA="$2"
    COVARIATES_DATA="${3:-}"
    
    # Create output directory
    mkdir -p "$OUTPUT_BASE"
    
    # Run pipeline steps
    check_required_files
    convert_data
    import_qiime2
    run_birdman
    run_alternative_models
    visualize_results
    export_results
    
    # Final message
    print_step "Analysis Complete!"
    echo "Results are in: $OUTPUT_BASE"
    echo ""
    echo "To view the visualization:"
    echo "  qiime tools view $OUTPUT_BASE/birdman_plot.qzv"
    echo ""
    echo "For detailed results, see:"
    echo "  $OUTPUT_BASE/analysis_summary.txt"
}

# Run main workflow
main "$@"
