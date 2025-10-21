# LargeLarryModelers
## Systems Biology Data Analysis for BENG 211

### Course Context
This repository contains exploratory data analysis (EDA) materials developed for **BENG 211: Systems Biology and Bioengineering I: Biological Components** at UC San Diego, taught by Dr. Benjamin Smarr.

**Course Description:**  
We explore Systems Biology from the molecular level, focusing on components of biological systems, their biochemical properties and functions, and the analytical tools for deciphering and mining these components. This course emphasizes understanding biological molecules in the context of dimensionality and systems-level interactions.

**Learning Objectives:**
1. Become familiar with the properties and functions of biological components
2. Understand interactions and regulations of biological components  
3. Learn techniques for measuring and perturbing biological systems
4. **Perform quantitative analysis on 'omics' data** ← *This project addresses this objective*

### Project Purpose
**The primary purpose of this EDA notebook is to provide students with a comprehensive example of how to conduct their own biological data analysis projects.** Students can use this as a template and reference to understand:
- How to approach longitudinal biological datasets
- Proper data preprocessing and cleaning techniques
- Effective visualization strategies for biological time series
- Statistical analysis methods for 'omics data
- How to interpret and communicate biological insights from data
- **Essential Python data science tools** (pandas, matplotlib, seaborn, scipy, statsmodels)
- **Data manipulation and analysis workflows** for biological research
- **Statistical modeling and machine learning applications** in systems biology
- **Reproducible data science practices** and code organization
- **Computational thinking** for biological problem-solving

This hands-on example demonstrates the complete workflow from raw data to meaningful biological conclusions, serving as a practical guide for students' independent research projects.

### Historical Dataset Context
This analysis uses a remarkable longitudinal dataset collected by **LS beginning in 1993** - representing over **30 years of continuous biological monitoring**. This dataset is one of the most comprehensive personal biological tracking datasets ever assembled, making it an invaluable resource for understanding long-term biological patterns and systems-level changes in human physiology.

The data collection began in 1993 and continues to the present day (2025), spanning multiple decades of biomarker measurements, gut microbiome analysis, and metabolic tracking. This longitudinal approach allows students to observe:
- Long-term trends in biological systems
- Seasonal and cyclical patterns in human physiology
- The evolution of measurement technologies over time
- Real-world data collection challenges and solutions

### Dataset Description

#### Original Raw Data Files (LS Collection)
**LS Biomarkers Timeseries.csv** - The original biomarker data in its raw collected format
- Contains measurements from 1993 to 2025 (30+ years of data)
- Wide format with dates as columns and biomarkers as rows
- Includes clinical measurements from multiple laboratories (YFH, UCSD, Scripps, Hood, etc.)
- Contains metadata about measurement units, normal ranges, and testing facilities

**LS Gut Microbiome.csv** - Original microbiome sequencing data
- Raw microbial abundance data from gut microbiome analysis
- Contains both 16S rRNA sequencing and metagenomic data
- Species-level abundance measurements over time

**LS_Gut_Microbiome_Summary.csv** - Summary statistics for the original microbiome data

#### Processed Data Files (For Analysis)
The table conversion notebook transforms the raw LS data into analysis-ready formats:

**LLM_Biomarkers_Wide.csv** - Machine learning-ready biomarker format
- Each row represents a time point
- Each column represents a different biomarker measurement
- Optimized for statistical analysis and modeling
- Key measurements include:
  - Weight tracking over 30+ years
  - Lipid panel (HDL, LDL, cholesterol, triglycerides)
  - Comprehensive metabolic panel (electrolytes, liver function)
  - Inflammatory markers (CRP, cytokines)
  - Metabolic markers (glucose, insulin, SCFA)
  - Vitamin and mineral levels
  - Hormone measurements

**LLM_Biomarkers_Long.csv** - Comprehensive long-format biomarker data
- Each row represents a single measurement
- Columns: Date, Biomarker, Value, Units, Normal_Range, Laboratory
- Ideal for time series analysis and visualization
- Facilitates grouping and filtering operations

**LLM_Gut_Microbiome_Wide.csv** - ML-ready microbiome format
- Rows represent time points, columns represent microbial species
- Abundance values normalized for comparative analysis
- Includes diversity metrics (Shannon, Simpson indices)

**LLM_Gut_Microbiome_Long.csv** - Comprehensive long-format microbiome data
- Each row represents abundance of one species at one time point
- Columns: Date, Species, Abundance, Sample_Type, Sequencing_Method
- Enables temporal analysis of microbial community dynamics

**LLM_Gut_Microbiome_Summary.csv** - Statistical summaries and metadata
- Summary statistics for each microbial species
- Temporal trends and stability metrics
- Diversity indices over time
- Quality control metrics

### Notebooks

#### `LLM_EDA1.ipynb` - Primary Exploratory Data Analysis
**This notebook serves as your comprehensive guide for conducting biological data analysis.** It demonstrates:

**Data Import and Preprocessing:**
- Loading multiple data formats (wide vs. long)
- Handling missing values in biological datasets
- Date/time formatting for longitudinal analysis
- Data quality assessment and validation

**Descriptive Statistics:**
- Summary statistics for 30+ years of biomarker data
- Data completeness analysis across different time periods
- Identification of measurement patterns and anomalies

**Time Series Analysis:**
- Temporal patterns in weight, SCFA, and other biomarkers
- Rolling window analysis to identify long-term trends
- Seasonal pattern detection in biological measurements
- Correlation analysis between different biomarkers over time

**Advanced Visualizations:**
- Multi-panel time series plots with dual y-axes
- Correlation heatmaps for biomarker relationships
- Distribution plots for understanding data characteristics
- Trend analysis with confidence intervals

**Biological Insights:**
- Host-microbiome interaction patterns
- Metabolic health indicators over decades
- Impact of life events on biological measurements
- Systems-level understanding of human physiology

#### `LLM_tableConversion.ipynb` - Data Format Conversion Utility
**This notebook demonstrates essential data preprocessing skills** that students will need for their own projects:

**Data Transformation Techniques:**
- Converting between wide and long data formats
- Handling complex multi-header CSV files
- Parsing date formats across different time periods
- Standardizing measurement units and laboratory values

**Data Cleaning Operations:**
- Removing blank cells and handling missing data
- Standardizing biomarker names and categories
- Creating consistent categorical variables (e.g., laboratory sources)
- Quality control checks and data validation

**Feature Engineering:**
- Creating biomarker panels (Lipid, Metabolic, etc.)
- Generating summary statistics and derived metrics
- Building lookup dictionaries for data mapping
- Calculating temporal features and rolling statistics

**Real-world Data Challenges:**
- Handling inconsistent data collection protocols
- Managing data from multiple laboratories with different standards
- Dealing with evolving measurement technologies over 30 years
- Creating robust pipelines for ongoing data collection

### Getting Started

#### Prerequisites
- Python 3.7+
- Jupyter Notebook or VS Code with Jupyter extension

#### Installation
1. Clone this repository
2. Install required packages:
```bash
pip install -r requirements.txt
```

#### Required Python Packages
- `numpy` - Numerical computing and array operations
- `pandas` - Data manipulation and analysis (essential for biomedical data)
- `matplotlib` - Plotting and visualization
- `seaborn` - Statistical data visualization
- `statsmodels` - Statistical modeling and time series analysis
- `scipy` - Scientific computing and advanced statistics

#### Running the Analysis
1. **Start with the Table Conversion notebook** (`LLM_tableConversion.ipynb`) to understand data preprocessing
2. **Proceed to the EDA notebook** (`LLM_EDA1.ipynb`) for comprehensive analysis
3. Run cells sequentially to reproduce the analysis
4. Modify parameters and explore different aspects of the data
5. Use the patterns and techniques as templates for your own projects

### Educational Objectives
This project helps students develop practical skills in:

**Quantitative Analysis of 'Omics Data:**
- Working with real longitudinal biological datasets
- Understanding data structure and quality in biological research
- Applying statistical methods appropriate for biological measurements

**Data Science for Biology:**
- Data preprocessing and cleaning techniques
- Time series analysis for biological systems
- Correlation analysis and pattern recognition
- Data visualization for biological insights

**Systems Biology Thinking:**
- Connecting molecular-level measurements to system-level understanding
- Understanding temporal dynamics in biological systems
- Recognizing patterns and relationships in complex biological data

### Python Data Science Tools & Skills
A core component of this project is learning to use essential **Python data science tools** that are fundamental to modern biological research and systems biology:

**pandas - Data Manipulation & Analysis:**
- Loading and parsing complex biological datasets (CSV, Excel, databases)
- Data cleaning, filtering, and transformation operations
- Handling missing values and data quality issues
- Merging and joining datasets from multiple sources
- Grouping and aggregating biological measurements
- Time series data manipulation and resampling

**matplotlib & seaborn - Data Visualization:**
- Creating publication-quality plots for biological data
- Time series visualizations for longitudinal studies
- Correlation heatmaps and scatter plot matrices
- Multi-panel figures with subplots
- Custom color schemes and styling for biological presentations
- Interactive plotting for data exploration

**scipy & statsmodels - Statistical Analysis:**
- Descriptive statistics for biological measurements
- Correlation and regression analysis
- Time series decomposition and trend analysis
- Statistical testing (t-tests, ANOVA, non-parametric tests)
- Confidence intervals and error estimation
- Hypothesis testing for biological research

**numpy - Numerical Computing:**
- Array operations for efficient data processing
- Mathematical functions for biological calculations
- Linear algebra operations for multivariate analysis
- Random sampling and simulation methods

**scikit-learn - Machine Learning (Advanced Applications):**
- Dimensionality reduction techniques (PCA, t-SNE)
- Clustering methods for biological pattern discovery
- Classification and regression models for biomarker prediction
- Cross-validation and model evaluation
- Feature selection and importance analysis

**Integration Skills:**
- Building complete data analysis pipelines
- Reproducible research workflows
- Version control for data science projects
- Documentation and code organization
- Collaborative data science practices

These tools form the **standard toolkit for computational biology and bioinformatics**, preparing students for advanced coursework and research in systems biology, bioengineering, and data-driven biological discovery.

### Learning Outcomes for Student Projects
After studying this analysis, students will be able to:

**Technical Skills:**
- Load and preprocess complex biological datasets
- Convert between different data formats (wide/long)
- Handle missing data and quality control issues
- Create effective visualizations for biological time series
- Perform statistical analysis on longitudinal biological data
- Apply correlation analysis to identify biological relationships

**Analytical Thinking:**
- Develop hypotheses about biological systems from data
- Interpret statistical results in biological context
- Recognize patterns and anomalies in biological measurements
- Connect molecular measurements to physiological outcomes

**Project Implementation:**
- Design data collection strategies for biological research
- Create reproducible analysis workflows
- Document and communicate biological insights effectively

### For Student Projects
**Use this repository as your starting template for biological data analysis projects.** The workflows demonstrated here can be adapted for:
- Analysis of your own experimental data
- Public biological datasets (GEO, SRA, etc.)
- Clinical trial data analysis
- Environmental biology datasets
- Comparative biology studies

### Data Scientific Significance
This 30+ year longitudinal dataset represents:
- One of the longest continuous personal biological monitoring studies
- Integration of traditional clinical biomarkers with modern 'omics technologies
- Real-world example of precision medicine data collection
- Demonstration of how measurement technologies have evolved
- Proof-of-concept for personalized biological monitoring

### Contact
- **Instructor**: Dr. Benjamin Smarr (bsmarr@ucsd.edu)
- **Me**: Conan Minihan (jrminihan@ucsd.edu)

### License
This educational material is provided for academic use in BENG 211.
