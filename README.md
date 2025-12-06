# Exploratory analysis of relative microbial abundance along a SCFA gradient in Patient LS

This repository includes code for the analysis of SCFAs & amplicon data in the Patient LS dataset.

**Course:** BENG 211 - Fall 2025  
**Authors:** Cale [@caleoseymour](https://github.com/caleoseymour), Leo [@l1joseph](https://github.com/l1joseph), Ziheng [@ziw165](https://github.com/ziw165), Yashwin [@YASHWINLM](https://github.com/YASHWINLM)

## Hypothesis and aims

Patient LS is an adult male who has undertaken continuous self-monitoring from 1993 until 2025. Data collected include both gut SCFA and microbiota quantification. In this project, we explore the relationship between gut microbiota and SCFA abundance using data from Patient LS. We hypothesized that patterns in gut ecology correlate with fluctuations in SCFA abundance irrespective of age, body weight, or cholesterol. To test this hypothesis, we implement a series of linear models via Maaslin3, including covariates representing the calendar year of sampling and the weight of Patient LS at the time of specimen collection.

## Methods

Data from Patient LS were obtained from the repository [ConanMinihan/LargeLarryModelers](https://github.com/ConanMinihan/LargeLarryModelers).

1. Data were harmonized such that 16S rRNA amplicon were matched with SCFA collection. 102 samples were paired.
2. Covariates were selected using PERMANOVA. Year and weight both correlate to microbial community composition.
3. Maaslin3 was used to implement linear models testing each SCFA vs. the abundance of each microbial OTU, accounting for year and weight.
4. Select representative genomes were probed for genes potentially involved in SCFA metabolism. 

## Contents

| Path                        | FileType  | Description                                 |
|-----------------------------|-----------|---------------------------------------------|
| `Data/`                     | Directory | Project source data                         |
| `KEGG_analysis/`            | Directory |Code & data for functional analysis.         |
| `beta-diversity/`           | Directory |Figures for beta-diversity analysis.         |
| `maaslin3/`                 | Directory |Results from maaslin3 diffabund analysis.    |
| `reports/`                  | Directory |knitr report files for each Rmarkdown file.  |
| `resources/`                | Directory |Taxonomic ref data for data cleansing.       |
| `sanitized-data/`           | Directory |Directory containing cleansed datasets.      |
| `scripts/`                  | Directory |Scripts for miscellanous ecology functions.  |
| `main.R`                    | R script  |'source' this to run all analyses.           |
| `betadiversity.Rmd`         | Rmarkdown |Code used for covariate analysis.            |
| `differential-abundance.Rmd`| Rmarkdown |Code used for diffabund analysis.            |
| `diversity.Rmd`             | Rmarkdown |Preliminary analysis of overall ecology.     |
| `import-and-preprocess.Rmd` | Rmarkdown |Code used for data cleansing/harmonization.  |
| `figure1.png`               | PNG       |Resulting figure.                            |

---

*Last updated: 6 December 2025*
