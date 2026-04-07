# R scripts for LASSO-Cox and ESTIMATE analyses in HBV-HCC

This repository contains the R scripts used for LASSO-Cox regression analysis and ESTIMATE score analysis in our study of HBV-related hepatocellular carcinoma (HBV-HCC).

## Files included

- `LASSO_Cox_analysis.R`: R script for LASSO-Cox model construction, coefficient extraction, risk score calculation, and LASSO figure generation
- `ESTIMATE_score_analysis.R`: R script for ESTIMATE score integration, risk group stratification, violin plot visualization, and Wilcoxon test analysis

## Script 1: LASSO-Cox analysis

### Main functions
- load and preprocess survival and gene expression data
- perform LASSO-Cox regression with 10-fold cross-validation
- extract selected genes at `lambda.min` and `lambda.1se`
- identify a 4-gene model nearest to log(lambda) = -2.173
- calculate individual risk scores
- export coefficients and model results
- generate LASSO coefficient and cross-validation plots

### Required input file
- `TCGA_LIHC_LASSO_input_13genes_survival.csv`

### Main output files
- `LASSO_lambda_nonzero_table.csv`
- `LASSO_selected_genes_lambda_min.csv`
- `LASSO_selected_genes_lambda_1se.csv`
- `LASSO_selected_genes_target_4gene_nearest_log_minus_2.173.csv`
- `RiskScore_target_4gene.csv`
- `RiskScore_formula_target_4gene.txt`
- `LASSO_Cox_C_D_style.pdf`
- `LASSO_Cox_C_D_style.png`

## Script 2: ESTIMATE score analysis

### Main functions
- read the 13-gene tumor expression matrix
- calculate risk scores using EHMT2 and INMT coefficients
- classify samples into high-risk and low-risk groups
- merge risk group data with ESTIMATE scores
- generate violin plots for StromalScore, ImmuneScore, and ESTIMATEScore
- perform Wilcoxon rank-sum tests and export p-values

### Required input files
- `TCGA_LIHC_13gene_tumor_expression.csv`
- `estimate_scores.txt`

### Main output files
- `TCGA_LIHC_LASSO_riskScore_from_13gene.csv`
- `TCGA_LIHC_LASSO_ESTIMATE_merged.csv`
- `TCGA_LIHC_ESTIMATE_violin_3panel.pdf`
- `TCGA_LIHC_ESTIMATE_violin_3panel.png`
- `TCGA_LIHC_ESTIMATE_wilcox_pvalues.csv`

## Required R packages

### For LASSO-Cox analysis
- glmnet
- survival

### For ESTIMATE score analysis
- ggplot2
- ggpubr
- dplyr
- patchwork

## Notes
Please ensure that all input files are placed in the working directory before running the scripts.

A fixed random seed was used in the LASSO-Cox script to improve reproducibility.

## Citation
If you use these scripts, please cite the associated study.
