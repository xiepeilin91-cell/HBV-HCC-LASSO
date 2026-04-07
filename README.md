# LASSO-Cox analysis script for HBV-HCC

This repository contains the R script used for LASSO-Cox regression analysis in our study of HBV-related hepatocellular carcinoma (HBV-HCC).

## File included

- `LASSO_Cox_analysis.R`: main R script for LASSO-Cox model construction, risk score calculation, and figure generation

## Required R packages

The script uses the following R packages:

- glmnet
- survival

## Input file

The script requires the following input file in the working directory:

- `TCGA_LIHC_LASSO_input_13genes_survival.csv`

## Main analysis workflow

This script performs the following steps:

1. Loads the input survival and gene expression data
2. Checks required columns
3. Applies log2 transformation to gene expression values
4. Performs LASSO-Cox regression with 10-fold cross-validation
5. Extracts selected genes at `lambda.min` and `lambda.1se`
6. Identifies a 4-gene model nearest to log(lambda) = -2.173
7. Calculates individual risk scores
8. Exports coefficient tables and risk score results
9. Generates LASSO coefficient and cross-validation plots

## Output files

Running the script will generate the following files:

- `LASSO_lambda_nonzero_table.csv`
- `LASSO_selected_genes_lambda_min.csv`
- `LASSO_selected_genes_lambda_1se.csv`
- `LASSO_selected_genes_target_4gene_nearest_log_minus_2.173.csv`
- `RiskScore_target_4gene.csv`
- `RiskScore_formula_target_4gene.txt`
- `LASSO_Cox_C_D_style.pdf`
- `LASSO_Cox_C_D_style.png`

## Notes

Please ensure that the input file is placed in the working directory before running the script.

A fixed random seed (`set.seed(1234)`) was used in the cross-validation process to improve reproducibility.

## Citation

If you use this script, please cite the associated study.
