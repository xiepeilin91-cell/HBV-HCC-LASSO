# =========================================================
# LASSO-Cox script with figure style similar to the user's panel C/D
# Input file: TCGA_LIHC_LASSO_input_13genes_survival.csv
# =========================================================
needed_pkgs <- c("glmnet", "survival")
for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(glmnet)
library(survival)

input_file <- "TCGA_LIHC_LASSO_input_13genes_survival.csv"
if (!file.exists(input_file)) stop(paste0("Input file not found: ", input_file))

dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

gene_list <- c("EHMT2","DNMT3B","EZH2","PRDM9","INMT","DNMT1",
               "PHF19","BHMT","NSD2","GNMT","CSKMT","TARBP1","ZNF711")
required_cols <- c("OS.time", "OS.status", gene_list)
missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))

keep <- required_cols
if ("patient_id" %in% colnames(dat)) keep <- c("patient_id", keep)
if ("SampleID" %in% colnames(dat)) keep <- c(keep, "SampleID")
dat <- dat[, unique(keep)]
dat <- dat[complete.cases(dat[, required_cols]), ]

dat$OS.time <- as.numeric(dat$OS.time)
dat$OS.status <- as.numeric(dat$OS.status)
for (g in gene_list) dat[[g]] <- as.numeric(dat[[g]])
if (!all(dat$OS.status %in% c(0, 1))) stop("OS.status must be coded as 0/1.")

# Important: log2 transform before LASSO
dat[, gene_list] <- log2(dat[, gene_list] + 1)

x <- as.matrix(dat[, gene_list])
storage.mode(x) <- "double"
y <- Surv(dat$OS.time, dat$OS.status)

set.seed(1234)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10,
                   type.measure = "deviance", standardize = TRUE, maxit = 100000)
fit <- glmnet(x, y, family = "cox", alpha = 1,
              standardize = TRUE, maxit = 100000)

extract_nonzero_coef <- function(model, s_value) {
  beta <- as.matrix(coef(model, s = s_value))
  res <- data.frame(Gene = rownames(beta), Coefficient = as.numeric(beta[,1]), stringsAsFactors = FALSE)
  res <- res[res$Coefficient != 0, , drop = FALSE]
  rownames(res) <- NULL
  res
}

beta_all <- as.matrix(coef(fit, s = fit$lambda))
nonzero_counts <- colSums(beta_all != 0)
lambda_table <- data.frame(lambda = fit$lambda, log_lambda = log(fit$lambda), nonzero = nonzero_counts)
write.csv(lambda_table, "LASSO_lambda_nonzero_table.csv", row.names = FALSE)

res_lambda_min <- extract_nonzero_coef(cvfit, cvfit$lambda.min)
res_lambda_1se <- extract_nonzero_coef(cvfit, cvfit$lambda.1se)
write.csv(res_lambda_min, "LASSO_selected_genes_lambda_min.csv", row.names = FALSE)
write.csv(res_lambda_1se, "LASSO_selected_genes_lambda_1se.csv", row.names = FALSE)

cand4 <- lambda_table[lambda_table$nonzero == 4, , drop = FALSE]
if (nrow(cand4) == 0) stop("No lambda with exactly 4 non-zero genes was found.")
target_log_lambda <- -2.173
idx4 <- which.min(abs(cand4$log_lambda - target_log_lambda))
target_lambda_exact <- cand4$lambda[idx4]
target_log_lambda_exact <- cand4$log_lambda[idx4]
res_target_4gene <- extract_nonzero_coef(fit, target_lambda_exact)
write.csv(res_target_4gene, "LASSO_selected_genes_target_4gene_nearest_log_minus_2.173.csv", row.names = FALSE)

selected_genes <- res_target_4gene$Gene
selected_coef <- res_target_4gene$Coefficient
risk_score <- as.numeric(as.matrix(dat[, selected_genes, drop = FALSE]) %*% selected_coef)
risk_df <- data.frame(
  patient_id = if ("patient_id" %in% colnames(dat)) dat$patient_id else NA,
  SampleID = if ("SampleID" %in% colnames(dat)) dat$SampleID else NA,
  OS.time = dat$OS.time,
  OS.status = dat$OS.status,
  RiskScore = risk_score,
  stringsAsFactors = FALSE
)
write.csv(risk_df, "RiskScore_target_4gene.csv", row.names = FALSE)
formula_text <- paste0("Risk score = ", paste(sprintf("%.6f × %s", selected_coef, selected_genes), collapse = " + "))
formula_text <- gsub("\\+ -", "- ", formula_text)
writeLines(c(paste0("Chosen lambda = ", signif(target_lambda_exact, 8)),
             paste0("Chosen log(lambda) = ", signif(target_log_lambda_exact, 8)),
             formula_text), con = "RiskScore_formula_target_4gene.txt")

make_lasso_plot <- function(device = c("pdf", "png"), file) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(file, width = 11, height = 5.2, family = "serif")
  } else {
    png(file, width = 2200, height = 1040, res = 220)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() })
  par(mfrow = c(1, 2), mar = c(5, 5, 2.3, 1.5), mgp = c(2.6, 0.8, 0),
      cex.lab = 1.05, cex.axis = 0.95, font.lab = 1, lwd = 1.2, bty = "o")
  plot(fit, xvar = "lambda", label = TRUE, main = "", xlab = "Log Lambda", ylab = "Coefficients", lwd = 1)
  mtext("C", side = 3, line = 1.0, adj = -0.16, cex = 1.5, font = 2)
  plot(cvfit, main = "", xlab = expression(Log(lambda)), ylab = "Partial Likelihood Deviance")
  abline(v = log(cvfit$lambda.min), lty = 2, lwd = 1)
  abline(v = log(cvfit$lambda.1se), lty = 2, lwd = 1)
  mtext("D", side = 3, line = 1.0, adj = -0.16, cex = 1.5, font = 2)
}

make_lasso_plot("pdf", "LASSO_Cox_C_D_style.pdf")
make_lasso_plot("png", "LASSO_Cox_C_D_style.png")

cat("\n===== LASSO-Cox completed =====\n")
cat("Samples used:", nrow(dat), "\n")
cat("lambda.min:", cvfit$lambda.min, "\n")
cat("lambda.1se:", cvfit$lambda.1se, "\n")
cat("Chosen exact lambda for 4-gene model:", target_lambda_exact, "\n")
cat("Chosen exact log(lambda):", target_log_lambda_exact, "\n\n")
print(res_target_4gene)
cat(formula_text, "\n")
