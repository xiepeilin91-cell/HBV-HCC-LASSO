rm(list = ls())
options(stringsAsFactors = FALSE)

############################
# 0. Load packages
############################
pkg_needed <- c("ggplot2", "ggpubr", "dplyr", "patchwork")
pkg_new <- pkg_needed[!(pkg_needed %in% installed.packages()[, "Package"])]
if(length(pkg_new) > 0){
  install.packages(pkg_new)
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

############################
# 1. File settings
############################
expr_file <- "TCGA_LIHC_13gene_tumor_expression.csv"
estimate_file <- "estimate_scores.txt"

# LASSO coefficients
coef_EHMT2 <- 1.189
coef_INMT  <- -0.282

############################
# 2. Read expression matrix
############################
expr_raw <- read.csv(expr_file, header = TRUE, check.names = FALSE)

cat("Expression matrix dim:\n")
print(dim(expr_raw))
cat("Expression matrix columns:\n")
print(colnames(expr_raw))

############################
# 3. Check required columns
############################
required_expr_cols <- c("SampleID", "patient_id", "EHMT2", "INMT")
missing_expr_cols <- setdiff(required_expr_cols, colnames(expr_raw))

if(length(missing_expr_cols) > 0){
  stop(paste(
    "表达矩阵缺少这些列：",
    paste(missing_expr_cols, collapse = ", "),
    "\n请先确认你的13基因表达矩阵里是否真的包含 INMT。"
  ))
}

############################
# 4. Extract EHMT2 and INMT
############################
expr_data <- expr_raw %>%
  dplyr::select(SampleID, patient_id, EHMT2, INMT)

colnames(expr_data) <- c("ID", "PatientID", "EHMT2", "INMT")

expr_data$ID <- as.character(expr_data$ID)
expr_data$PatientID <- as.character(expr_data$PatientID)
expr_data$EHMT2 <- as.numeric(expr_data$EHMT2)
expr_data$INMT  <- as.numeric(expr_data$INMT)

expr_data <- expr_data %>%
  distinct(PatientID, .keep_all = TRUE)

cat("expr_data dim:\n")
print(dim(expr_data))
print(head(expr_data))

############################
# 5. Calculate risk score
############################
expr_data <- expr_data %>%
  mutate(
    riskScore = coef_EHMT2 * EHMT2 + coef_INMT * INMT
  )

risk_cutoff <- median(expr_data$riskScore, na.rm = TRUE)

expr_data <- expr_data %>%
  mutate(
    risk = ifelse(riskScore >= risk_cutoff, "high_risk", "low_risk")
  )

expr_data$risk <- factor(expr_data$risk, levels = c("high_risk", "low_risk"))

cat("Risk group counts:\n")
print(table(expr_data$risk))

write.csv(expr_data, "TCGA_LIHC_LASSO_riskScore_from_13gene.csv", row.names = FALSE)

############################
# 6. Read ESTIMATE table
############################
estimate <- read.delim(
  estimate_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(estimate)[1] <- "ID"

required_estimate_cols <- c("ID", "StromalScore", "ImmuneScore", "ESTIMATEScore")
missing_estimate_cols <- setdiff(required_estimate_cols, colnames(estimate))

if(length(missing_estimate_cols) > 0){
  stop(paste("ESTIMATE表缺少这些列：", paste(missing_estimate_cols, collapse = ", ")))
}

estimate$ID <- as.character(estimate$ID)
estimate$PatientID <- substr(estimate$ID, 1, 12)

estimate <- estimate %>%
  dplyr::select(PatientID, ID, StromalScore, ImmuneScore, ESTIMATEScore) %>%
  distinct(PatientID, .keep_all = TRUE)

cat("estimate dim:\n")
print(dim(estimate))
print(head(estimate))

############################
# 7. Merge risk data with ESTIMATE
############################
plot_data <- inner_join(expr_data, estimate, by = "PatientID")

cat("plot_data dim:\n")
print(dim(plot_data))
print(head(plot_data))

if(nrow(plot_data) == 0){
  stop("合并后 plot_data 为空，请检查 patient_id 与 ESTIMATE 表前12位ID 是否一致。")
}

write.csv(plot_data, "TCGA_LIHC_LASSO_ESTIMATE_merged.csv", row.names = FALSE)

############################
# 8. Plot settings
############################
my_colors <- c(
  "high_risk" = "#0B78C5",
  "low_risk"  = "#F2C300"
)

make_violin_plot <- function(data, score_col, title_text){
  ggplot(data, aes(x = risk, y = .data[[score_col]], fill = risk)) +
    geom_violin(trim = FALSE, color = NA, alpha = 1) +
    geom_boxplot(
      width = 0.25,
      fill = "white",
      color = "black",
      outlier.shape = 16,
      outlier.size = 1.5
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      size = 5
    ) +
    scale_fill_manual(values = my_colors) +
    labs(
      title = title_text,
      x = "risk",
      y = NULL,
      fill = "risk"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0, face = "plain", size = 18),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 13),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13),
      panel.grid.major = element_line(color = "#D9D9D9"),
      panel.grid.minor = element_line(color = "#EEEEEE")
    )
}

############################
# 9. Draw plots
############################
p1 <- make_violin_plot(plot_data, "StromalScore", "StromalScore")
p2 <- make_violin_plot(plot_data, "ImmuneScore", "ImmuneScore")
p3 <- make_violin_plot(plot_data, "ESTIMATEScore", "ESTIMATEScore")

final_plot <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)

ggsave("TCGA_LIHC_ESTIMATE_violin_3panel.pdf", final_plot, width = 18, height = 7)
ggsave("TCGA_LIHC_ESTIMATE_violin_3panel.png", final_plot, width = 18, height = 7, dpi = 300)

print(final_plot)

############################
# 10. Export p-values
############################
p_stromal  <- wilcox.test(StromalScore ~ risk, data = plot_data)$p.value
p_immune   <- wilcox.test(ImmuneScore ~ risk, data = plot_data)$p.value
p_estimate <- wilcox.test(ESTIMATEScore ~ risk, data = plot_data)$p.value

p_table <- data.frame(
  Score = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
  P_value = c(p_stromal, p_immune, p_estimate)
)

write.csv(p_table, "TCGA_LIHC_ESTIMATE_wilcox_pvalues.csv", row.names = FALSE)

cat("All done.\n")
