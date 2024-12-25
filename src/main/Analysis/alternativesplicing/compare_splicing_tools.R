library("data.table")
library(ROCR)
library(ggplot2)
load("data/alternativesplicing/diff_psi_test.RData")
psi.files <- paste0("data/alternativesplicing/psi/sample", 1:10, ".psi") #list of files


source("src/main/Analysis/alternativesplicing/LRS.R")
psi_data <- diff.splicing(psi.files, group)

source("src/main/Analysis/alternativesplicing/dexseq_script.R") # gives us dxr1

load("data/alternativesplicing/differential_exons.RData") # Gives us differential.skipped

dxr1$genomicData

dxr1_dt <- as.data.table(dxr1)
psi_data_dt <- as.data.table(psi_data)
skipped_dt <- as.data.table(differential.skipped)

# Extract the ranges from dxr1's genomicData for easier comparison
dxr1_dt[, exon_range := paste0(start(ranges(dxr1$genomicData)), "-", end(ranges(dxr1$genomicData)) + 1)]
psi_data_dt[, exon_range := exon]  # Assuming `psi_data` has exon in "start-end" format

# Perform the merge by `groupID` and `exon_range`
merged_dt <- merge(
  psi_data_dt[, .(gene, exon_range, p_lrs = padj)],  # Only necessary columns
  dxr1_dt[, .(groupID, exon_range, pd_dseq = padj)],  # Only necessary columns
  by.x = c("gene", "exon_range"),  # Matching by gene and exon_range
  by.y = c("groupID", "exon_range"),  # Matching by groupID and exon_range

  all = FALSE
)

# Determine ground truth flags based on `differential.skipped`
# Join to check if genomicData overlaps with skipped ranges
skipped_dt[, skip_range := paste0(start, "-", end + 1)]
merged_dt[, ground_truth := exon_range %in% skipped_dt$skip_range]

# na ehandling
# merged_dt$pd_dseq[is.na(merged_dt$pd_dseq)] <- 1
merged_dt <- merged_dt[!is.na(merged_dt$pd_dseq),]
pred_lrs <- prediction(1 - merged_dt$p_lrs, merged_dt$ground_truth)
pred_dseq <- prediction(1 - merged_dt$pd_dseq, merged_dt$ground_truth)

# Create performance objects for ROC curves
perf_lrs <- performance(pred_lrs, "tpr", "fpr")
perf_dseq <- performance(pred_dseq, "tpr", "fpr")

# Calculate AUC for pd_dseq
auc_dseq <- performance(pred_dseq, measure = "auc")
auc_dseq_value <- auc_dseq@y.values[[1]]

# Calculate AUC for p_lrs
auc_lrs <- performance(pred_lrs, measure = "auc")
auc_lrs_value <- auc_lrs@y.values[[1]]

# Plot ROC curve for pd_dseq
par(pty = "s")  # Force square plotting region
png("src/main/Analysis/alternativesplicing/roc-comparison.png", width = 6, height = 6, units = "in", res = 300)

roc_dseq <- performance(pred_dseq, measure = "tpr", x.measure = "fpr")
plot(roc_dseq, col = "blue", main = "ROC Curves", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i', asp = 1)

# Add ROC curve for p_lrs
roc_lrs <- performance(pred_lrs, measure = "tpr", x.measure = "fpr")
plot(roc_lrs, col = "red", lwd = 2, add = TRUE,
     xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
# Add the y = x line (representing random classifier)
abline(a = 0, b = 1, col = "gray", lty = 2)  # Gray dashed line

# Add grid lines for better visibility
grid(lwd = 1, col = "gray", lty = 1)  # Add major grid lines
grid(lwd = 0.5, col = "lightgray", lty = 3, minor = TRUE)  # Add minor grid lines

# Add legend with AUC values
legend("bottomright",
       legend = c(paste("pd_dseq (AUC =", round(auc_dseq_value, 3), ")"),
                  paste("p_lrs (AUC =", round(auc_lrs_value, 3), ")")),
       col = c("blue", "red"), lwd = 2)
dev.off()
# Filter for adjusted p-value < 0.05 for pd_dseq
filtered_dseq <- merged_dt[pd_dseq < 0.05]

# Calculate sensitivity and FDR for pd_dseq
sensitivity_dseq <- sum(filtered_dseq$ground_truth) / nrow(filtered_dseq)
fdr_dseq <- sum(!filtered_dseq$ground_truth) / nrow(filtered_dseq)

# Filter for adjusted p-value < 0.05 for p_lrs
filtered_lrs <- merged_dt[p_lrs < 0.05]

# Calculate sensitivity and FDR for p_lrs
sensitivity_lrs <- sum(filtered_lrs$ground_truth) / nrow(filtered_lrs)
fdr_lrs <- sum(!filtered_lrs$ground_truth) / nrow(filtered_lrs)

# Print sensitivity and FDR for both
cat("Sensitivity for pd_dseq:", sensitivity_dseq, "\n")
cat("FDR for pd_dseq:", fdr_dseq, "\n")
cat("Sensitivity for p_lrs:", sensitivity_lrs, "\n")
cat("FDR for p_lrs:", fdr_lrs, "\n")