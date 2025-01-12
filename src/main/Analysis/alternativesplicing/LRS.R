library("data.table")
library("dplyr")
library("ggplot2")
library("pheatmap")

LRS <- function(inclusion, total, group) { group_1_included <- 0
  group_2_included <- 0
  group_1_total <- 0
  group_2_total <- 0
  GROUP_ONE_ID <- 1
  # Initialize the variables
  for (i in seq_along(group)) { if (group[i] == GROUP_ONE_ID) { group_1_included <- group_1_included + inclusion[i]
    group_1_total <- group_1_total + total[i] } else { group_2_included <- group_2_included + inclusion[i]
    group_2_total <- group_2_total + total[i] } }
  # Model parameters basic/full
  p_0 <- (group_1_included + group_2_included) / (group_1_total + group_2_total)
  p_1 <- group_1_included / group_1_total
  p_2 <- group_2_included / group_2_total

  # Calculate the log likelihoods
  log_L_reduced <- 0
  log_L_full <- 0
  for (i in seq_along(group)) { g <- group[i]
    total_count <- total[i]
    inclusion_count <- inclusion[i]
    log_L_full <- log_L_full + dbinom(inclusion_count, total_count, ifelse(g == GROUP_ONE_ID, p_1, p_2), log = TRUE)
    log_L_reduced <- log_L_reduced + dbinom(inclusion_count, total_count, p_0, log = TRUE) }
  # calculate the likelihood ratio statistic
  lrs <- -2 * (log_L_reduced - log_L_full)
  # The LRS values follow a Chi-square distribution with df degrees of freedom  # df = #params.full.model − #params.reduced.model
  p_value <- pchisq(lrs, df = 1, lower.tail = FALSE)

  return(list(p0 = p_0, p1 = p_1, p2 = p_2, llreduced = log_L_reduced, llfull = log_L_full, lrs = lrs, pvalue = p_value)) }

process.event <- function(row, group) {
  inclusion <- numeric(length(group))
  total <- numeric(length(group))
  for (i in seq_along(group)) {
    if (is.na(row[[paste0("num_incl_reads_", i)]])) {
      inclusion[i] <- 0
    } else {
      inclusion[i] <- row[[paste0("num_incl_reads_", i)]]
    }

    if (is.na(row[[paste0("num_total_reads_", i)]])) {
      total[i] <- 0
    } else {
      total[i] <- row[[paste0("num_total_reads_", i)]]
    }
  }
  return(LRS(inclusion, total, group))
}

diff.splicing <- function(psi.files, group) { counts_dt <- fread(psi.files[1], sep = "\t", header = TRUE)
  # rename the incl excl total and psi columns and append _1 to the column names but don't rename exon and gene
  setnames(counts_dt, c("gene", "exon", paste0(names(counts_dt)[3:ncol(counts_dt)], "_1")))
  for (i in seq_along(psi.files)[-1]) { temp_data <- fread(psi.files[i], sep = "\t", header = TRUE)
    setnames(temp_data, c("gene", "exon", paste0(names(temp_data)[3:ncol(temp_data)], "_", i)))
    counts_dt <- merge(counts_dt, temp_data, by = c("gene", "exon"), all = TRUE) }
  results <- lapply(seq_len(nrow(counts_dt)), function(idx) { row <- counts_dt[idx,]
    processed_info <- process.event(row, group)
    list(gene = row$gene, exon = row$exon, processed = processed_info) })
  # Convert results to a data.table or data.frame for better readability
  processed_results <- rbindlist(lapply(results, as.data.frame), fill = TRUE)
  setnames(processed_results, sub("^processed\\.", "", names(processed_results)))
  # now p.adjust the p values with bH and call the column padj
  processed_results$padj <- p.adjust(processed_results$pvalue, method = "BH")
  return(processed_results) }


# psi.files <- paste0("data/alternativesplicing/psi/sample", 1:10, ".psi")
# load("data/alternativesplicing/diff_psi_test.RData")
# test<-diff.splicing(psi.files, group)
# fwrite(test, "data/alternativesplicing/psi_data.csv")
# test[pvalue < 0.05, .N]
# test[padj < 0.05, .N]
