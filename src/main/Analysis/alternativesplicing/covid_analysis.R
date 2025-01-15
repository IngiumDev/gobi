library("data.table")
library("dplyr")
library("ggplot2")
library("DESeq2")
library("pheatmap")

# https://web.stanford.edu/class/bios221/book/08-chap.html
info_dt <- fread("data/alternativesplicing/Corona/sample.list")
# split condition column into two by _
info_dt <- info_dt[, c("condition", "hours") := tstrsplit(condition, "_", fixed = TRUE)]
info_dt$hours <- as.factor(sub("h", "", info_dt$hours))
info_dt$condition <- as.factor(info_dt$condition)

covid.analysis <- function(counts_file, dataset_name, output_path) {
  star_counts <- fread(counts_file)
  # remove .bam from colnames so we can merge
  colnames(star_counts) <- gsub(".bam", "", colnames(star_counts))
  # def a function with a coutns as input

  star_deseq_m <- DESeqDataSetFromMatrix(
    countData = as.matrix(star_counts, rownames = "Geneid"),
    colData = info_dt,
    design = ~condition)
  star_deseq <- DESeq(star_deseq_m, parallel = TRUE)
  res_star <- results(star_deseq)
  # p-value histogram
  ggplot(as(res_star, "data.frame"), aes(x = pvalue)) +
    geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0) +
    ggtitle(paste0("p-value distribution for ", dataset_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(output_path, "pvalue_histogram_", dataset_name, ".pdf"), limitsize = FALSE)
  # MA plot: fold change versus mean of size-factor normalized counts. Logarithmic scaling is used for both axes. By default, points are colored red if the adjusted p-value is less than 0.1. Points which fall out of the
  # -axis range are plotted as triangles.
  # Set scipen to a high value to avoid scientific notation
  options(scipen = 999)

  # Create the MA plot
  png(filename = paste0(output_path, "MA_plot_", dataset_name, ".png"), width = 8, height = 6, units = "in", res = 300)
  plotMA(star_deseq, ylim = c(-10, 10), main = paste0("MA Plot for ", dataset_name, " DESeq Analysis"))
  dev.off()

  # Reset scipen to default if needed
  options(scipen = 0)
  # Figure 8.6: PCA plot. The  samples are shown in the 2D plane spanned by their first two principal components.
  star_rlog <- rlogTransformation(star_deseq)
  plotPCA(star_rlog, intgroup = c("condition", "hours")) + coord_fixed()

  pca_plot <- plotPCA(star_rlog, intgroup = c("condition", "hours")) +
    coord_fixed() +
    aes(color = hours, shape = condition) +
    ggtitle(paste0("PCA Plot for ", dataset_name, " Analysis"))
  # Save the PCA plot as a PDF
  ggsave(filename = paste0(output_path, "PCA_plot_", dataset_name, ".pdf"),
         plot = pca_plot,
         device = "pdf",
         width = 8, height = 6)

  star_lrd <- rlog(star_deseq, blind = FALSE)
  select <- head(order(rowVars(assay(star_lrd), useNames = FALSE), decreasing = TRUE), 30)
  mat <- assay(star_lrd)[select,]

  mat <- mat - rowMeans(mat)

  annotation_col <- as.data.frame(colData(star_lrd)[, c("condition", "hours")])

  heatmap <- pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    main = paste0("Top 30 Variable Genes (rlog transformed) for ", dataset_name)
  )
  pdf(file = paste0(output_path, "Heatmap_", dataset_name, ".pdf"), width = 8, height = 6)
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  dev.off()

  twoFactor <- star_deseq_m
  design(twoFactor) <- formula(~hours + condition)
  twoFactor <- DESeq(twoFactor, parallel = TRUE)
  res2 = results(twoFactor)
  trsf = function(x) ifelse(is.na(x), 0, (-log10(x))^(1 / 6))
  comparison_plot <- ggplot(tibble(pOne = res_star$pvalue,
                                   pTwo = res2$pvalue),
                            aes(x = trsf(pOne), y = trsf(pTwo))) +
    geom_hex(bins = 75) +
    coord_fixed() +
    xlab("Single factor analysis (condition)") +
    ylab("Two factor analysis (type + condition)") +
    geom_abline(col = "orange") +
    ggtitle("Comparison of Single vs. Two Factor Analysis")

  # Save the plot as a PNG file
  ggsave(filename = paste0(output_path, "Comparison_Plot_", dataset_name, ".png"),
         plot = comparison_plot,
         dpi = 300,
         width = 8, height = 6)
  compareRes = table(
    `simple analysis` = res_star$padj < 0.1,
    `two factor` = res2$padj < 0.1)
  addmargins(compareRes)
  vsp = varianceStabilizingTransformation(star_deseq)
  j = 1
  # Graph of variance-stabilizing transformation for the data of one of the samples, and for comparison also of the
  #  transformation. The variance-stabilizing transformation has finite values and finite slope even for counts close to zero, whereas the slope of
  #  becomes very steep for small counts and is undefined for counts of zero. For large counts, the two transformation are essentially the same.
  ggplot(
    tibble(
      counts = rep(assay(star_deseq)[, j], 2),
      transformed = c(
        assay(vsp)[, j],
        log2(assay(star_deseq)[, j])
      ),
      transformation = rep(c("VST", "log2"), each = nrow(star_deseq))
    ),
    aes(x = counts, y = transformed, col = transformation)) +
    geom_line() +
    xlim(c(0, 600)) +
    ylim(c(0, 9))
  # Create the graph with a title
  vst_plot <- ggplot(
    tibble(
      counts = rep(assay(star_deseq)[, j], 2),
      transformed = c(
        assay(vsp)[, j],
        log2(assay(star_deseq)[, j])
      ),
      transformation = rep(c("VST", "log2"), each = nrow(star_deseq))
    ),
    aes(x = counts, y = transformed, col = transformation)
  ) +
    geom_line() +
    xlim(c(0, 600)) +
    ylim(c(0, 9)) +
    ggtitle("Variance-Stabilizing Transformation vs. Log2 Transformation")

  # Save the plot as a PDF file
  ggsave(filename = paste0(output_path, "VST_Comparison_", dataset_name, ".pdf"),
         plot = vst_plot,
         device = "pdf",
         width = 8, height = 6)
  par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))

  # Define the myMA function with a filename parameter
  myMA <- function(h, v, theta = 0.5) {
    plotMA(star_deseq, lfcThreshold = theta, altHypothesis = h, ylim = c(-2.5, 2.5))
    abline(h = v * theta, col = "dodgerblue", lwd = 2)
  }

  # Open a PNG device to save all plots in one image
  png(filename = paste0(output_path, "All_MA_plots_", dataset_name, ".png"),
      width = 12, height = 10, units = "in", res = 300)

  # Set the layout to 2 rows and 2 columns
  par(mfrow = c(2, 2))

  # Plot each variant of the MA plot
  myMA("greaterAbs", c(-1, 1))
  myMA("lessAbs", c(-1, 1))
  myMA("greater", 1)
  myMA("less", -1)

  # Close the device (save the plot)
  dev.off()
}

covid.analysis("data/alternativesplicing/Corona/gene.counts.star", "star", "src/main/Analysis/alternativesplicing/corona/")
covid.analysis("data/alternativesplicing/Corona/gene.counts.hisat", "hisat", "src/main/Analysis/alternativesplicing/corona/")