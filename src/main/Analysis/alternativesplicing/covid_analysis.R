library("data.table")
library("dplyr")
library("ggplot2")
library("DESeq2")
library("pheatmap")


info_dt <- fread("data/alternativesplicing/Corona/sample.list")
# split condition column into two by _
info_dt <- info_dt[, c("condition", "hours") := tstrsplit(condition, "_", fixed = TRUE)]
info_dt$hours <- as.factor(sub("h", "", info_dt$hours))
info_dt$condition <- as.factor(info_dt$condition)

star_counts <- fread("data/alternativesplicing/Corona/gene.counts.star")
hisat2_counts <- fread("data/alternativesplicing/Corona/gene.counts.hisat")
# remove .bam from colnames so we can merge
colnames(star_counts) <- gsub(".bam", "", colnames(star_counts))
colnames(hisat2_counts) <- gsub(".bam", "", colnames(hisat2_counts))
# def a function with a coutns as input
star_deseq_m <- DESeqDataSetFromMatrix(
  countData = as.matrix(star_counts, rownames = "Geneid"),
  colData = info_dt,
  design = ~condition)
star_deseq <- DESeq(star_deseq_m, parallel = TRUE)
res_star <- results(star_deseq)
#plotting
ggplot(as(res_star, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
plotMA(star_deseq, ylim = c(-2, 2))
star_rlog <- rlogTransformation(star_deseq)
plotPCA(star_rlog, intgroup = c("condition", "hours")) + coord_fixed()
select = order(rowMeans(assay(star_rlog)), decreasing = TRUE)[1:30]
pheatmap(assay(star_rlog)[select,],
         scale = "row",
         annotation_col = as.data.frame(
           colData(star_rlog)[, c("condition", "hours")]))
twoFactor <- star_deseq_m
design(twoFactor) <- formula(~hours + condition)
twoFactor <- DESeq(twoFactor, parallel = TRUE)
res2 = results(twoFactor)
trsf = function(x) ifelse(is.na(x), 0, (-log10(x))^(1 / 6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
       aes(x = trsf(pOne), y = trsf(pTwo))) +
  geom_hex(bins = 75) +
  coord_fixed() +
  xlab("Single factor analysis (condition)") +
  ylab("Two factor analysis (type + condition)") +
  geom_abline(col = "orange")
compareRes = table(
  `simple analysis` = res$padj < 0.1,
  `two factor` = res2$padj < 0.1)
addmargins(compareRes)
vsp = varianceStabilizingTransformation(star_deseq)
j = 1
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
par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))

myMA = function(h, v, theta = 0.5) {
  plotMA(star_deseq, lfcThreshold = theta, altHypothesis = h,
         ylim = c(-2.5, 2.5))
  abline(h = v * theta, col = "dodgerblue", lwd = 2)
}

myMA("greaterAbs", c(-1, 1))
myMA("lessAbs", c(-1, 1))
myMA("greater", 1)
myMA("less", -1)


# round 2
star_deseq_m <- DESeqDataSetFromMatrix(
  countData = as.matrix(hisat2_counts, rownames = "Geneid"),
  colData = info_dt,
  design = ~condition)
star_deseq <- DESeq(star_deseq_m, parallel = TRUE)
res_star <- results(star_deseq)
#plotting
ggplot(as(res_star, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
plotMA(star_deseq, ylim = c(-2, 2))
star_rlog <- rlogTransformation(star_deseq)
plotPCA(star_rlog, intgroup = c("condition", "hours")) + coord_fixed()
select = order(rowMeans(assay(star_rlog)), decreasing = TRUE)[1:30]
pheatmap(assay(star_rlog)[select,],
         scale = "row",
         annotation_col = as.data.frame(
           colData(star_rlog)[, c("condition", "hours")]))
twoFactor <- star_deseq_m
design(twoFactor) <- formula(~hours + condition)
twoFactor <- DESeq(twoFactor, parallel = TRUE)
res2 = results(twoFactor)
trsf = function(x) ifelse(is.na(x), 0, (-log10(x))^(1 / 6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
       aes(x = trsf(pOne), y = trsf(pTwo))) +
  geom_hex(bins = 75) +
  coord_fixed() +
  xlab("Single factor analysis (condition)") +
  ylab("Two factor analysis (type + condition)") +
  geom_abline(col = "orange")
compareRes = table(
  `simple analysis` = res$padj < 0.1,
  `two factor` = res2$padj < 0.1)
addmargins(compareRes)
vsp = varianceStabilizingTransformation(star_deseq)
j = 1
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
par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))

myMA = function(h, v, theta = 0.5) {
  plotMA(star_deseq, lfcThreshold = theta, altHypothesis = h,
         ylim = c(-2.5, 2.5))
  abline(h = v * theta, col = "dodgerblue", lwd = 2)
}

myMA("greaterAbs", c(-1, 1))
myMA("lessAbs", c(-1, 1))
myMA("greater", 1)
myMA("less", -1)