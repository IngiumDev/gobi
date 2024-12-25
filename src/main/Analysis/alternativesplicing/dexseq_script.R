#'
#' Simple script to load data in DEXSeq and perform
#' Differential Exon Usage analysis
#'

library(DEXSeq)

## list the htseq-count files
# TODO possibly adjust path
countFiles = list.files("data/alternativesplicing/dexseq", full = T)
names(countFiles) <- gsub("_dexseq.txt", "", basename(countFiles))

## prepare the sample annotation
group <- rep(1, length(countFiles))
# sample names ending in 6..10 are in group 2
group[grepl("[06-9]$", names(countFiles))] <- 2
sampleTable = data.frame(condition = factor(group))
rownames(sampleTable) = names(countFiles)

## load the data
# TODO possibly adjust path to annotation
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = sampleTable,
  design = ~sample + exon + condition:exon,
  flattenedfile = "data/alternativesplicing/annotation_dexseq_b37.gff")

## estimate the size factors
dxd = estimateSizeFactors(dxd)
## fit function to compute the variance in dependence on the mean
dxd = estimateDispersions(dxd)
## check the fit of the variance (dispersion) model
plotDispEsts(dxd)
## run the test for differential exon usage
dxd = testForDEU(dxd)
## obtain the results
dxr1 = DEXSeqResults(dxd)