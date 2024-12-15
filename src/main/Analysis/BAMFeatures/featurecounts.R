if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)
library(data.table)

# Specify paths to input files
bam_file <- "data/bamfeatures/complete_bams/nookaew_cm.bam"  # Your BAM file
annotation_file <- "data/gtf/Saccharomyces_cerevisiae.R64-1-1.75.gtf"   # Your GTF file
output_file <- "data/bamfeatures/complete_bams/ebna_hisat.counts"           # Output file for counts

# Run featureCounts
counts <- featureCounts(
  files = bam_file,
  annot.ext = annotation_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,      # Equivalent to -p (paired-end reads)
  countMultiMappingReads = FALSE,  # Default behavior, ignore multimappers
  nthreads = 8,             # Use 4 threads
  allowMultiOverlap = TRUE,
  primaryOnly = TRUE,
  strandSpecific = 0,
  requireBothEndsMapped = TRUE,
  countReadPairs = TRUE,
  countChimericFragments = FALSE,
)
write.table(
  counts$counts,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
raw_counts <- counts$counts  # Raw counts matrix
gene_lengths <- counts$annotation$Length  # Lengths of genes (from featureCounts output)
library_size <- colSums(raw_counts)  # Total reads per library (column)
gene_names <- counts$annotation$GeneID  # Gene names

# Convert gene lengths from bp to kb
gene_lengths_kb <- gene_lengths / 1000

# Calculate RPKM
rpkm <- t(t(raw_counts) / (gene_lengths_kb * library_size / 1e6))
t<-as.data.table(rpkm)
rownames(t)<-gene_names
colnames(t)<-"RPKM"
# Sort the rows based on the 'rpkm' column in descending order
top_10_genes <- t[order(-rpkm), .(gene = .I)][1:10, gene]
top_10_genes <- rownames(t)[order(-t$RPKM)][1:10]

# Print the top 10 row names (gene identifiers)
print(cat(paste(top_10_genes, collapse = ", "), "\n"))


bam_file <- "data/bamfeatures/complete_bams/ebna_hisat.bam"  # Your BAM file
annotation_file <- "data/gtf/Homo_sapiens.GRCh37.75.gtf"   # Your GTF file
output_file <- "data/bamfeatures/complete_bams/ebna_hisat.counts"           # Output file for counts

# Run featureCounts
counts <- featureCounts(
  files = bam_file,
  annot.ext = annotation_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,      # Equivalent to -p (paired-end reads)
  countMultiMappingReads = FALSE,  # Default behavior, ignore multimappers
  nthreads = 8,             # Use 4 threads
  allowMultiOverlap = TRUE,
  primaryOnly = TRUE,
  strandSpecific = 1,
  requireBothEndsMapped = TRUE,
  countReadPairs = TRUE,
  countChimericFragments = FALSE,
)
write.table(
  counts$counts,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
raw_counts <- counts$counts  # Raw counts matrix
gene_lengths <- counts$annotation$Length  # Lengths of genes (from featureCounts output)
library_size <- colSums(raw_counts)  # Total reads per library (column)
gene_names <- counts$annotation$GeneID  # Gene names

# Convert gene lengths from bp to kb
gene_lengths_kb <- gene_lengths / 1000

# Calculate RPKM
rpkm <- t(t(raw_counts) / (gene_lengths_kb * library_size / 1e6))
t<-as.data.table(rpkm)
rownames(t)<-gene_names
colnames(t)<-"RPKM"
# Sort the rows based on the 'rpkm' column in descending order
top_10_genes <- t[order(-rpkm), .(gene = .I)][1:10, gene]
top_10_genes <- rownames(t)[order(-t$RPKM)][1:10]

# Print the top 10 row names (gene identifiers)
print(cat(paste(top_10_genes, collapse = ", "), "\n"))

bam_file <- "data/bamfeatures/complete_bams/hes_star.bam"  # Your BAM file
annotation_file <- "data/gtf/Homo_sapiens.GRCh37.75.gtf"   # Your GTF file
output_file <- "data/bamfeatures/complete_bams/ebna_hisat.counts"           # Output file for counts

# Run featureCounts
counts <- featureCounts(
  files = bam_file,
  annot.ext = annotation_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,      # Equivalent to -p (paired-end reads)
  countMultiMappingReads = FALSE,  # Default behavior, ignore multimappers
  nthreads = 8,             # Use 4 threads
  allowMultiOverlap = TRUE,
  primaryOnly = TRUE,
  strandSpecific = 2,
  requireBothEndsMapped = TRUE,
  countReadPairs = TRUE,
  countChimericFragments = FALSE,
)
write.table(
  counts$counts,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
raw_counts <- counts$counts  # Raw counts matrix
gene_lengths <- counts$annotation$Length  # Lengths of genes (from featureCounts output)
library_size <- colSums(raw_counts)  # Total reads per library (column)
gene_names <- counts$annotation$GeneID  # Gene names

# Convert gene lengths from bp to kb
gene_lengths_kb <- gene_lengths / 1000

# Calculate RPKM
rpkm <- t(t(raw_counts) / (gene_lengths_kb * library_size / 1e6))
t<-as.data.table(rpkm)
rownames(t)<-gene_names
colnames(t)<-"RPKM"
# Sort the rows based on the 'rpkm' column in descending order
top_10_genes <- t[order(-rpkm), .(gene = .I)][1:10, gene]
top_10_genes <- rownames(t)[order(-t$RPKM)][1:10]

# Print the top 10 row names (gene identifiers)
print(cat(paste(top_10_genes, collapse = ", "), "\n"))