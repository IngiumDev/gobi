if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")
# Load required packages
library(fgsea)

# Load required packages
library(ontologyIndex)
library(fgsea)

# -------------------------------
# Step 1: Load and Process the OBO File
# -------------------------------
obo_file <- "data/GOEnrichment/go.obo"  # Path to your OBO file
obo <- get_ontology(obo_file, propagate_relationships = c("is_a"), extract_tags = "everything")

# Filter GO terms to only those in the biological_process namespace.
bp_go_ids <- names(obo$namespace)[obo$namespace == "biological_process"]

# -------------------------------
# Step 2: Parse the GAF File
# -------------------------------
gaf_file <- "data/GOEnrichment/goa_human.gaf.gz"  # Path to your GAF file (gzipped)
# Read in the GAF file; lines starting with "!" are comments and are ignored.
gaf <- read.delim(gaf_file, header = FALSE, sep = "\t", comment.char = "!", stringsAsFactors = FALSE)

# The assignment specifies that only columns 3-5 are needed:
# Column 3: gene id, Column 4: association qualifier, Column 5: GO id.
gaf_sub <- gaf[gaf$V4 == "" & gaf$V5 %in% bp_go_ids, c(3, 5)]
colnames(gaf_sub) <- c("gene", "go")

# -------------------------------
# Step 3: Propagate Annotations and Build Gene Sets
# -------------------------------
# Initialize an empty list for gene sets.
geneSets <- list()

# For each annotation, propagate the gene association to all ancestors (including itself)
# that belong to the biological_process namespace.
# Initialize geneSets before the loop if not already done.
geneSets <- list()

for (i in seq_len(nrow(gaf_sub))) {
  gene <- gaf_sub$gene[i]
  go_term <- gaf_sub$go[i]

  # Get ancestors (including the term itself)
  all_terms <- unique(c(go_term, get_ancestors(obo, go_term)))

  # Filter ancestors to only those in the biological_process namespace.
  all_terms <- all_terms[sapply(all_terms, function(term) {
    !is.null(obo$namespace[[term]]) && obo$namespace[[term]] == "biological_process"
  })]

  # Add the gene to each relevant GO term's gene set.
  for (term in all_terms) {
    if (is.null(geneSets[[term]])) {
      geneSets[[term]] <- gene
    } else {
      geneSets[[term]] <- unique(c(geneSets[[term]], gene))
    }
  }
}  # <- This closing bracket ends the for loop.

# Add the gene to each of the GO term's gene set.
geneSets <- lapply(all_terms, function(term) {
  if (is.null(geneSets[[term]])) {
    geneSets[[term]] <- gene
  } else {
    geneSets[[term]] <- unique(c(geneSets[[term]], gene))
  }
  geneSets[[term]]
})

# Optionally, remove any gene sets that have no genes (if applicable)
geneSets <- geneSets[sapply(geneSets, length) > 0]

# --- Step 1: Read Differential Expression Data ---
# Replace "simul_exp_go_bp_ensembl.tsv" with your diffexp file path.
diffexp_file <- "data/GOEnrichment/simul_exp_go_bp_ensembl.tsv"
diffexp_data <- read.delim(diffexp_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")

# Inspect the data (assumed columns: id, fc, signif)
head(diffexp_data)

# --- Step 2: Create a Ranked Gene List ---
# Use all genes along with their log2 fold changes (fc) for ranking.
# This is preferred for ranking-based enrichment analysis.
geneList <- diffexp_data$fc
names(geneList) <- diffexp_data$id

# Sort the gene list in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# --- Step 3: Prepare Gene Sets and Run fgsea ---
# Assume that you have parsed your OBO file and mapping file to build a list of gene sets.
# For example, geneSets is a named list where each element is a vector of gene ids.
# Example:
# geneSets <- list("GO:0008150" = c("GeneA", "GeneB", "GeneC"),
#                  "GO:0009987" = c("GeneD", "GeneE", "GeneF"))
# Adjust minSize and maxSize as specified (e.g., 50 and 500).
fgseaRes <- fgsea(
  pathways = geneSets,
  stats = geneList,
  minSize = 50,
  maxSize = 500,
  nperm = 10000  # You can adjust the number of permutations as needed.
)

# fgsea typically returns an adjusted p-value field (padj).
# If you need to calculate additional BH corrections on other test p-values,
# you can do so as follows:

# For example, if your fgsea output had additional columns:
if ("hg_pval" %in% colnames(fgseaRes)) {
  fgseaRes$hg_fdr <- p.adjust(fgseaRes$hg_pval, method = "BH")
}
if ("fej_pval" %in% colnames(fgseaRes)) {
  fgseaRes$fej_fdr <- p.adjust(fgseaRes$fej_pval, method = "BH")
}
if ("ks_pval" %in% colnames(fgseaRes)) {
  fgseaRes$ks_fdr <- p.adjust(fgseaRes$ks_pval, method = "BH")
}

# View the top fgsea results
head(fgseaRes)

# --- Step 4: Read in Your Enrichment Results and Apply BH Correction ---
# Replace "your_results.tsv" with the path to your own enrichment result file.
results_file <- "your_results.tsv"
my_results <- read.delim(results_file, header = TRUE, stringsAsFactors = FALSE)

# Check the structure to confirm expected columns (e.g., term, name, size, is_true,
# noverlap, hg_pval, fej_pval, ks_stat, ks_pval, etc.)
head(my_results)

# Apply BH correction to the p-values in your results
my_results$hg_fdr <- p.adjust(my_results$hg_pval, method = "BH")
my_results$fej_fdr <- p.adjust(my_results$fej_pval, method = "BH")
my_results$ks_fdr <- p.adjust(my_results$ks_pval, method = "BH")

# Verify that the corrected values are added
head(my_results)