# Load required packages
library(ontologyIndex)  # for reading and traversing the GO OBO file
library(data.table)     # for fast file I/O and manipulation
library(fgsea)          # for performing GSEA

# --- Step 1. Parse the GO ontology (OBO file) ---
# Read the GO ontology and automatically propagate "is_a" relationships.
# (Note: The 'get_ontology' function will include ancestors in the returned list.)
go_obo_path <- "data/GOEnrichment/go.obo"  # set your OBO file path
go <- get_ontology(go_obo_path, propagate_relationships = "is_a")
# (This object 'go' now contains elements such as 'ancestors' for each term.)
# For example, go$ancestors[["GO:0006955"]] returns all ancestors of that term.
# [ [oai_citation_attribution:0‡rdrr.io](https://rdrr.io/github/DanWiebe/fsgor/man/prepare_annotaion.html)]

# --- Step 2. Read and filter the GAF mapping file ---
# Here we use fread with a system pipe to read the gzipped file directly.
# The GAF file has comment lines starting with '!', so we skip those.
gaf_path <- "data/GOEnrichment/goa_human.gaf"  # Adjust the file path as needed
gaf <- fread(gaf_path, sep = "\t", header = FALSE)
# For this assignment, we need columns 3 (gene id), 4 (association qualifier), and 5 (GO term)# For this assignment, we need columns 3 (gene id), 4 (association qualifier), and 5 (GO term)
gaf_subset <- gaf[, .(gene = V3, qualifier = V4, go_term = V5)]
# Keep only rows without any qualifier modifier (i.e. where qualifier is empty)
gaf_subset <- gaf_subset[qualifier == ""]

# --- Step 3. Propagate gene annotations using the GO DAG ---
# Build a gene-to-GO mapping from the GAF file.
gene2go <- split(gaf_subset$go_term, gaf_subset$gene)

# For each gene, propagate its GO terms up to all ancestors.
gene2go_propagated <- lapply(gene2go, function(terms) {
  # For each annotated term, combine it with its ancestors (if any)
  all_terms <- unique(unlist(lapply(terms, function(term) {
    # Some terms may not be in the ontology (e.g. if obsolete) so check first.
    if (!is.null(go$ancestors[[term]])) {
      c(term, go$ancestors[[term]])
    } else {
      term
    }
  })))
  return(all_terms)
})

# Now, invert the mapping to create gene sets (GO term to genes)
go2genes <- list()
for (gene in names(gene2go_propagated)) {
  for (term in gene2go_propagated[[gene]]) {
    go2genes[[term]] <- c(go2genes[[term]], gene)
  }
}
# Remove duplicate gene entries from each GO term set
go2genes <- lapply(go2genes, unique)

# Optionally, you might wish to restrict to a specific GO namespace (e.g., biological_process).
# For that, use the go$namespace vector, e.g.,
# go2genes <- go2genes[names(go2genes) %in% names(go$namespace)[go$namespace == "biological_process"]]

# --- Step 4. Load the differential expression data ---
# The diffexp file is assumed to be a tab-delimited file with columns:
# id (gene id), fc (log2 fold change), and signif (boolean flag)
diffexp_path <- "data/GOEnrichment/simul_exp_go_bp_ensembl.tsv"  # adjust path accordingly
# Read the file without the comment.char argument
# Read the file into a character vector
diffexp_lines <- readLines(diffexp_path)

# Filter out lines starting with '#'
diffexp_lines <- diffexp_lines[!grepl("^#", diffexp_lines)]

# Write the filtered lines to a temporary file
temp_file <- tempfile()
writeLines(diffexp_lines, temp_file)

# Read the filtered file into a data table
diffexp <- fread(temp_file, sep = "\t", header = TRUE)

# Remove the temporary file
unlink(temp_file)

# Create a named numeric vector of fold changes (ranking metric)
ranks <- diffexp$fc
names(ranks) <- diffexp$id
# Sort in decreasing order (largest fold changes first)
ranks <- sort(ranks, decreasing = TRUE)
# Create a named numeric vector of fold changes (ranking metric)
ranks <- diffexp$fc
names(ranks) <- diffexp$id
# Sort in decreasing order (largest fold changes first)
ranks <- sort(ranks, decreasing = TRUE)

# --- Step 5. Run GSEA with fgsea ---
# Here, we filter gene sets by size. In this example, we only consider sets with between 50 and 500 genes.
fgseaRes <- fgsea(pathways = go2genes, stats = ranks, minSize = 50, maxSize = 500, nperm = 10000)
# Order results by adjusted p-value
fgseaRes <- fgseaRes[order(fgseaRes$padj),]

# --- Step 6. Compare to standard-of-truth ---
# Assume that the diffexp file (or a separate file) contains lines starting with '#' that indicate the true enriched GO IDs.
# Here we use an example placeholder vector; replace with your actual standard-of-truth GO ids.
true_terms <- c("GO:0006955", "GO:0002376")  # example GO IDs marked as truly enriched
# Add a column to flag if a GO term is in the true list
fgseaRes[, is_true := pathway %in% true_terms]

# --- Step 7. Write out results ---
output_path <- "go_enrichment_results.tsv"  # adjust as needed
fwrite(fgseaRes, file = output_path, sep = "\t")

# --- Optional: Print usage info if no command-line arguments are provided ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript my_gsea_analysis.R -obo <path/to/go.obo> -gaf <path/to/goa.gaf.gz> -enrich <path/to/diffexp.tsv> -o <output.tsv>\n")
}