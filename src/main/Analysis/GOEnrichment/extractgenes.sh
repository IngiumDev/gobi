#!/bin/bash
# extract_genes.sh
# This script extracts gene IDs for Enrichr analysis from a TSV file.
# The input file should have a header and three columns:
#   1. gene id
#   2. log2 fold change
#   3. signif (with values like "true" or "false")
#
# Usage:
#   ./extract_genes.sh input_file.tsv output_file.txt

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_file.tsv> <output_file.txt>"
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Skip the header and extract gene ids (first column) where the third column equals "true"
awk -F"\t" 'NR>1 && tolower($3)=="true" {print $1}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Extracted gene IDs have been saved to $OUTPUT_FILE."
