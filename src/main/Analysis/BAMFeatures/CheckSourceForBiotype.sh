#!/bin/bash
#count the lines where the source (Column 2) is not equal to the gene_biotype, and the third column is specifically gene
for gtf in *.gtf; do
    count=$(awk -F'\t' '$3 == "gene" && $2 != "" && $0 ~ /gene_biotype/ {
        split($9, attrs, ";");
        for (i in attrs) {
            if (attrs[i] ~ /gene_biotype/) {
                gsub(/gene_biotype |"| /, "", attrs[i]);
                if ($2 != attrs[i]) { print $0; }
            }
        }
    }' "$gtf" | wc -l);
    echo "$gtf: $count";
done