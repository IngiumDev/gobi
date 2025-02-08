#!/bin/bash

# Change directory: go up four levels so that the base directory contains "target" and "data"
cd ../../../../

# Define the mapping types
mapping_types=("ensembl" "go")

# Define the namespaces: abbreviations and corresponding root terms
namespaces_abbr=("bp" "cc" "mf")
namespaces_root=("biological_process" "cellular_component" "molecular_function")

# Loop over mapping types and namespaces
for mapping in "${mapping_types[@]}"; do

  # Set the mapping file based on the mapping type
  if [ "$mapping" = "ensembl" ]; then
    mapping_file="data/GOEnrichment/goa_human_ensembl.tsv"
  elif [ "$mapping" = "go" ]; then
    mapping_file="data/GOEnrichment/goa_human.gaf.gz"
  fi

  # Iterate over the three namespaces
  for i in "${!namespaces_abbr[@]}"; do
    ns_abbr="${namespaces_abbr[$i]}"
    ns_root="${namespaces_root[$i]}"


    # Define output file names that indicate the mapping and namespace
    out_file="data/GOEnrichment/Runs/out_${mapping}_${ns_abbr}.tsv"
    overlap_file="data/GOEnrichment/Runs/overlap_${mapping}_${ns_abbr}.tsv"
    analysis_file="data/GOEnrichment/Runs/analysis-out-${mapping}-${ns_abbr}.tsv"

    # Run the jar with the constructed parameters
    java -jar target/GOEnrichmentRunner.jar \
      -obo data/GOEnrichment/go.obo \
      -root "$ns_root" \
      -mappingtype "$mapping" \
      -mapping "$mapping_file" \
      -minsize 50 \
      -maxsize 500 \
      -enrich data/GOEnrichment/simul_exp_go_bp_ensembl.tsv \
      -o "$out_file" \
      -overlapout "$overlap_file" \
      -analysis "$analysis_file"
  done
done