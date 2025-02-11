#!/bin/bash

# Change directory: go up four levels so that the base directory contains "target" and "data"
cd ../../../../

# Define the main output folder
#out_dir="data/GOEnrichment/Runs"
out_dir="data/GOEnrichment/timingtest"
# Define the mapping types
mapping_types=("ensembl" "go")

# Define the namespaces: abbreviations and corresponding root terms
namespaces_abbr=("bp" "cc" "mf")
namespaces_root=("biological_process" "cellular_component" "molecular_function")

# Define the jar files to test
jars=("jars/GOEnrichmentRunner-big.jar" "jars/GOEnrichmentRunner-bignonparallel.jar")

# Number of iterations for averaging runtime
iterations=5

# Calculate the total number of iterations
total_iterations=$(( ${#mapping_types[@]} * ${#namespaces_abbr[@]} * ${#jars[@]} * iterations ))
current_iteration=0

# Loop over mapping types
for mapping in "${mapping_types[@]}"; do

  # Set the mapping file based on the mapping type
  if [ "$mapping" = "ensembl" ]; then
    mapping_file="data/GOEnrichment/goa_human_ensembl.tsv"
  elif [ "$mapping" = "go" ]; then
    mapping_file="data/GOEnrichment/goa_human.gaf.gz"
  fi

  # Loop over namespaces
  for i in "${!namespaces_abbr[@]}"; do
    ns_abbr="${namespaces_abbr[$i]}"
    ns_root="${namespaces_root[$i]}"

    # Loop over each jar file
    for jar in "${jars[@]}"; do
      # Extract a simple jar identifier (e.g., GOEnrichmentRunner-big or GOEnrichmentRunner-bignonparallel)
      jar_base=$(basename "$jar" .jar)

      # Run the configuration multiple times
      for run in $(seq 1 $iterations); do
        # Increment the overall iteration counter and print progress
        current_iteration=$(( current_iteration + 1 ))
        echo "Iteration: ${current_iteration}/${total_iterations}"

        # Define output files unique for each run using the main output folder variable
        out_file="${out_dir}/out_${mapping}_${ns_abbr}_${jar_base}_run${run}.tsv"
        overlap_file="${out_dir}/overlap_${mapping}_${ns_abbr}_${jar_base}_run${run}.tsv"
        log_file="${out_dir}/log_${mapping}_${ns_abbr}_${jar_base}_run${run}.txt"

        # Run the jar with the constructed parameters.
        # Note: The -analysis argument has been removed.
        java -Dlogback.configurationFile=logback.xml -jar "$jar" \
          -obo data/GOEnrichment/go.obo \
          -root "$ns_root" \
          -mappingtype "$mapping" \
          -mapping "$mapping_file" \
          -minsize 50 \
          -maxsize 500 \
          -enrich data/GOEnrichment/simul_exp_go_bp_ensembl.tsv \
          -o "$out_file" \
          -overlapout "$overlap_file" \
          > "$log_file" 2>&1

      done
    done
  done
done