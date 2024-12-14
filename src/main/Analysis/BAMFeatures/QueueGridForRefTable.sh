#!/bin/bash

# Check if the reference table is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <reference_table>"
  exit 1
fi
source /etc/sge61.sh
# Variables
ref_table=$1
working_directory=$2
email=$3
jarpath=$4
gtfbasepath="/mnt/biosoft/praktikum/genprakt-data/gtfs"
bamfeaturesbasepath="/mnt/biosoft/praktikum/genprakt-data/BamFeatures/complete_bams"
outputpath=$5

# Read the reference table and process each line
tail -n +2 "$ref_table" | while IFS=$'\t' read -r bam gtf strandness reference_solution; do
  # Extract the base names
  bamname=$(basename "$bam" .bam)
  gtfname=$(basename "$gtf")

  # Construct the java command
  java_cmd="java -Xmx10g -jar $jarpath -gtf $gtfbasepath/$gtfname -bam $bamfeaturesbasepath/$bam -o $outputpath/$bamname.annot"

  # Add the -frstrand parameter if strandness is not empty
  if [ -n "$strandness" ]; then
    java_cmd="$java_cmd -frstrand $strandness"
  fi

  # Submit the grid job
  qsub -b y -wd "$working_directory" -M "$email" -m bea -N "$bamname" -p -200 -l vf=13000M "$java_cmd"
done