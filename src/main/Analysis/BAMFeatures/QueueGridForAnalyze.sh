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
analysisscriptpath=$4
outputpath=$5
lengthbasepath=$6
annot_path=$7


# Read the reference table and process each line
tail -n +2 "$ref_table" | while IFS=$'\t' read -r bam gtf strandness; do
  # Extract the base names
  bamname=$(basename "$bam" .bam)
  gtfname=$(basename "$gtf" .gtf)

  # Construct the java command
  python_command="python3.12 $analysisscriptpath --input $annot_path/$bamname.annot --annotation $bamname --lengths $lengthbasepath/$gtfname.lengths --output $outputpath"


  # Submit the grid job
#  echo "python_command: $python_command"
  qsub -b y -wd "$working_directory" -M "$email" -m bea -N "$bamname" -p -200 -l vf=13000M "$python_command"
done