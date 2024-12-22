#!/bin/bash
source /etc/sge61.sh
# Variables
bampath="/mnt/biosoft/praktikum/genprakt/AlternativeSplicing/bam"
gtf="/mnt/biosoft/praktikum/genprakt/AlternativeSplicing/annotation_b37.gtf"
working_directory=$1
email=$2
jarpath=$3
outputpath=$4
for bam_file in "$bampath"/*.bam; do
  if [ -f "$bam_file" ]; then
    bamname=$(basename "$bam_file" .bam)

    java_cmd="java -Xmx10g -jar $jarpath -gtf $gtf -bam $bam_file -o $outputpath/$bamname.psi"
    qsub -b y -wd "$working_directory" -M "$email" -m bea -N "$bamname" -p -200 -l vf=13000M "$java_cmd"
  fi
done