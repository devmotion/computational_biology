#!/bin/bash

# Merge BAM files BB151124_I of single cell sequencing.
# Authors: Rebekka Mueller, David Widmann
# 2016-01-20
# v0.1

# Function that outputs correct usage of script
display_usage() {
  echo -e "Usage:"
  echo "$0 [<input> <output>]"
  echo "$0 -h | --help"
  echo "$0 --version"
  echo -e "\nOptions:"
  echo -e "-h --help  Show this screen."
  echo -e "--version  Show version."
  echo -e "\nArguments:"
  echo -e "<input>  Input folder with BAM files BB151124_I [default: .]."
  echo -e "<output>  Output folder for merged BAM file [default: .]."
}

# Display help message
if [[ "$@" == "-h" || "$@" == "--help" ]]; then
  echo -e "Merge BAM files BB151124_I of single cell sequencing.\n"
  display_usage
  exit 0
fi

# Display version
if [[ "$@" == "--version" ]]; then
  echo "v0.1 (2016-01-20)"
  exit 0
fi

# Set input and output directory
INPUT="${1:-.}"
OUTPUT="${2:-.}"

# Create output directory if necessary
mkdir -p "$OUTPUT"

# Merge single-cell BAM files without controls 49, 68, and 78 and bad samples 57, 64, 80, 85, 92, and 94
echo "Merging BAM files of best quality ..."
samtools merge "$OUTPUT/BB151124_I_merged_best.bam" "./BB151124_I_0"{{50..56},{58..63},{65..67},{69..77},79,{81..84},{86..91},93,95,96}".bam"

# Index merged BAM file (already sorted)
echo "Indexing merged BAM file ..."
samtools index "$OUTPUT/BB151124_I_merged_best.bam"

echo "Done. Merged BAM file written to \"$OUTPUT/BB151124_I_merged_best.bam\"."
