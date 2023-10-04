#!/bin/bash

# Get function arguments:
function usage {
        echo "Usage: $(basename $0) [-F] [-O]" 2>&1
        echo '   -F   file holding path to all fastq files to be trimmed'
        echo '   -O   output directory for trimmed fastq files'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="F:O:"

while getopts ${optstring} arg;
do
  case "${arg}" in
    F) FASTQ_LIST="${OPTARG}" ;;
    O) OUTPUT_DIR="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# Trim fastq
while read i
do
			/Users/nataliethornton/biotools/TrimGalore-0.6.5/trim_galore \
					--path_to_cutadapt /Users/nataliethornton/anaconda3/envs/wgs/bin/cutadapt \
					--output_dir "$OUTPUT_DIR" "$i"
					#--quality 20 --fastqc --output_dir "$OUTPUT_DIR" "$i"

done < "$FASTQ_LIST"


