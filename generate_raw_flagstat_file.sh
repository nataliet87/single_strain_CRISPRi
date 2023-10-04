#!/bin/bash

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-I] [-O]" 2>&1
        echo '   -I   path to bam filelist'
        echo '   -O   path to output flagstat report file'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi


# Define list of arguments expected in the input
optstring="I:O:"

while getopts ${optstring} arg; do
  case "${arg}" in
#    I) bam_dir="${OPTARG}" ;;
    I) filelist="${OPTARG}" ;;
    O) report="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done


#for i in "$bam_dir"/*.sorted.bam

while read i
do
    #		echo "$i" >> "$bam_dir"/flagstat_eval.txt
		echo "$i" >> "$report"
		samtools flagstat "$i" >> "$report"

done < "$filelist"


