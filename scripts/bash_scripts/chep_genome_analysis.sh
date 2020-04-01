#!/bin/bash
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
SNAKEFILE="${SCRIPTPATH}/../scripts/snakemake_scripts/Snakefile_chep_genome_analysis"
ARGSTRING=""

#CONFIGFILE="${PWD}/config.yaml"
CONFIGFILE="test"

INCONFIG=0
while (( "$#" )); do
  if [ "${INCONFIG}" -eq 0 ]
  then
    if [ "$1" == "--configfile" ]
    then
        INCONFIG=1
    else
        ARGSTRING+=" $1"
    fi
  else
    CONFIGFILE="$1"
    INCONFIG=0
  fi
  shift
done

echo "Running the command: snakemake --configfile ${CONFIGFILE} --snakefile ${SNAKEFILE} ${ARGSTRING}"
snakemake --configfile "${CONFIGFILE}" --snakefile "${SNAKEFILE}" "${ARGSTRING}"
