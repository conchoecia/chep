#!/bin/bash
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
SNAKEFILE="${SCRIPTPATH}/../scripts/snakemake_scripts/Snakefile_chep_genome_analysis"
ARGSTRING=""

CONFIGFILE="${PWD}/config.yaml"

if [ "$#" -eq 0 ]
then
   snakemake --snakefile "${SNAKEFILE}" --config help=help
   exit 1
fi

INCONFIG=0
while (( "$#" )); do
  echo "$1"
  if [ "$1" == "--help" ]
  then
      snakemake --snakefile "${SNAKEFILE}" --config help=help
      exit 1
  fi
  if [ "$1" == "-h" ]
  then
      snakemake --snakefile "${SNAKEFILE}" --config help=help
      exit 1
  fi
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
