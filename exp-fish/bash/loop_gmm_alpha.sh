#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUT_SUFFIX="gmm.out"
RUN_SCRIPT="run_fish.sh"

if [ -d "$output" ]; then rm -Rf $output; fi
mkdir output
for alpha in $(seq 0.05 0.05 0.5); do
  export ALPHA=${alpha}
  for G in $(seq 1. 0.2 4.0); do
    export SA=$G
    export NAME=bash
    export OUTNAME=${SA}.${ALPHA}
  sbatch --job-name=${NAME} \
  --output=output/${SA}.${ALPHA} \
  --exclude=gonzo \
  ${RUN_SCRIPT}

  done
done

