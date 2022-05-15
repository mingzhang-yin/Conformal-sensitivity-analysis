#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUT_SUFFIX="gmm.out"
RUN_SCRIPT="bash/exp2/run_fct.sh"
Dt="homo hete"

for DTYPE in ${Dt}; do
  export DTYPE=${DTYPE}
  for G in $(seq 1. 0.5 4.0); do
    export SA=$G
    export NAME=bash
    export OUTNAME=${OUT_SUFFIX}
  sbatch --job-name=${NAME} \
  --output=output/${OUTNAME}.${SA}.${DTYPE} \
  --exclude=gonzo \
  ${RUN_SCRIPT}

  done
done

