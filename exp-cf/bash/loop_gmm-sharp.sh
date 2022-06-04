#!/bin/bash
# call this at exp1 folder level

TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUT_SUFFIX="gmm.out"
RUN_SCRIPT="bash/run_gmm-sharp.sh"
Dt="homo hete"
if [ -d "$output" ]; then rm -Rf $output; fi
mkdir output
for DTYPE in ${Dt}; do
  export DTYPE=${DTYPE}
  for G in $(seq 1. 0.5 4.0); do
    for ((trial=1;trial<=20;trial+=1)); do
      export SA=$G
      export NAME="${DTYPE}_${G}_${trial}"
      export TRIAL=$trial
      export OUTNAME=${OUT_SUFFIX}
    sbatch --job-name=${NAME} \
    --output=output/${trial}.${SA}.${DTYPE} \
    --exclude=gonzo\
    ${RUN_SCRIPT}
    done
  done
done

