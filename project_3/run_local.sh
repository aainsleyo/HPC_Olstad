#!/bin/bash

NP=4                  # number of MPI processes (1 manager + 3 workers)
EXE=./mainsim

NS=(100 200 500 750 1000 1500)      # matrix sizes to test
NDATS=(500 1000 1500)      # number of realizations to test

for N in "${NS[@]}"; do
  for NDAT in "${NDATS[@]}"; do

    OUTDIR="results/n${N}_ndat${NDAT}"
    mkdir -p "$OUTDIR"

    echo "Running n=$N ndat=$NDAT ..."

    # Write n and ndat into a temporary input file read by init()
    echo "$N"    > input.txt
    echo "$NDAT" >> input.txt

    mpirun -np $NP $EXE

    # Move output files into the results directory
    mv eigs  "$OUTDIR/eigs"
    mv stats "$OUTDIR/stats"

    echo "  -> saved to $OUTDIR"
  done
done

echo "All runs complete."
