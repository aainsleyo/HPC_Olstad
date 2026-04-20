#!/bin/bash
#SBATCH --job-name=rmt_project3
#SBATCH --account=def-veen       #  supervisors account
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # 1 manager + 3 workers
#SBATCH --mem=8G
#SBATCH --time=1:00:00            # adjust for larger sizes
#SBATCH --output=slurm-%j.out

module load intel
module load intelmpi
module load mkl

EXE=./mainsim

NS=(100 500)
NDATS=(500 1000)

for N in "${NS[@]}"; do
  for NDAT in "${NDATS[@]}"; do

    OUTDIR="results/n${N}_ndat${NDAT}"
    mkdir -p "$OUTDIR"

    echo "Running n=$N ndat=$NDAT ..."

    echo "$N"    > input.txt
    echo "$NDAT" >> input.txt

    srun $EXE

    mv eigs  "$OUTDIR/eigs"
    mv stats "$OUTDIR/stats"

  done
done

echo "All runs complete."
