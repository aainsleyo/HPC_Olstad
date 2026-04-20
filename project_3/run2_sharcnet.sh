#!/bin/bash
#SBATCH --job-name=rmt_project3
#SBATCH --account=def-veen
#SBATCH --nodes=1
#SBATCH --ntasks=17               # 1 manager + 16 workers
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=slurm-%j.out

module load intel
module load intelmpi
module load mkl

EXE=./mainsim
NS=(1000 2000 4000 6000)
NDATS=(5000 10000)

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
    # capture wall time from slurm output into a dedicated file
    grep "wtime" slurm-${SLURM_JOB_ID}.out | tail -1 >> "$OUTDIR/walltime"
  done
done

echo "All runs complete."
