# Project 3 — Random Matrix Eigenvalue Analysis

## Overview

Project 3 Workflow:
- Fortran simulation (matrix generation + eigenvalue extraction)
- Python analysis (plotting + statistical testing)
- Shell scripts for automation on local and HPC systems

## Code Structure

### Core Simulation (Fortran)
- `main.f90` - main driver program
- `manager.f90` - manages workload distribution
- `worker.f90` - performs matrix computations and eigenvalue extraction

## Shell Scripts

### Local Testing
- `run_local.sh`
  - Runs smaller matrix sizes on a personal computer
  - Saves results into the `results/` folder
  - e careful: outputs may be overwritten

### HPC (SHARCNET) Runs
- `run_sharcnet.sh`
  - Test run to verify correctness of HPC pipeline outputs

- `run_sharcnet2.sh`
  - Main production script for large matrix sizes
  - Can be modified to change:
    - matrix sizes
    - number of iterations
  - Includes an automation loop for running multiple parameter configurations

## Data & Results

### Output Location

Each run produces:
- `Eigs/` - eigenvalue data
- `Stats/` - computed summary statistics
- `slurm*.out` - HPC job logs

Each matrix size + iteration configuration is stored in its own folder.

### External Data
- Connor’s reference computations are stored via a shared Google Drive link (external dataset for comparison).

## Python Analysis

The file:
- `plotting.py`

Workflow:
- Histogram analysis of $\lambda_{\max}$
- Tracy–Widom normalization:
- Kernel density estimation (KDE)
- Gaussian comparison fits
- Statistical testing:
  - Shapiro–Wilk test
  - Kolmogorov–Smirnov tests (Normal vs Tracy–Widom)
  - Skewness and kurtosis

Outputs:
- `.pdf` and `.png` plots
- `output.log` with full statistics

## How to Run

Make executable:
```bash id="run1"
chmod +x run_analysis.sh
