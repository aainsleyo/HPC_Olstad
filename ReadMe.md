# Matrix Multiplication (Fortran)

## Overview

This program performs matrix–matrix multiplication in Fortran using a naïve triple-loop algorithm. It is intended for performance testing and comparison across compiler optimization levels and matrix sizes in an HPC setting.

The program computes:
[
A = B \times C
]

and reports the elapsed CPU time.

---

## What the Program Does

* Dynamically allocates three `n × n` matrices
* Fills the input matrices with random values
* Performs matrix multiplication using a triple loop
* Times only the multiplication step
* Prints a checksum to prevent compiler optimization from removing the computation

---

## Compiling

Load the Intel compiler environment:

```bash
source /opt/intel/oneapi/setvars.sh
```

Compile without optimization:

```bash
ifx matrix.f90 -o matrix
```

Compile with optimization:

```bash
ifx -O2 matrix.f90 -o matrix_O2
ifx -O3 matrix.f90 -o matrix_O3
```

---

## Running

```bash
./matrix
```

or:

```bash
./matrix_O2
./matrix_O3
```

The program prints:

* Elapsed CPU time for the multiplication
* A checksum of the result matrix

---

## Changing Matrix Size

Edit the line in the source code:

```fortran
n = 4096
```

Common test sizes:

* 512
* 1024
* 2048
* 4096

Recompile after changing `n`.

---

## Notes

* The loop order is chosen to improve cache efficiency for Fortran’s column-major layout.
* Compiler optimization flags alone may not significantly change performance for this naïve algorithm.
* The checksum ensures the computation is not optimized away.

---
