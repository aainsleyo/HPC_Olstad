Read me:

The main program is contained in the main.f90, manager.f90, and worker.f90
The shell script, run_local.sh, runs smaller matrix sizes on the personal computer and saved results to the "results" folder. Be careful about overwrite.
The shell script, run_sharcnet.sh, runs a test to make sure results come out correctly from sharcnet.
The shell script, run_sharcnet2.sh, runs the larger matrix sizes, and can be edited, just like run_sharcnet.sh, to have different matrix sizes, and different iteration lengths, but it has an automation loop that sets different parameters.

Results are saved to the "Results folder", Connor's computations are from a shared google drive link. 
In the results folder are the slurm.out files and folders of each of the ran matrix sizes with iteration count. Inside of these folders are Stats and Eigs. 
