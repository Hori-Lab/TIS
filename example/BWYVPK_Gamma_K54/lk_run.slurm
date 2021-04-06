#!/bin/bash
#SBATCH -J Mg1.0
#SBATCH -o ./slurm.out
#SBATCH -e ./slurm.err
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 240:00:00
#SBATCH --mem=500m

### If your code uses a lot of stack memory (often happens by Fortran)
ulimit -s unlimited
export OMP_NUM_THREADS=4
export KMP_STACKSIZE=250m

### Move there
#cd /home/user/my_project5

### Run the code
../../md_HTN ./BWYVPK_Gamma.inp 1> log 2> err
