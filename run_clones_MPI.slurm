#!/bin/bash
#SBATCH -J run_clones
#SBATCH --nodes=100
#SBATCH --ntasks=2400
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=118000
#SBATCH --time=10:00:00
#SBATCH --constraint=HSW24
#SBATCH --output run_clones.output

module purge
module load intel/17.2
module load fftw3/3.3.5
module load openmpi/intel/2.0.1
module load fftw/3.3.6-pl2

srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./run_clones_MPI
