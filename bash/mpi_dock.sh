#!/bin/bash
#BSUB -a openmpi
#BSUB -n 100
#BSUB -W 100:00
#BSUB -J mpi_fast
#BSUB -o /home/mdk24/job/mpi_fast.out
#BSUB -e /home/mdk24/job/mpi_fast.err
#BSUB -q mpi
module load dev/openmpi-1.8.6
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/mdk24/wy_program/pdb_sth
mpirun -np 2 python dockingbroadcast.py