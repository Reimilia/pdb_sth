#!/bin/bash
#BSUB -n 2
#BSUB -W 10:00
#BSUB -J clean
#BSUB -o /home/yw174/job/clean.out
#BSUB -e /home/yw174/job/clean.err
#BSUB -q short
export OMP_NUM_THREADS=1
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/program/pdb_sth/bash
python clean.py