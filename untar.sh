#!/bin/bash
#BSUB -n 1
#BSUB -W 8:00
#BSUB -J untar
#BSUB -o /home/yw174/job/test.out
#BSUB -e /home/yw174/job/test.err
#BSUB -q priority
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/program/pdb_sth
python addH.py