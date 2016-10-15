#!/bin/bash
#BSUB -n 2
#BSUB -W 10:00
#BSUB -J jobArray[2001-3000]
#BSUB -o /home/yw174/job/out/%J.%I.out
#BSUB -e /home/yw174/job/err/%J.%I.err
#BSUB -q short
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/pdb_data/addHdata
babel -h ${LSB_JOBINDEX}.pdb ${LSB_JOBINDEX}.pdb
babel -d ${LSB_JOBINDEX}.pdb ${LSB_JOBINDEX}.pdb
obminimize -cg -ff MMFF94 -h -n 500 ${LSB_JOBINDEX}.pdb > ${LSB_JOBINDEX}_hydro.pdb