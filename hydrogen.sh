#!/bin/bash
#BSUB -n 4              #each  job run on 1 core
#BSUB -W 5:00          #job run 10 hour
#BSUB -J jobArray[1-14483] #job array list goes 1,2,3...1000
#BSUB -o /home/yw174/job/out/%J.%I.out        #lsf output file
#BSUB -e /home/yw174/job/err/%J.%I.err       #lsf error file
#BSUB -q mcore         #submit to "short" queue
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/pdb_data/addHdata
babel -h ${LSB_JOBINDEX}.pdb ${LSB_JOBINDEX}.pdb
babel -d ${LSB_JOBINDEX}.pdb ${LSB_JOBINDEX}.pdb
obminimize -cg -ff MMFF94 -h -n 500 ${LSB_JOBINDEX}.pdb > ${LSB_JOBINDEX}_hydro.pdb