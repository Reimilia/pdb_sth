#!/bin/bash
#BSUB -n 1
#BSUB -W 168:00
#BSUB -J job_monitor
#BSUB -o /home/yw174/job/monitor.out
#BSUB -e /home/yw174/job/monitor.err
#BSUB -q priority
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/program/pdb_sth/bash
python job_monitor.py