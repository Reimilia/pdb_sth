#!/bin/bash
export OMP_NUM_THREADS=1
export PATH=$PATH:/home/yw174/usr/babel/bin/
source /home/yw174/python_env/wy/bin/activate
cd /home/yw174/pdb_data/addHdata
obminimize -cg -ff MMFF94 -h -n 500 1.pdb > 1_hydro.pdb