import sys
import os
import mpi4py.MPI as MPI
import numpy as np
from main import bindingDB_pdb_tar_generator,docking_csv_generator
from fileparser import do_one_pdb,initiate_report,quick_split
from Config import PDB_tar,pdb_PREFIX
from job_dispatcher import *
'''
    This program use MPI to fulfill multiprocessing need with a dirty way.
    Just divide the PDB list into pieces and broadcast them to all processes
'''

#
#  Global variables for MPI
#

# instance for invoking MPI relatedfunctions
comm = MPI.COMM_WORLD
# the node rank in the whole community
comm_rank = comm.Get_rank()
# the size of the whole community, i.e.,the total number of working nodes in the MPI cluster
comm_size = comm.Get_size()



def get_file_list():
    filename='/home/yw174/pdb_data/input.txt'
    with open(filename,'rb') as f:
        list_dir =[]
        for each in f.readlines():
            list_dir.append(each)
    return list_dir


if __name__ == '__main__':
    #dirty and lazy way to multiprocessing
    #just split files to each processors.
    #and you will see multiprocessing
    if comm_rank == 0:
        # the No.0 one hand out issues
        file_list = get_file_list()
        file_num = len(file_list)
        sys.stderr.write("%d files\n" % len(file_list))
    else:
        file_num =0
    # broadcast filelist
    file_list = comm.bcast(file_list if comm_rank == 0 else None, root=0)
    file_num = len(file_list)
    local_files_offset = np.linspace(0, file_num, comm_size + 1).astype('int')

    # receive own part
    local_files = file_list[local_files_offset[comm_rank]:local_files_offset[comm_rank + 1]]
    sys.stderr.write("%d/%d processor gets %d/%d data \n" % (comm_rank, comm_size, len(local_files), file_num))
    #report_name = 'report.csv'

# Do seperately
for file_name in local_files:
    pdb = file_name.split('_')[0]
    resid = file_name.split('_')[1].rstrip('\n')
    print 'try to do %s_%s' % (pdb, resid)
    docking_csv_generator(pdb,resid)