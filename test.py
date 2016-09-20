from fileparser import do_one_pdb
from Config import pdb_PREFIX
import os,sys,io
import gzip

PDB_tar=['1hck']

#do_one_pdb('1hck')
do_one_pdb('1eob')

for pdb in PDB_tar:

    with gzip.open(os.path.join(pdb_PREFIX,'{}.pdb.gz'.format(pdb))) as f:
        for line in f:
            #print line
            pass