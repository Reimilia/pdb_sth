from mapping import copy_pdbfile
from Config import pdb_PREFIX,PDB_tar
from source import PDB_protein_tar
import os

if __name__ == '__main__':
    index = 1
    for pdb in PDB_protein_tar:
        filepos = os.path.join(pdb_PREFIX,pdb+'.pdb.gz')
        copy_pdbfile(filepos,os.path.join('/home/yw174/pdb_data/addHdata_whole',str(index)+'.pdb'))
        index+=1
        print 'Untar ' + pdb+ ' successfully!'
