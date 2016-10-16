from mapping import copy_pdbfile
from Config import pdb_PREFIX,PDB_tar
import os

if __name__ == '__main__':
    for pdb in PDB_tar:
        filepos = os.path.join(pdb_PREFIX,pdb+'.pdb.gz')
        copy_pdbfile(filepos,os.path.join(pdb_PREFIX,pdb+'.pdb'))
        print 'Untar ' + pdb+ ' successfully!'

    os.system('rm /home/yw174/pdb_data/pdb_raw/*.gz')
    print 'All unneeded file is deleted'