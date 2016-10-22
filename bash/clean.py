from Config import temp_pdb_PREFIX
import os

if __name__=='__main__':
    filedir = os.listdir(temp_pdb_PREFIX)
    for files in filedir:
        if files != 'fake-ligand.pdb':
            os.system('rm -r '+ files)
            print files + ' has been wiped up'