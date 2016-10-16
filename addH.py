'''
Add H in PDB protein files
'''
import os
import commands
from prody import *
from Config import PDB_tar,pdb_PREFIX,temp_pdb_PREFIX

from_dir = pdb_PREFIX
to_dir = temp_pdb_PREFIX


def add_hydrogens(filedir,pdbfilename,pdb):
    real_dir = os.path.join(filedir,pdbfilename)

    #cmd = 'babel -h {0} {0}'.format(real_dir)
    #os.system(cmd)
    #cmd = 'babel -d {0} {0}'.format(real_dir)
    #os.system(cmd)

    cmd = 'obminimize -cg -ff MMFF94 -h -n 500 {0}.pdb > {0}_hydro.pdb'.format(pdbfilename.split('.')[:-1])
    stat, out = commands.getstatusoutput(cmd)
    # If anything goes wrong , return False
    if stat == 256:
        print out
        return False

    return True

def split_receptors(pdbname,src,tardir):
    print src
    try:
        parse= parsePDB(src)
    except:
        return 0
    print 'here'

    if parse.select('nucleic') is not None:
        return 1
    protein = parse.select('protein')
    if protein is None:
        return 2


    return 3

def write_it(pdbname,src,tardir,index):
    parse = parsePDB(src)
    protein = parse.select('protein')
    if not os.path.exists(tardir):
        os.makedirs(tardir)
    writePDB(os.path.join(tardir,str(index)+'.pdb'),protein)
    return True

if __name__=='__main__':

    A =[]
    N =[]
    P =[]

    for i in range(len(PDB_tar)):
        pdb=PDB_tar[i]
        flag= split_receptors(pdb,os.path.join(from_dir,pdb+'.pdb.gz'),to_dir)
        if flag==3:
            print pdb
            A.append(pdb)
            write_it(pdb,os.path.join(from_dir,pdb+'.pdb.gz'),to_dir,len(A))
        if flag==2:
            P.append(pdb)
        if flag==1:
            N.append(pdb)

    with open('source.py','w') as f:
        f.write('PDB_protein_tar='+str(A)+'\n')
        f.write('Nucleic_tar='+str(N)+'\n')
        f.write('Unknown='+str(P)+'\n')

