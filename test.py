'''
This file shows examples on how to use the scripts in project
'''

from fileparser import do_one_pdb
from Config import *
import os,sys,io
import gzip
from vector_gen import pdb_container,fake_pdb_container
from mapping import *


#PDB_tar=['1uwj']

#do_one_pdb('100d')
#A= pdb_container('100d',filepos=os.path.join(pdb_PREFIX,'1j8q.pdb.gz'))
#print A.set_vina_benchmark('147')

#a= fake_pdb_container('aa2ar',filepos=fake_src_PREFIX+'aa2ar/'+fake_pdb_name)

#filenames =  os.listdir(fake_hetero_PREFIX)

#for filename in filenames:
#    if filename.split('.')[-1]=='pdb':
#        a.append_vectors(os.path.join(fake_hetero_PREFIX,filename))

#a.self_generating()

'''
P =0
N =0
C =0

ct=0

for pdb in PDB_tar:
    pdb= pdb.lower()
    A= pdb_container(pdb,filepos=os.path.join(pdb_PREFIX,pdb+'.pdb.gz'),OUT=False)
    pdbtype= A.get_pdb_type()
    if pdbtype=='Protein':
        P+=1
    if pdbtype=='Nucleic':
        N+=1
    if pdbtype=='Protein_Nucleic_Complex':
        C+=1
    ct+=1
    print str(ct)+'/14388'

print P
print N
print C
'''

#prepare_receptor('/home/wy/Documents/BCH_coding/pdb_data_extracter/1bib.pdb','1bib')
#I= pdb_container('1bib',filepos='/home/wy/Documents/BCH_coding/pdb_data_extracter/result/1bib/1bib.pdb',OUT=True)

#real_dir= os.path.join(pdb_PREFIX,'1bib.pdb.gz')
#repair_pdbfile(real_dir,'1bib',OVERWRITE=True)

#A= pdb_container('/media/wy/data/all/aa2ar/aa2ar_decoys_1_docked.pdb')

import pandas as pd

chunks = pd.read_csv('result/filter_1avd.csv',iterator = True)
chunk = chunks.get_chunk(5)
print chunk
