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

#def quick_split(pdb):
#    pdb = pdb.lower()
#    fake_pdb_container(pdb,filepos=os.path.join(pdb_PREFIX,pdb+'.pdb.gz'))

#quick_split('1avd')

A = pdb_container('1avd',filepos='/media/wy/data/pdb_raw/1avd.pdb.gz')
A.add_ligands('/media/wy/data/fast/1avd/1avd_248_ligand.pdb',suffix='fast')
