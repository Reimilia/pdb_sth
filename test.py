'''
This file shows examples on how to use the scripts in project
'''

from fileparser import do_one_pdb
from Config import *
import os,sys,io
import gzip
from vector_gen import pdb_container,fake_pdb_container
from prody import *
from mapping import *


N=11549

prefix = '/home/yw174/pdb_data/addHdata'

os.remove('error.txt')

Succ = 0
Fail = 0

for i in range(N):
    try:
        file = os.path.join(prefix, '{}_hydro.pdb'.format(i+1))
        parsePDB(file)
        print str(i+1) + ' is OK'
        Succ+=1
    except:
        print str(i+1) + ' fails'
        with open('error.txt','a') as w:
            w.write(str(i+1)+'\n')
        Fail+=1

print 'Succ: {}, Fail: {}'.format(Succ,Fail)

'''

A= pdb_container('1avd',filepos='media/wy/data/pdb_raw/1avd.pdb.gz')

filename= 'media/wy/data/fast/1avd/1avd_248_ligand.pdb'
A.add_ligands(filename,'fast')
'''