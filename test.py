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


#PDB_tar=['1uwj']

#do_one_pdb('100d')
#A= pdb_container('100d',filepos=os.path.join(pdb_PREFIX,'1j8q.pdb.gz'))
#print A.set_vina_benchmark('147')

#a= fake_pdb_container('aa2ar',filepos=fake_src_PREFIX+'aa2ar/'+fake_pdb_name)

#filenames =  os.listdir(fake_hetero_PREFIX)

#for filename in filenames:
#    if filename.split('.')[-1]=='pdb':
#        a.append_vectors(os.path.join(fake_hetero_PREFIX,filename))
'''
def quick_split(pdb):
    pdb = pdb.lower()
    fake_pdb_container(pdb,filepos=os.path.join(pdb_PREFIX,pdb+'.pdb.gz'))

quick_split('1avd')
'''
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

