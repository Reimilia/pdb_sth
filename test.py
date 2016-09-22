from fileparser import do_one_pdb
from Config import pdb_PREFIX,fake_src_PREFIX,fake_pdb_name,fake_hetero_PREFIX
import os,sys,io
import gzip
from vector_gen import pdb_container,fake_pdb_container
from Config import pdb_PREFIX

PDB_tar=['1uwj']

#do_one_pdb('1j8q')
#do_one_pdb('1vq6')
A= pdb_container('1j8q',filepos=os.path.join(pdb_PREFIX,'1j8q.pdb.gz'))
print A.set_vina_benchmark('147')
'''
a= fake_pdb_container('aa2ar',filepos=fake_src_PREFIX+'aa2ar/'+fake_pdb_name)

filenames =  os.listdir(fake_hetero_PREFIX)

for filename in filenames:
    if filename.split('.')[-1]=='pdb':
        a.append_vectors(os.path.join(fake_hetero_PREFIX,filename))

a.self_generating()
'''