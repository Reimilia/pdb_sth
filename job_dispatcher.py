'''
This is the class for docking result processing
JOBname
BASEFile dir
Benchmarkfile dir


'''
from Config import pdb_PREFIX,temp_pdb_PREFIX,result_PREFIX
from vector_gen import pdb_container
import os

class dock_dispatcher:

    def __init__(self,**kwargs):

        if 'jobname' in kwargs:
            self.JOBname = kwargs['jobname']
        else:
            self.JOBname = 'Unknown'

        if 'filedir' in kwargs:
            self.file_dir = kwargs['filedir']
        else:
            print 'No specific docking result root is found!'
            return

        if 'benchmark' in kwargs:
            self.benchmark_dir = kwargs['benchmark']
        else:
            self.benchmark_dir = None

    def do_one_ligand(self,pdb_name,ligand_name,pdb_filepos=None,ligand_filepos=None):
        if pdb_filepos is None:
            pdb_filepos = os.path.join(pdb_PREFIX,pdb_name+'.pdb.gz')

        A = pdb_container(pdb_name, filepos=pdb_filepos)
        partial_name = pdb_name+'/'+pdb_name+'_'+ligand_name+'_ligand'

        if ligand_filepos is None:
            ligand_filepos = os.path.join(self.file_dir, partial_name+'.mol2')
        if self.benchmark_dir is not None:
            benchmark_file = os.path.join(self.benchmark_dir, partial_name+'.pdb')
        else:
            benchmark_file = None
        result_dir= os.path.join(result_PREFIX,self.JOBname)
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        A.add_ligands(ligand_filepos, suffix=self.JOBname, benchmark_file=benchmark_file)
