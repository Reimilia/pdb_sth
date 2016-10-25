'''
This script is only intended for Orchestra job_dispatch to do autodock_vina for each ligand.
'''


from Config import pdb_PREFIX,temp_pdb_PREFIX,result_PREFIX
from vector_gen import pdb_container
import os,sys
import thread

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
        try:
            A.add_ligands(ligand_filepos, suffix=self.JOBname, benchmark_file=benchmark_file)
            print 'Done %s_%s' %(pdb_name,ligand_name)
        except:
            print 'Can\'t find files or something is wrong! at %s_%s'%(pdb_name,ligand_name)

benchmark_dir = '/n/scratch2/xl198/data/H/addH'

job_dir = {
    'fast' : '/n/scratch2/xl198/data/H/wp_fast',
    'rigor': '/n/scratch2/xl198/data/H/wp_rigorous',
    'rigor_so': '/n/scratch2/xl198/data/H/so_rigorous',
    'random': ''
}
job_dispatcher = {
    'fast' : dock_dispatcher(jobname='fast',filedir= job_dir['fast'],benchmark= benchmark_dir),
    'rigor': dock_dispatcher(jobname='rigor',filedir= job_dir['rigor'],benchmark= benchmark_dir),
    'rigor_so': dock_dispatcher(jobname='rigor_so',filedir= job_dir['rigor_so'],benchmark= benchmark_dir),
    'random': dock_dispatcher(jobname='random',filedir= job_dir['random'],benchmark= benchmark_dir),
}



if __name__ == '__main__':
    try:
        job_type = sys.argv[1]
        pdb = sys.argv[2]
        resid = sys.argv[3]
    except:
        print 'Argument expected to be three!'
    else:
        if job_type not in job_dispatcher:
            print 'Incorrect job type!'
        else:
            print 'Submitting job {} to {}_{}'.format(job_type,pdb,resid)
            thread.sleep(5)
            job_dispatcher[job_type].do_one_ligand(pdb,resid)

    '''
    path= os.path.join(result_PREFIX,'experiment')
    fast_dir = '/n/scratch2/xl198/data/H/wp_fast'
    benchmark_dir = '/n/scratch2/xl198/data/H/addH'
    fast_job = dock_dispatcher(jobname='fast',filedir= fast_dir,benchmark= benchmark_dir)
    filenames = os.listdir(path)
    for filename in filenames:
        pdb = filename.split('_')[0]
        resid= filename.split('_')[1]
        fast_job.do_one_ligand(pdb,resid)
    '''
