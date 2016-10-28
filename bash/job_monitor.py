import os,sys
import time

def get_file_list():
    filename='/home/yw174/pdb_data/input.txt'
    with open(filename,'rb') as f:
        list_dir =[]
        for each in f.readlines():
            list_dir.append(each)
    return list_dir


def submit_one_job(job_suffix,jobtype,pdb,resid):

    print jobtype,pdb,resid
    string = jobtype+ '_' + pdb+ '_' +resid
    print string
    with open('dock.sh', 'w') as w:
        w.write('# !/bin/bash\n')
        w.write('# BSUB -n 1\n')
        w.write('# BSUB -W 5:00\n')
        w.write('# BSUB -J dock_%s\n'%(string))
        w.write('# BSUB -o /home/yw174/job/dock_decorated/%s/out_%s\n'%(jobtype,string))
        w.write('# BSUB -e /home/yw174/job/dock_decorated/%s/err_%s\n'%(jobtype,string))
        w.write('# BSUB -q short\n')
        w.write('export OMP_NUM_THREADS=1\n')
        w.write('export PATH=$PATH:/home/yw174/usr/babel/bin/\n')
        w.write('source /home/yw174/python_env/wy/bin/activate\n')
        w.write('cd /home/yw174/program/pdb_sth\n')
        cmd = 'python job_dispatcher.py %s %s %s'%(jobtype,pdb,resid)
        w.write(cmd + '\n')
    os.system('chmod 777 dock.sh')
    os.system('bsub < dock.sh')
    os.remove('dock.sh')

if __name__=='__main__':
    index = 0
    filelist = get_file_list()
    while True:
        time.sleep(60)
        command=os.popen('bjobs | grep short | wc -l')
        ls = int(command.read())
        while ls<100:
            filename=filelist[index]
            index+=1
            pdb = filename.split('_')[0]
            resid = filename.split('_')[1]
            submit_one_job('fast','fast',pdb, resid)
            submit_one_job('rigor','rigor',pdb, resid)
            submit_one_job('rigor_so', 'rigor_so', pdb, resid)
            ls+=3


