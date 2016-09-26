'''
Writing external shells here.
This time the factory mode is changed to file-initiated.
I.e. all functions receive input as files
'''

import sys,io,os
from functools import wraps
import time
import re
from Autodock_Config import autodock_store_dir,pythonsh_dir

WORK_DIR = os.getcwd()
CURRENT_DIR = os.getcwd()+'/mapping'
#os.chdir(CURRENT_DIR)

def set_new_folder(PDBname,storedir):
    #os.chdir(storedir)
    if not os.path.exists(os.path.join(storedir,PDBname)):
        os.mkdir(os.path.join(storedir,PDBname))
    #os.chdir(os.getcwd())

def fn_timer(function):
    '''
    This is the decorator used for time counting issue
    Need not understand this one. It has nothing to do with generating files
    :param function:
    :return: no return. just print and record the time the decorated program ran.
    '''
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" %
               (function.func_name, str(t1 - t0))
               )
        return result

    return function_timer


def prepare_receptor(filename,pdbname):
    if filename.split('.')[-1]!='pdb':
        print 'Error! when prepare receptor'
        return
    real_dir = os.path.join(autodock_store_dir,pdbname)
    real_filepos= os.path.join(real_dir,filename.split('/')[-1])+'qt'
    if not os.path.exists(real_filepos):
        os.chdir(CURRENT_DIR)
        os.system(os.path.join(pythonsh_dir, 'pythonsh') + ' prepare_receptor4.py -r {0} -o {1}'
              .format(filename,real_filepos))
        os.chdir(WORK_DIR)
    #print 'Ok'


def prepare_ligand(filename,pdbname):

    if filename.split('.')[-1]!='pdb':
        print 'Error! when prepare ligand'
        return
    real_dir = os.path.join(autodock_store_dir, pdbname)
    real_filepos = os.path.join(real_dir, filename.split('/')[-1]) + 'qt'
    if not os.path.exists(real_filepos):
        os.chdir(CURRENT_DIR)
        os.system(os.path.join(pythonsh_dir, 'pythonsh') + ' prepare_ligand4.py -l {0} -o {1}'
                  .format(filename, real_filepos))
        os.chdir(WORK_DIR)
    #print 'Ok'

@fn_timer
def do_auto_grid(receptor,ligand,center=None):
    rname = receptor.split('/')[-1]
    lname = ligand.split('/')[-1]
    pdbname = rname.split('_')[0]


    if not os.path.exists(receptor) or not os.path.exists(ligand):
        print receptor
        print ligand
        raise TypeError('fileformat')
        return
    if rname.split('.')[-1]=='pdb':
        prepare_receptor(receptor,pdbname)
        rname+='qt'
    else:
        if rname.split('.')[-1]!='pdbqt':
            raise TypeError('fileformat')
            return

    if lname.split('.')[-1] == 'pdb':
        prepare_ligand(ligand,pdbname)
        lname+='qt'
    else:
        if lname.split('.')[-1] != 'pdbqt':
            raise TypeError('fileformat')
            return

    naming = "".join(rname.split('.')[:-1])
    real_dir = os.path.join(autodock_store_dir,pdbname)
    glg_output_dir = os.path.join(real_dir,naming)

    rloc = os.path.join(real_dir,rname)
    lloc = os.path.join(real_dir,lname)

    os.chdir(CURRENT_DIR)

    if center is None:
        os.system(os.path.join(pythonsh_dir, 'pythonsh') + ' prepare_gpf4.py -l {} -r {} -o {}.gpf -p spacing=1.0 -p npts=\"20,20,20\" '
              .format(lloc,rloc,glg_output_dir))
    else:
        os.system(os.path.join(pythonsh_dir,'pythonsh') +
                  ' prepare_gpf4.py -l {} -r {} -o {}.gpf -p spacing=1.0 -p npts=\"20,20,20\" -p gridcenter=\"{},{},{}\" '
                  .format(lloc,rloc ,glg_output_dir, center[0],center[1],center[2]))

    #Suppose autogrid and autodock has installed
    os.chdir(real_dir)
    os.system('autogrid4 -p {0}.gpf -l {0}.glg'.format(glg_output_dir))
    os.chdir(WORK_DIR)

    #print 'Ok'


@fn_timer
def do_auto_dock(receptor,ligand,center=None):
    rname = receptor.split('/')[-1]
    lname = ligand.split('/')[-1]
    pdbname = rname.split('_')[0]

    if not os.path.exists(receptor) or not os.path.exists(ligand):
        raise TypeError('fileformat')
        return

    do_auto_grid(receptor,ligand,center)

    if rname.split('.')[-1] == 'pdb':
        prepare_receptor(receptor,pdbname)
        rname += 'qt'
    else:
        if rname.split('.')[-1] != 'pdbqt':
            raise TypeError('fileformat')
            return

    if lname.split('.')[-1] == 'pdb':
        prepare_ligand(ligand,pdbname)
        lname += 'qt'
    else:
        if lname.split('.')[-1] != 'pdbqt':
            raise TypeError('fileformat')
            return

    os.chdir(CURRENT_DIR)

    naming = "".join(rname.split('.')[:-1])
    real_dir = os.path.join(autodock_store_dir, pdbname)
    dlg_output_dir = os.path.join(real_dir, naming)

    rloc = os.path.join(real_dir, rname)
    lloc = os.path.join(real_dir, lname)

    os.system(os.path.join(pythonsh_dir,'pythonsh') + ' prepare_dpf42.py -l {} -r {} -o {}.dpf'
              .format(lloc, rloc, dlg_output_dir ))

    # Suppose autogrid and autodock has installed
    os.chdir(real_dir)
    #os.system(os.path.join('autodock4 -p {0}.dpf -l {0}.dlg'.format(dlg_output_dir)))
    os.chdir(WORK_DIR)

    print 'Ok'


@fn_timer
def do_auto_vina_score(receptor,ligand,center,Box=20):
    # receptor_file_loc = os.path.join('data/',self.PDBname+'_{}_2.pdb'.format(ResId))
    rname = receptor.split('/')[-1]
    lname = ligand.split('/')[-1]
    pdbname = rname.split('_')[0]

    if not os.path.exists(receptor) or not os.path.exists(ligand):
        raise TypeError('fileformat')
        return
    if rname.split('.')[-1]=='pdb':
        prepare_receptor(receptor,pdbname)
        rname+='qt'
    else:
        if rname.split('.')[-1]!='pdbqt':
            raise TypeError('fileformat')
            return

    if lname.split('.')[-1] == 'pdb':
        prepare_ligand(ligand,pdbname)
        lname+='qt'
    else:
        if lname.split('.')[-1] != 'pdbqt':
            raise TypeError('fileformat')
            return

    real_dir = os.path.join(autodock_store_dir, pdbname)

    rloc = os.path.join(real_dir, rname)
    lloc = os.path.join(real_dir, lname)

    os.chdir(real_dir)
    # write config files
    with open('vina_config.txt', 'w') as f:
        f.write('    receptor = {}\n'.format(rname))
        f.write('    ligand = {}\n'.format(lname))
        f.write('    center_x = {}\n'.format(center[0]))
        f.write('    center_y = {}\n'.format(center[1]))
        f.write('    center_z = {}\n'.format(center[2]))
        f.write('    size_x = {}\n'.format(Box))
        f.write('    size_y = {}\n'.format(Box))
        f.write('    size_z = {}\n'.format(Box))
        f.close()

    # Now do docking:
    # Suppose vina is installed
    command = os.popen('vina --config vina_config.txt --score_only')

    os.chdir(WORK_DIR)

    ls = command.read()
    for line in ls.split('\n'):
        if 'Affinity' in line:
            real_num = re.compile(r"[-+]?\d+\.\d+")
            score = real_num.search(line.split(':')[1])
            if score:
                return float(score.group())
            else:
                return 'NA'

                # scp =  re.split('=|\n', ls)[2]

if __name__=='__main__':
    #Example on how to finish auto docking process

    set_new_folder('1j8q','/home/wy/Documents/BCH_coding/pdb_data_extracter/result')

    #protein only
    do_auto_dock('/home/wy/Documents/BCH_coding/pdb_data_extracter/data/1j8q_147_pure.pdb',
                 '/home/wy/Documents/BCH_coding/pdb_data_extracter/data/fake-ligand.pdb',center=[21.36,10.47,81.86])

    #ligand only
    do_auto_dock('/home/wy/Documents/BCH_coding/pdb_data_extracter/data/1j8q_147_ligand.pdb',
                 '/home/wy/Documents/BCH_coding/pdb_data_extracter/data/fake-ligand.pdb', center=[21.36, 10.47, 81.86])

    #protein-ligand complex
    do_auto_dock('/home/wy/Documents/BCH_coding/pdb_data_extracter/data/1j8q_147_complex.pdb',
                 '/home/wy/Documents/BCH_coding/pdb_data_extracter/data/fake-ligand.pdb', center=[21.36, 10.47, 81.86])