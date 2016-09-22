__author__= 'wy'

from prody import *
import random
import numpy as np
import os,re
import logging

'''
Core part for generating vectors
With the help of prody library.
Need installation first
Using:
  sudo pip install -r requirements.txt

'''

# Tag for hetero atoms
HETERO_PART =1
# Tag for protein atoms
PROTEIN_PART = 2

# score endurance with confidence
CONFIDENCE = 0.85

class pdb_container:
    '''
    For real pdb-ligand data
    It will separate each ligand (except ions)
    '''
    def __init__(self,PDB,filepos=None,OUT=True,**kwargs):
        self.PDBname= PDB
        self.heterodict = {}
        self.ct=0
        self.sequence = ''

        # filepos is to determine whether we download pdb files from wwPDB
        # or use what we have
        # Using downloaded is better
        try:
            if filepos is not None:
                parse = parsePDB(filepos)
            else:
                parse = parsePDB(PDB)
        except:
            #raise IOError
            logging.warning('PDB {} is ignored due to file-not-found error'.format(PDB))
            return

        #Generating sequence here
        storage = []
        for chain in parse.getHierView():
            for seq in storage:
                if chain.getSequence==seq:
                    continue
            self.sequence = self.sequence + repr(chain) +'|' + chain.getSequence()
            storage.append(chain.getSequence())

        # magic for selection desired atom group (see instruction of prody)
        hetero = parse.select('(hetero and not water) or resname ATP or resname ADP')
        other = parse.select('protein or nucleic')

        self.pure_protein = parse.select('protein')

        # bad one (no hetero or protein)
        if hetero is None or self.pure_protein is None or other is None:
            return
        if OUT:
            writePDB('data/{}_pure.pdb'.format(PDB),other)




        #Make vectors for every single hetero parts
        #Their values will be stored in a dict
        for pick_one in HierView(hetero).iterResidues():
            # less than 3 atoms may be not ok
            if pick_one.numAtoms() <= 3:
                continue

            ResId = str(pick_one.getResindex())

            # Extract this ligand from protein (as input for openbabel)
            if filepos is not None:
                filename = 'data/{}_{}.pdb'.format(filepos.split('/')[-1].split('.')[0], ResId)
            else:
                filename = 'data/{}_{}.pdb'.format(PDB, ResId)

            if not os.path.isfile(filename):
                if not os.path.exists('data'):
                    os.mkdir('data')

            if OUT and not os.path.exists(filename):
                writePDB(filename, pick_one)

            # Get coordinate of center
            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)
            #print middle

            scale = max(max(xyz[:, 0]) - middle[0], middle[0] - min(xyz[:, 0]),
                        max(xyz[:, 1]) - middle[1], middle[1] - min(xyz[:, 1]),
                        max(xyz[:, 2]) - middle[2], middle[2] - min(xyz[:, 2]))

            # assert scale <= 10
            if scale>10:
                logging.warning('Warning! {} has a ligand out of box scale with {} atom distance to center'.format(PDB, scale))
                #Now shifting the boxes:
                max_scale = max(max(xyz[:, 0]) - min(xyz[:, 0]),
                        max(xyz[:, 1])- min(xyz[:, 1]),
                        max(xyz[:, 2]) - min(xyz[:, 2]))
                if max_scale>20:
                    logging.error('Assertion failed, {} has a ligand out of box completely with scale'.format(PDB,scale))
                    continue
                #Try to move to the new center
                middle =  [(max(xyz[:, 0])+ min(xyz[:, 0]))/2,(max(xyz[:, 1])+min(xyz[:, 1]))/2,(max(xyz[:, 2])+min(xyz[:, 2]))/2]

            xx, yy, zz = np.meshgrid(np.linspace(middle[0] - 9.5, middle[0] + 9.5, 20),
                                     np.linspace(middle[1] - 9.5, middle[1] + 9.5, 20),
                                     np.linspace(middle[2] - 9.5, middle[2] + 9.5, 20))

            # print xx
            vector = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

            num_vector = [0] * 8000
            for atom in pick_one.iterAtoms():
                x, y, z = atom.getCoords()
                x_pos = int(round(x - vector[0][0]))
                #assert 0 <= x_pos <= 19
                y_pos = int(round(y - vector[0][1]))
                #assert 0 <= y_pos <= 19
                z_pos = int(round(z - vector[0][2]))
                #assert 0 <= z_pos <= 19
                if 0<= x_pos <= 19 and 0<=y_pos <= 19 and 0<=z_pos<=19:
                    #Simply change here to fulfill the mark as 'H_1'
                    num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getName()+'_'+str(HETERO_PART)

            # quick,dirty way to find atoms of protein in cubic boxes
            defSelectionMacro('inbox','abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0],middle[1],middle[2]))
            nearby= self.pure_protein.select('inbox')

            if nearby is not None:
                for atom in nearby.iterAtoms():
                    x, y, z = atom.getCoords()
                    x_pos = int(round(x - vector[0][0]))
                    #assert 0 <= x_pos <= 19
                    y_pos = int(round(y - vector[0][1]))
                    #assert 0 <= y_pos <= 19
                    z_pos = int(round(z - vector[0][2]))
                    #assert 0 <= z_pos <= 19
                    if 0<=x_pos<=19 and 0<=y_pos<=19 and 0<=z_pos<=19 and num_vector[x_pos * 400 + y_pos * 20 + z_pos]==0:
                        #Simply change here to fulfill the mark as 'C_2'
                        num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getName()+'_'+str(PROTEIN_PART)
                    else:
                        logging.warning('Coorinate {} {} {} found at {}'.format(x_pos,y_pos,z_pos,self.PDBname))

                #This is for checking the correctness when we add atoms in proteins.
                filename2= 'data/{}_{}_2.pdb'.format(PDB, ResId)
                writePDB(filename2, nearby)

            #Save into the dict for future locating
            self.heterodict[ResId] = {
                'raw_vector': num_vector,
                'center': middle,
                'selectmarco': 'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0],middle[1],middle[2]),
                'filename': filename,
                'id': ResId,
                'ligand': pick_one,
                'vina_score' : 0
            }
    def find_similar_target(self,sdf_filedir,**kwargs):
        '''
        Find the ligands that is highly possible to be the same compounds
        the default confidence is 0.85
        the score was based on tanimoto scoring method
        :param sdf_filedir: where the source sdf file is. In theory, if we are using openbabel
                            it is ok even if the file is not sdf, but it should only contain single
                            molecules, other wise this function cannot get right result
        :param kwargs:
        :return:
        '''

        assert isinstance(sdf_filedir,str)
        if not os.path.exists(sdf_filedir) or sdf_filedir.split('.')[-1]!='sdf':
            raise IOError('Please use a right location, {} is not a legal file name of sdf file'.format(sdf_filedir))

        possible_ones=[]

        for k,v in self.heterodict.items():
            try:
                command = os.popen('babel -d {0}/{1} {0}/{2} -ofpt -xfFP4'.format(os.getcwd(),sdf_filedir,v['filename']))
                ls= command.read()
                #print ls
                cp = re.split('=|\n', ls)[2]
                print cp
            except:
                with open('error.txt','a') as f:
                    f.write(self.PDBname+'\n')
                logging.warning('Babel encountered a problem at pdb {} ligand {}'.format(self.PDBname, v['filename']))
                cp = 0


            #print cp
            if float(cp) >= 0.85:
                v['cp']=float(cp)
                possible_ones.append(v)

        return possible_ones


    def pick_one(self, ResId, **kwargs):
        return self.heterodict[ResId] or None

    def list_ResId(self):
        return self.heterodict.keys()

    def set_vina_benchmark(self,ResId,BoxSize=20):
        '''
        Get one score from one ligang on one protein data.
        The score will be append to the corresponding dict
        :param ResId: The string that will be the key in heterodict
        :param BoxSize: The boxsize for docking, need consistent with the vector size
        :return: return the score and append it to the heterodict
        '''
        if self.heterodict.get(ResId,None) is None:
            print 'No specific one found!'
            return
        from Config import pythonsh_dir

        # write receptor pdbqt files:
        # receptor_file_loc = os.path.join('data/',self.PDBname+'_pure.pdb')
        receptor_file_loc = os.path.join('data/',self.PDBname+'_{}_2.pdb'.format(ResId))
        if not os.path.exists(receptor_file_loc+'qt'):
            os.system(os.path.join(pythonsh_dir, 'pythonsh') + ' prepare_receptor4.py -r {0} -o {0}qt'.format(receptor_file_loc))

        # write ligand pdbqt files
        ligand_file_loc = self.heterodict[ResId]['filename']
        os.system(os.path.join(pythonsh_dir,'pythonsh')+' prepare_ligand4.py -l {0} -o {0}qt'.format(ligand_file_loc))

        # write config files
        with open('vina_config.txt','w') as f:
            f.write('    receptor = {}qt\n'.format(receptor_file_loc))
            f.write('    ligand = {}qt\n'.format(ligand_file_loc))
            f.write('    center_x = {}\n'.format(self.heterodict[ResId]['center'][0]))
            f.write('    center_y = {}\n'.format(self.heterodict[ResId]['center'][1]))
            f.write('    center_z = {}\n'.format(self.heterodict[ResId]['center'][2]))
            f.write('    size_x = {}\n'.format(BoxSize))
            f.write('    size_y = {}\n'.format(BoxSize))
            f.write('    size_z = {}\n'.format(BoxSize))
            f.close()

        #Now do docking:
        #Suppose vina is installed
        command = os.popen('vina --config vina_config.txt --score_only')

        ls = command.read()
        for line in ls.split('\n'):
            if 'Affinity' in line:
                real_num = re.compile(r"[-+]?\d+\.\d+")
                score = real_num.search(line.split(':')[1])
                if score:
                    self.heterodict[ResId]['vina_score'] =float(score.group())
                    return float(score.group())
                else:
                    logging.error('Cannot give a score at {} and {}'.format(self.PDBname,self.heterodict[ResId]['ligand'].getResname()))
                    with open('error.txt','a') as f:
                        f.write(self.PDBname+' '+self.heterodict[ResId]['ligand'].getResname()+'\n')
                    return 0;

        #scp =  re.split('=|\n', ls)[2]

    def set_all_vina_benchmark(self,Box=20):
        '''
        Just get all vina_score generated
        :param BoxSize: box size for docking
        :return:
        '''
        for k in self.heterodict.keys():
            self.set_vina_benchmark(k,BoxSize=Box)

    '''
    def __iter__(self):
        if self.heterodict is None:
            return {}
        return self.heterodict
    '''

    def __repr__(self):
        print self.PDBname+'({} hetero parts found)'.format(len(self.heterodict.keys()))


class fake_pdb_container:
    '''
    For fake docking results.
    Proteins do not have info about ligands
    Need append manually
    '''
    def __init__(self,PDB,filepos=None):
        self.PDBname = PDB
        self.heterodict = {}
        self.ct=0
        self.sequence = ''

        # filepos is to determine whether we download pdb files from wwPDB
        # or use what we have
        # Using downloaded is better
        try:
            if filepos is not None:
                parse = parsePDB(filepos)
            else:
                parse = parsePDB(PDB)
        except:
            # raise IOError
            logging.warning('PDB {} is ignored due to file-not-found error'.format(PDB))
            return

        self.protein = parse.select('not water')
        print self.protein.numAtoms()

        # Generating sequence here
        storage = []
        for chain in parse.getHierView():
            for seq in storage:
                if chain.getSequence == seq:
                    continue
            self.sequence = self.sequence + repr(chain) + '|' + chain.getSequence()
            storage.append(chain.getSequence())


    def self_generating(self):
        from Config import fake_result_PREFIX,NAME,key
        import csv
        filename= 'fake_{}.csv'.format(self.PDBname)
        writer = file(os.path.join(fake_result_PREFIX+filename), 'wb')
        w = csv.writer(writer)
        w.writerow(['Name', NAME, 'Target PDB', 'ResIndex', 'Similarity'] + key + ['Vector','Sequence'])

        for k,v in self.heterodict.items():
            line = ['sth', 'NA', self.PDBname, 'NA', 'NA' ] + ['NA']*len(key) + [v['raw_vector'],self.sequence]
            w.writerow(line)
        writer.flush()
        writer.close()

    def append_vectors(self,hetero_file):
        '''
        Append each docked result as a vector to the dict
        :param hetero_file: file position
        :return: nothing , but will generate a vector into the dict
        '''

        #need to split the files
        TEMP = 'temp.pdb'
        o = open(hetero_file,'r')
        one_pdb=''
        for line in o:
            one_pdb+=line
            if 'END' in line:
                #write a temporial file
                with open(TEMP,'wb') as w:
                    w.write(one_pdb)
                    w.close()
                one_pdb = ''
                pdb = parsePDB(TEMP)
                if pdb.numAtoms() <= 3:
                    continue


                # Get coordinate of center
                xyz = pdb.getCoords()
                middle = calcCenter(pdb)

                scale = max(max(xyz[:, 0]) - middle[0], middle[0] - min(xyz[:, 0]),
                            max(xyz[:, 1]) - middle[1], middle[1] - min(xyz[:, 1]),
                            max(xyz[:, 2]) - middle[2], middle[2] - min(xyz[:, 2]))

                # assert scale <= 10
                if scale > 10:
                    logging.warning(
                        'Warning! {} has a ligand out of box scale with {} atom distance to center'.format(self.PDBname,scale))
                    # Now shifting the boxes:
                    max_scale = max(max(xyz[:, 0]) - min(xyz[:, 0]),
                                    max(xyz[:, 1]) - min(xyz[:, 1]),
                                    max(xyz[:, 2]) - min(xyz[:, 2]))
                    if max_scale > 20:
                        logging.error(
                            'Assertion failed, {} has a ligand out of box completely with scale'.format(self.PDBname, scale))
                        continue
                    # Try to move to the new center
                    middle = [(max(xyz[:, 0]) + min(xyz[:, 0])) / 2, (max(xyz[:, 1]) + min(xyz[:, 1])) / 2,
                              (max(xyz[:, 2]) + min(xyz[:, 2])) / 2]

                xx, yy, zz = np.meshgrid(np.linspace(middle[0] - 9.5, middle[0] + 9.5, 20),
                                         np.linspace(middle[1] - 9.5, middle[1] + 9.5, 20),
                                         np.linspace(middle[2] - 9.5, middle[2] + 9.5, 20))

                # print xx
                vector = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

                num_vector = [0] * 8000
                for atom in pdb.iterAtoms():
                    x, y, z = atom.getCoords()
                    x_pos = int(round(x - vector[0][0]))
                    # assert 0 <= x_pos <= 19
                    y_pos = int(round(y - vector[0][1]))
                    # assert 0 <= y_pos <= 19
                    z_pos = int(round(z - vector[0][2]))
                    # assert 0 <= z_pos <= 19
                    if 0 <= x_pos <= 19 and 0 <= y_pos <= 19 and 0 <= z_pos <= 19:
                        # Simply change here to fulfill the mark as 'H_1'
                        num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getName()+ '_' + str(HETERO_PART)

                # quick,dirty way to find atoms of protein in cubic boxes
                defSelectionMacro('inbox',
                                  'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0], middle[1],
                                                                                                middle[2]))
                nearby = self.protein.select('inbox')

                if nearby is not None:
                    for atom in nearby.iterAtoms():
                        x, y, z = atom.getCoords()
                        x_pos = int(round(x - vector[0][0]))
                        # assert 0 <= x_pos <= 19
                        y_pos = int(round(y - vector[0][1]))
                        # assert 0 <= y_pos <= 19
                        z_pos = int(round(z - vector[0][2]))
                        # assert 0 <= z_pos <= 19
                        if 0 <= x_pos <= 19 and 0 <= y_pos <= 19 and 0 <= z_pos <= 19 and num_vector[
                                                    x_pos * 400 + y_pos * 20 + z_pos] == 0:
                            # Simply change here to fulfill the mark as 'C_2'
                            num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getName() + '_' + str(PROTEIN_PART)
                        else:
                            logging.warning('Coorinate {} {} {} found at {}'.format(x_pos, y_pos, z_pos, self.PDBname))

                            # This is for checking the correctness when we add atoms in proteins.
                            # filename2= 'data/{}_{}_2.pdb'.format(PDB, ResId)
                            # writePDB(filename2, pick_one+nearby)

                # Save into the dict for future locating
                self.heterodict[str(self.ct)] = {
                    'raw_vector': num_vector,
                    'center': middle,
                    'filename': hetero_file,
                    'id': hetero_file.split('/')[-1].split('.')[0]+'_'+str(self.ct)
                }
                self.ct+=1

