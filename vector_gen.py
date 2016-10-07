__author__= 'wy'

from prody import *
import random
import numpy as np
import os,re
import logging
from mapping import *

'''
Core part for generating vectors and split source pdb files with
With the help of prody library.
Need installation first
Using:
  sudo pip install -r requirements.txt

For local installation, see: http://prody.csb.pitt.edu/downloads/

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
    def get_pdb_type(self):
        '''
        Nucleic Protein or Complex
        :return:
        '''
        if self.pure_protein is not None:
            if self.pure_nucleic is None:
                return 'Protein'
            else:
                return 'Protein_Nucleic_Complex'
        else:
            if self.pure_protein is None:
                return 'Nucleic'
            else:
                return 'Unknown or empty'

    def __init__(self,PDB,filepos=None,OUT=True,**kwargs):
        '''

        :param PDB: name of PDB
        :param filepos: directory of where PDB file stores
        :param OUT: if true, splitted files will be output in './data' folder
        :param kwargs: for further extension
        '''
        self.PDBname= PDB
        self.heterodict = {}
        self.ct=0
        self.sequence = ''
        self.pure_protein= None
        self.pure_nucleic= None

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
            #print chain
            for seq in storage:
                if chain.getSequence()==seq:
                    continue
            self.sequence = self.sequence + repr(chain) +'|' + chain.getSequence()
            storage.append(chain.getSequence())

        # magic for selection desired atom group (see instruction of prody)
        # parse = parse.select('not hydrogen')

        hetero = parse.select('(hetero and not water) or resname ATP or resname ADP')
        other = parse.select('protein or nucleic and not (resname ATP or resname ADP)')
        #print parse.numAtoms(), hetero.numAtoms(), other.numAtoms()

        self.pure_protein = parse.select('protein')
        self.pure_nucleic = parse.select('nucleic')
        self.other = other

        # bad one (no hetero or protein)
        if hetero is None or other is None:
            return
        if OUT:
            writePDB('data/{}_receptor.pdb'.format(PDB),other)
            writePDB('data/{}_complex.pdb'.format(PDB),parse)




        #Make vectors for every single hetero parts
        #Their values will be stored in a dict
        for pick_one in HierView(hetero).iterResidues():
            # less than 3 atoms may be not ok
            if pick_one.numAtoms() <= 3:
                continue

            ResId = str(pick_one.getResindex())

            # Extract this ligand from protein (as input for openbabel)
            if filepos is not None:
                filename = 'data/{}_{}_ligand.pdb'.format(filepos.split('/')[-1].split('.')[0], ResId)
            else:
                filename = 'data/{}_{}_ligand.pdb'.format(PDB, ResId)

            if not os.path.isfile(filename):
                if not os.path.exists('data'):
                    os.mkdir('data')

            if OUT and not os.path.exists(filename):
                writePDB(filename, pick_one)

            # Get coordinate of center
            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)
            # in pi degree , the rotation of the box (if needed)
            rotation= [0,0,0]

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

            #print middle
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

            # This place might have some potential problem
            # for ADP or ATP , they might either be part of nucleic and the ligand
            # This will cause a severe bug when calculating autovina score
            # TODO fix this issue
            nearby= other.select('inbox')

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
                        print atom.getName()
                        logging.warning('Coorinate {} {} {} found at {}'.format(x_pos,y_pos,z_pos,self.PDBname))

                if OUT:
                    #Output the pure protein part in the box and the ligand-protein complex part
                    filename2= 'data/{}_{}_receptor.pdb'.format(PDB, ResId)
                    writePDB(filename2, nearby)
                    filename2= 'data/{}_{}_complex.pdb'.format(PDB,ResId)
                    writePDB(filename2, nearby+pick_one)

            #Save into the dict for future locating

            # Do autogrid part:
            fake_ligand_filename = os.path.join(temp_pdb_PREFIX, 'fake-ligand.pdb')
            naming = '{}_{}'.format(PDB,ResId)
            ligand_filename = os.path.join(temp_pdb_PREFIX, naming + '_ligand.pdb')
            receptor_filename = os.path.join(temp_pdb_PREFIX, naming + '_receptor.pdb')
            complex_filename = os.path.join(temp_pdb_PREFIX, naming + '_complex.pdb')
            do_auto_grid(receptor_filename, fake_ligand_filename, center=middle)
            do_auto_grid(ligand_filename, fake_ligand_filename, center=middle)
            do_auto_grid(complex_filename, fake_ligand_filename, center=middle)

            self.heterodict[ResId] = {
                'raw_vector': num_vector,
                'center': middle,
                'rotation': rotation,
                'selectmarco': 'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0],middle[1],middle[2]),
                'naming': '{}_{}'.format(PDB,ResId),
                'filename': filename,
                'id': ResId,
                'ligand': pick_one,
                'vina_score' : 'NA',
                'gridmap_protein': fetch_gridmaps(naming+'_receptor'),
                'gridmap_ligand': fetch_gridmaps(naming+'_ligand'),
                'gridmap_complex': fetch_gridmaps(naming+'_complex')
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

    def add_ligand(self,ligand_pdb_file,index,OUT=True):
        '''
        Add ligands on to pdb. The result should be generated by docking
        :param ligand_pdb_file:
        :return:
        '''
        try:
            parse = parsePDB(ligand_pdb_file)

        except:
            #raise IOError
            logging.warning('cannot add ligang file on PDB {}'.format(self.PDBname))
            return

        hetero = parse.select('hetero')
        if hetero is None:
            logging.error('no ligands were found!')
            return

        # Make vectors for every single hetero parts
        # Their values will be stored in a dict
        for pick_one in HierView(hetero).iterResidues():
            # less than 3 atoms may be not ok
            if pick_one.numAtoms() <= 3:
                continue

            ResId = str(pick_one.getResindex())

            # Extract this ligand from protein (as input for openbabel)
            filename = 'data/{}_{}_ligand.pdb'.format(filepos.split('/')[-1].split('.')[0], ResId)
            else:
                filename = 'data/{}_{}_ligand.pdb'.format(PDB, ResId)

            if not os.path.isfile(filename):
                if not os.path.exists('data'):
                    os.mkdir('data')

            if OUT and not os.path.exists(filename):
                writePDB(filename, pick_one)

            # Get coordinate of center
            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)
            # in pi degree , the rotation of the box (if needed)
            rotation = [0, 0, 0]

            scale = max(max(xyz[:, 0]) - middle[0], middle[0] - min(xyz[:, 0]),
                        max(xyz[:, 1]) - middle[1], middle[1] - min(xyz[:, 1]),
                        max(xyz[:, 2]) - middle[2], middle[2] - min(xyz[:, 2]))

            # assert scale <= 10
            if scale > 10:
                logging.warning(
                    'Warning! {} has a ligand out of box scale with {} atom distance to center'.format(PDB, scale))
                # Now shifting the boxes:
                max_scale = max(max(xyz[:, 0]) - min(xyz[:, 0]),
                                max(xyz[:, 1]) - min(xyz[:, 1]),
                                max(xyz[:, 2]) - min(xyz[:, 2]))
                if max_scale > 20:
                    logging.error(
                        'Assertion failed, {} has a ligand out of box completely with scale'.format(PDB, scale))
                    continue
                # Try to move to the new center
                middle = [(max(xyz[:, 0]) + min(xyz[:, 0])) / 2, (max(xyz[:, 1]) + min(xyz[:, 1])) / 2,
                          (max(xyz[:, 2]) + min(xyz[:, 2])) / 2]

            # print middle
            xx, yy, zz = np.meshgrid(np.linspace(middle[0] - 9.5, middle[0] + 9.5, 20),
                                     np.linspace(middle[1] - 9.5, middle[1] + 9.5, 20),
                                     np.linspace(middle[2] - 9.5, middle[2] + 9.5, 20))

            # print xx
            vector = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

            num_vector = [0] * 8000
            for atom in pick_one.iterAtoms():
                x, y, z = atom.getCoords()
                x_pos = int(round(x - vector[0][0]))
                # assert 0 <= x_pos <= 19
                y_pos = int(round(y - vector[0][1]))
                # assert 0 <= y_pos <= 19
                z_pos = int(round(z - vector[0][2]))
                # assert 0 <= z_pos <= 19
                if 0 <= x_pos <= 19 and 0 <= y_pos <= 19 and 0 <= z_pos <= 19:
                    # Simply change here to fulfill the mark as 'H_1'
                    num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getName() + '_' + str(HETERO_PART)

            # quick,dirty way to find atoms of protein in cubic boxes
            defSelectionMacro('inbox',
                              'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0], middle[1],
                                                                                            middle[2]))

            # This place might have some potential problem
            # for ADP or ATP , they might either be part of nucleic and the ligand
            # This will cause a severe bug when calculating autovina score
            # TODO fix this issue
            nearby = other.select('inbox')

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
                        print atom.getName()
                        logging.warning('Coorinate {} {} {} found at {}'.format(x_pos, y_pos, z_pos, self.PDBname))

                if OUT:
                    # Output the pure protein part in the box and the ligand-protein complex part
                    filename2 = 'data/{}_{}_receptor.pdb'.format(PDB, ResId)
                    writePDB(filename2, nearby)
                    filename2 = 'data/{}_{}_complex.pdb'.format(PDB, ResId)
                    writePDB(filename2, nearby + pick_one)

            # Save into the dict for future locating

            # Do autogrid part:
            fake_ligand_filename = os.path.join(temp_pdb_PREFIX, 'fake-ligand.pdb')
            naming = '{}_{}'.format(PDB, ResId)
            ligand_filename = os.path.join(temp_pdb_PREFIX, naming + '_ligand.pdb')
            receptor_filename = os.path.join(temp_pdb_PREFIX, naming + '_receptor.pdb')
            complex_filename = os.path.join(temp_pdb_PREFIX, naming + '_complex.pdb')
            do_auto_grid(receptor_filename, fake_ligand_filename, center=middle)
            do_auto_grid(ligand_filename, fake_ligand_filename, center=middle)
            do_auto_grid(complex_filename, fake_ligand_filename, center=middle)

            self.heterodict[ResId] = {
                'raw_vector': num_vector,
                'center': middle,
                'rotation': rotation,
                'selectmarco': 'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0], middle[1],
                                                                                             middle[2]),
                'naming': '{}_{}'.format(PDB, ResId),
                'filename': filename,
                'id': ResId,
                'ligand': pick_one,
                'vina_score': 'NA',
                'gridmap_protein': fetch_gridmaps(naming + '_receptor'),
                'gridmap_ligand': fetch_gridmaps(naming + '_ligand'),
                'gridmap_complex': fetch_gridmaps(naming + '_complex')
            }



    def add_ligands(self,ligand_file):
        try:
            filename = ligand_file.split('/')[-1]
            with open(ligand_file,'rb') as f:
                output=''
                index=0
                for line in f :
                    output+=line
                    if 'ENDMOL' in line:
                        with open('data/{}'.format(filename),'wb') as w:
                            w.write(output)
                        self.add_ligand(self,'data/'+filename, index)
                        output = ''
                        index+=1


            return True
        except:
            return False

    def __repr__(self):
        print self.PDBname+'({} hetero parts found)'.format(len(self.heterodict.keys()))

    def clean_temp_file(self):
        pass



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

