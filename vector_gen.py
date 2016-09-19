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
    def __init__(self,PDB,filepos=None,**kwargs):
        self.PDBname= PDB
        self.heterodict = {}

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

        # magic for selection desired atom group (see instruction of prody)
        hetero = parse.select('not protein and not nucleic and not water')
        other = parse.select('protein or nucleic')

        print hetero.numAtoms()
        print other.numAtoms()

        # bad one (no hetero or protein)
        if hetero is None or other is None:
            return

        #Make vectors for every single hetero parts
        #Their values will be stored in a dict
        for pick_one in HierView(hetero).iterResidues():
            # less than 3 atoms may be not ok
            print 'here'
            if pick_one.numAtoms() <= 3:
                continue

            ResId = str(pick_one.getResindex())

            # Extract this ligand from protein (as input for openbabel)

            filename = 'data/{}_{}.pdb'.format(PDB, ResId)

            if not os.path.isfile(filename):
                if not os.path.exists('data'):
                    os.mkdir('data')

            if not os.path.exists(filename):
                writePDB(filename, pick_one)

            # Get coordinate of center
            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)

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
                    num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getElement()+'_'+str(HETERO_PART)

            # quick,dirty way to find atoms of protein in cubic boxes
            defSelectionMacro('inbox','abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0],middle[1],middle[2]))
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
                    if 0<=x_pos<=19 and 0<=y_pos<=19 and 0<=z_pos<=19:
                        #Simply change here to fulfill the mark as 'C_2'
                        num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getElement()+'_'+str(PROTEIN_PART)
                    else:
                        logging.warning('Coorinate {} {} {} found'.format(x_pos,y_pos,z_pos))

                #This is for checking the correctness when we add atoms in proteins.
                #filename2= 'data/{}_{}_2.pdb'.format(PDB, ResId)
                #writePDB(filename2, pick_one+nearby)

            #Save into the dict for future locating
            self.heterodict[ResId] = {
                'raw_vector': num_vector,
                'center': middle,
                'selectmarco': 'abs(x-{}) < 10 and abs(y-{}) < 10 and abs(z-{}) < 10'.format(middle[0],middle[1],middle[2]),
                'filename': filename,
                'id': pick_one.getResindex()
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
                command = os.popen('babel {0}/{1} {0}/{2} -ofpt'.format(os.getcwd(),sdf_filedir,v['filename']))
                cp = re.split('=|\n', command.read())[2]
            except:
                logging.warning('Babel encountered a problem at pdb {} ligand {}'.format(self.PDBname, v['filename']))
                cp = 0
                raise TypeError

            if float(cp) >= 0.85:
                v['cp']=float(cp)
                possible_ones.append(v)

        return possible_ones

    def self_generating(self):
        from Config import result_PREFIX,NAME,key
        import csv
        filename= 'fake_{}.csv'.format(self.PDBname)
        writer = file(os.path.join(result_PREFIX+filename), 'wb')
        w = csv.writer(writer)
        w.writerow(['Name', NAME, 'Target PDB', 'ResIndex', 'Similarity'] + key + ['Vector'])

        for k,v in self.heterodict.items():
            line = ['sth', 'NA', self.PDBname, 'NA', 'NA' ] + ['NA']*len(key) + [v['raw_vector']]
            w.writerow(line)
        writer.flush()
        writer.close()


    def pick_one(self, ResId, **kwargs):
        return self.heterodict[ResId] or None

    def list_ResId(self):
        return self.heterodict.keys()

    def __repr__(self):
        print self.PDBname+'({} hetero parts found)'.format(len(self.heterodict.keys()))

