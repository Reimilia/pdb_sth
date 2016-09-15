from prody import *
import random
import numpy as np
import os,re
import logging

from rdkit import Chem

HETERO_PART =1
PROTEIN_PART = 2

CONFIDENCE = 0.85

class vector_generator:
    def __init__(self, PDBname, **kwargs):
        parse = parsePDB(PDBname)
        hetero = parse.select('hetero and not water')
        self.name = PDBname
        self.vectors = []
        self.raw_vectors = []
        self.middles = []
        self.dictionary = []
        self.resname = []

        for pick_one in HierView(hetero).iterResidues():
            if pick_one.numAtoms() <= 3:
                continue


            self.raw_vectors.append(pick_one.getCoords())

            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)

            scale = max(max(xyz[:, 0]) - middle[0], middle[0] - min(xyz[:, 0]),
                        max(xyz[:, 1]) - middle[1], middle[1] - min(xyz[:, 1]),
                        max(xyz[:, 2]) - middle[2], middle[2] - min(xyz[:, 2]))
            assert scale <= 10

            xx, yy, zz = np.meshgrid(np.linspace(middle[0] - 9.5, middle[0] + 9.5, 20),
                                 np.linspace(middle[1] - 9.5, middle[1] + 9.5, 20),
                                 np.linspace(middle[2] - 9.5, middle[2] + 9.5, 20))


            #print xx
            vector = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

            self.middles.append(np.array([vector[0][0],vector[0][1],vector[0][2]]))
            #print vector[0][0],vector[0][1],vector[0][2]

            dic = {}
            pos_list= []

            num_vector = [0]*8000
            for atom in pick_one.iterAtoms():
                x, y, z = atom.getCoords()
                x_pos = int(round(x - vector[0][0]))
                assert 0<=x_pos <= 19
                y_pos = int(round(y - vector[0][1]))
                assert 0<=y_pos <= 19
                z_pos = int(round(z - vector[0][2]))
                assert 0<=z_pos <= 19
                num_vector[x_pos * 400 + y_pos * 20 + z_pos] = atom.getIndex()
                dic[str(atom.getIndex())] = atom.getCoords()
                pos_list.append([vector[0][0]+x_pos,vector[0][1]+y_pos,vector[0][2]+z_pos])
                print pos_list[-1]

            #Find the atom in protein in cubic boxes
            self.dictionary.append(dic)
            self.vectors.append(num_vector)
            self.resname.append(pick_one.getResindex())

    def get_one(self,index=None):
        random.seed(0)
        length = len(self.vectors)
        if not index:
            index = random.randint(0, length - 1)
        return self.vectors[index], self.filenames[index]


class pdb_container:
    def __init__(self,PDB,filepos=None,**kwargs):
        self.PDBname= PDB
        self.heterodict = {}

        try:
            if filepos is not None:
                parse = parsePDB(filepos)
            else:
                parse = parsePDB(PDB)
        except:
            #raise IOError
            logging.warning('PDB {} is ignored due to file-not-found error'.format(PDB))
            return
        hetero = parse.select('not protein and not nucleic and not water')
        other = parse.select('protein or nucleic')


        if hetero is None or other is None:
            return

        #Make index for every single hetero parts
        for pick_one in HierView(hetero).iterResidues():
            if pick_one.numAtoms() <= 3:
                continue

            ResId = str(pick_one.getResindex())
            filename = 'data/{}_{}.pdb'.format(PDB, ResId)
            writePDB(filename, pick_one)



            xyz = pick_one.getCoords()
            middle = calcCenter(pick_one)


            scale = max(max(xyz[:, 0]) - middle[0], middle[0] - min(xyz[:, 0]),
                        max(xyz[:, 1]) - middle[1], middle[1] - min(xyz[:, 1]),
                        max(xyz[:, 2]) - middle[2], middle[2] - min(xyz[:, 2]))
            #assert scale <= 10
            if scale>10:
                logging.warning('Warning! {} has a ligand out of box scale with {} atom distance to center'.format(PDB, scale))
                #ignore them for now
                continue

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
                    num_vector[x_pos * 400 + y_pos * 20 + z_pos] = HETERO_PART


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
                        num_vector[x_pos * 400 + y_pos * 20 + z_pos] = PROTEIN_PART
                    else:
                        logging.warning('Coorinate {} {} {} found'.format(x_pos,y_pos,z_pos))

                #filename2= 'data/{}_{}_2.pdb'.format(PDB, ResId)

                #writePDB(filename2, pick_one+nearby)

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
        :param sdf_filedir:l
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
                raise TypeError
                #sdf = Chem.SDMolSupplier(sdf_filedir)
                #w = Chem.SDWriter('data/{}.sdf'.format(v['filename']), 'w')
                #for mol in sdf:
                #    w.write(mol)
                #w.close()
                logging.warning('Babel encountered a problem at pdb {} ligand {}'.format(self.PDBname,v['filename']))
                cp = 0
            if float(cp) >= 0.85:
                v['cp']=float(cp)
                possible_ones.append(v)

        return possible_ones

    def pick_one(self, ResId, **kwargs):
        return self.heterodict[ResId] or None

    def list_ResId(self):
        return self.heterodict.keys()