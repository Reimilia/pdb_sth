from prody import *
import random
import numpy as np


class vector_generator:
    def __init__(self, PDBname, isfile = False):
        parse = parsePDB(PDBname)
        hetero = parse.select('hetero')
        self.name = PDBname
        self.vectors = []
        self.raw_vectors = []
        self.middles = []
        self.dictionary = []
        self.filenames = []

        for pick_one in HierView(hetero).iterResidues():
            if pick_one.numAtoms() <= 1:
                continue

            self.raw_vectors.append(pick_one.getCoords())
            subfilename= '{}_hetero_{}_{}.pdb'.format(PDBname,pick_one.getChid(),pick_one.getResindex())
            self.filenames.append(subfilename)

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
            newhetero = AtomGroup('Hetero')
            newhetero.setCoords(np.array(pos_list,dtype=float))
            #print newhetero.__repr__()
            newhetero.setNames(['H']*len(pos_list))
            newhetero.setResnums([1]*len(pos_list))
            newhetero.setResnames(['LIG']*len(pos_list))

            writePDB(subfilename,newhetero)

            self.dictionary.append(dic)
            self.vectors.append(num_vector)

    def get_one(self,index=None):
        random.seed(0)
        length = len(self.vectors)
        if not index:
            index = random.randint(0, length - 1)
        return self.vectors[index], self.filenames[index]
