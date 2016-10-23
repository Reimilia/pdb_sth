import prody as protein
import numpy as np
import os

'''
Background: It will be much easier to use 4*4 matrix instead of 3*3 matrix to do transformation.
This will help when the program needs to judge whether the ligand is in rotated boxes, (in case of loss of information)

To be simple , all matrix is right-associative i.e.
[x',y',z',1] = [x,y,z,1]*T

T is transformation matrix
'''
def rotation_matrix_by_x(theta):
    '''
    just return the rotation matrix with clockwise theta degree along x-axis
    :param coord:
    :param theta:
    :return:
    '''
    return np.matrix(
        [[1,0,0,0],
         [0,np.cos(theta),-np.sin(theta),0],
         [0,np.sin(theta),np.cos(theta),0],
         [0,0,0,1]])

def rotation_matrix_by_y(theta):
    '''
    just return the rotation matrix with clockwise theta degree along y-axis
    :param coord:
    :param theta:
    :return:
    '''
    return np.matrix(
        [[np.cos(theta),0,np.sin(theta),0],
         [0,1,0,0],
         [-np.sin(theta),0,np.cos(theta),0],
         [0,0,0,1]])

def rotation_matrix_by_z(theta):
    '''
    just return the rotation matrix with clockwise theta degree along z-axis
    :param coord:
    :param theta:
    :return:
    '''
    return np.matrix(
        [[np.cos(theta),-np.sin(theta),0,0],
         [np.sin(theta),np.cos(theta),0,0],
         [0,0,1,0]
         [0,0,0,1]])


def transition_matrix(coord):
    '''

    :param coord:
    :param shift:
    :return:
    '''
    return np.matrix(
        [[1,0,0,0],
         [0,1,0,0],
         [0,0,1,0],
         [coord[0],coord[1],coord[2],1]
         ]
    )

def get_rotation_marix_along_anyaxis(center,rotation):
    '''
    There is a vector as rotation axis, which starts at center and
    according to following steps to be that direction:
    0. It originally point along +z-axis
    1. It rotate clockwisely with rotation[0] degree along x-axis
    2. It rotate clockwisely with rotation[1] degree along y-axis
    3. It rotate clockwisely with rotation[2] degree along z-axis

    (P.S. 6-freedom)

    :param center: center vector (used in transition transformation)
    :param rotation: rotation vector (used in rotation transformation)
    :return:
    '''
    T = transition_matrix(-center)
    T_inv= transition_matrix(center)
    R_x = rotation_matrix_by_x(rotation[0])
    R_x_inv = rotation_matrix_by_x(-rotation[0])
    R_y = rotation_matrix_by_y(rotation[1])
    R_y_inv = rotation_matrix_by_y(-rotation[1])
    R_z = rotation_matrix_by_z(rotation[2])
    return reduce(np.dot, [T,R_x,R_y,R_z,R_y_inv,R_x_inv,T_inv])

def get_transformation_inv_fromoldtonew(new_xcoord,new_ycoord,new_zcoord,transition):
    '''

    :param new_xcoord: (1,0,0) -> (ux,uy,uz)
    :param new_ycoord: (0,1,0) -> (vx,vy,vz)
    :param new_zcoord: (0,0,1) -> (wx,wy,wz)
    :param transition: **very important** this is used to align the start (may be down-left coordinate after rotation)
    because we need vectors instead of absolute coordinate to calculate
    i.e. (a*ux+transition[0],a*uy+transition[1],a*uz+transition[2]) might represent one edge of a cubic box.
    :return: **Inverse matrix** of this transformation
    '''
    T_inv = transition_matrix(-transition)
    U= np.matrix([new_xcoord,new_ycoord,new_zcoord])
    return np.dot(T_inv*U.transpose())

def get_coord_after_transformation(coord,transformation_matrix):
    return np.dot(coord+[1],transformation_matrix)[0:2]



class Box:
    center= np.array([10,10,10])
    Boxsize= 20
    down_left = np.array([0,0,0])
    x_axis= np.array([1,0,0])
    y_axis= np.array([0,1,0])
    z_axis= np.array([0,0,1])

    # resolution of box
    Boxrange= 1
    Boxnum =  int(np.ceil(Boxsize / Boxrange))

    rotation = [0,0,0]
    rotated= False
    transition_matrix= None

    def __init__(self,**kwargs):
        if 'center' in kwargs:
            self.center= kwargs['center']
        if 'Boxsize' in kwargs:
            self.Boxsize = kwargs['Boxsize']
        if 'Boxrange' in kwargs:
            self.Boxrange = kwargs['Boxrange']
        self.down_left= self.center-[self.Boxsize/2,self.Boxsize/2,self.Boxsize/2]
        self.Boxnum=  int(np.ceil(self.Boxsize/ self.Boxrange))

    def change_size(self,Boxsize,Boxrange=1):
        self.Boxsize= Boxsize
        self.Boxrange= Boxrange
        self.Boxnum = int(np.ceil(Boxsize/Boxrange))

    def change_manually(self,x_vec,y_vec,z_vec,center):
        self.x_axis = x_vec
        self.y_axis = y_vec
        self.z_axis = z_vec
        self.center = center
        self.down_left = self.center- self.Boxsize/2 * (x_vec+y_vec+z_vec)

    def set_default(self,center,Boxsize,Boxrange):
        self.center = center
        self.Boxsize = Boxsize
        self.Boxrange = Boxrange

        self.down_left = self.center - [self.Boxsize / 2, self.Boxsize / 2, self.Boxsize / 2]
        self.Boxnum = int(np.ceil(self.Boxsize / self.Boxrange))
        self.x_axis = np.array([1, 0, 0])
        self.y_axis = np.array([0, 1, 0])
        self.z_axis = np.array([0, 0, 1])

        self.rotation = [0, 0, 0]
        self.rotated = False
        self.transition_matrix = np.ones(4)
        self.coordinate_shift_matrix = np.ones(4)

    def transform(self,rotation,**kwargs):
        '''

        :param rotation: degrees to rotate , always [x,y,z] x for 'along x-axis clockwisely', so on so forth
        :param kwargs: center : new center of box
                        transition : shift from old box center
                        note only one can be put as input
        :return:
        '''
        if 'center' in kwargs:
            if 'transition' in kwargs:
                print 'This will cause ambiguity. Give either new center or shift from the old center'
                return
            else:
                self.center= kwargs['center']
        else:
            if 'transition' in kwargs:
                self.center += kwargs['transition']
        T= get_rotation_marix_along_anyaxis(self.center,rotation)
        self.rotated = True
        if self.transition_matrix is None:
            self.transition_matrix = T
        else:
            self.transition_matrix *= T
        self.x_axis = get_coord_after_transformation(self.x_axis, T)
        self.y_axis = get_coord_after_transformation(self.y_axis, T)
        self.z_axis = get_coord_after_transformation(self.z_axis, T)

        self.coordinate_shift_matrix = get_transformation_inv_fromoldtonew(self.x_axis,self.y_axis,self.z_axis,self.down_left)

    def find_lattice_coord(self,coord):
        '''

        :param coord: [x,y,z] all be interger to represent grid (lattice) in box (down_lefy is [0,0,0])
        :return: real coords in xyz-coordinatesystem (Absolute one)
        '''
        if coord[0]<0 or coord[0]>=self.Boxsize:
            print 'x-coordinate out of border!'
            return
        if coord[1]<0 or coord[1]>=self.Boxsize:
            print 'y-coordinate out of border!'
            return
        if coord[2]<0 or coord[2]>=self.Boxsize:
            print 'z-coordinate out of border!'
            return

        print self.down_left + np.dot(coord,np.matrix([self.x_axis,self.y_axis,self.z_axis]))
        return self.down_left + np.dot(coord,np.matrix([self.x_axis,self.y_axis,self.z_axis]))

    def get_lattice_coord(self,coords):
        return get_coord_after_transformation(coords, self.transition_matrix)

    def in_cubic_box(self,coords):
        '''
        This is the function that used to judge whether the coords are in the cubix box or not
        :param coords:
        :param center:
        :param BOXsize:
        :return:
        '''
        relative_coords= get_coord_after_transformation(coords,self.coordinate_shift_matrix)
        print relative_coords
        return (0<=relative_coords[0]<self.Boxsize and 0<=relative_coords[1]<self.Boxsize and 0<=relative_coords[2]<self.Boxsize)



class vector_generator:
    receptor= None
    heterodict={}
    boxdict={}

    def __init__(self,receptor_filename):
        '''
        :param receptor_filename: source pdb file or .gz file that can be opened and parsed
        '''

        try:
            parse = protein.parsePDB(receptor_filename)
        except:
            print 'Cannot parse file %s, please check your input' % receptor_filename.split('/')[-1]
            return

        self.receptor= parse.select('protein')

        if parse.select('nucleic') is not None:
            print 'This program does not support nucleic cocrystal structure for now.'
            return

        hetero = parse.select('(hetero and not water) or resname ATP or resname ADP')

        for pick_one in protein.HierView(hetero).iterResidues():
            # less than 3 atoms may be not ok
            if pick_one.numAtoms() <= 3:
                continue

            ResId = pick_one.getResindex()
            self.heterodict[ResId]=pick_one

            # Set up a new box class here
            self.boxdict[ResId]= Box(center=protein.calcCenter(pick_one).getCoords(),Boxsize=20,Boxrange=1)

    def set_ligand_from_file(self,ligand_filename):
        '''
        :param ligand_filename: source pdb file or mol2 file that can be opened and parsed
         pdb might be better because mol2 will loss atom_name (most only indicate O/N/C)
        :return:
        '''
        try:
            suffix = ligand_filename.split('.')[-1]
        except:
            print 'ligand file parse error! check the address '+ ligand_filename
            return 'NA'

        try:
            if suffix!='pdb':
                #Note it will only try to parse first one
                os.system('babel -imol2 %s -opdb temp.pdb'% ligand_filename)
                ligand = protein.parsePDB('temp.pdb')
            else:
                ligand = protein.parsePDB(ligand_filename)
        except:
            print 'Only support pdb or mol2 format, please check your file. It will cast a parse Error!'
            return 'NA'

        ResId= ''.join(map(lambda xx:(hex(ord(xx))[2:]),os.urandom(16)))
        self.heterodict[ResId] = ligand
        return ResId

    def generate_vector_from_file(self,ligand_filename,try_threshold=200,shift_threshold=[1,1,1]):
        '''

        :param ligand_filename: where source file can be read and parsed , pdb or mol2 format only
        :param try_threshold: if random transition failed over this many times, the program will inform an error
        instead of return a valid vector. Note the first half trial times will only concern rotation. Then shift.
        :param shift_threshold: the absolute value of transition transform will not exceed this threshold
        :return: either a vector in random box 'on the fly' or get an Error (False, i.e.)
        '''
        ResId= self.set_ligand_from_file(ligand_filename)
        if ResId=='NA':
            return False, [0]*8000

        ligand = self.heterodict[ResId]

        Box = self.boxdict[ResId]
        ligand_center = protein.calcCenter(ligand).getCoords()

        for iteration_step in range(try_threshold/2):
            # Now try to do one rotation

            rotation=np.random.random_sample(3)*np.pi/4
            Box.transform(rotation)
            Tag= True
            for atom in ligand.iterAtoms():
                coord = atom.getCoords()
                if Box.in_cubic_box(coord)==False:
                    Tag= False
            if Tag == False:
                print 'Try %s : rotation failed to contain ligand.' % str(iteration_step)
                continue




        return False,[0]*8000