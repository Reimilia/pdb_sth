import prody as protein
import numpy as np

def rotation_by_x(coord, theta):
    '''
    Make the coordinate shift with rotation to x-axis
    :param coord:
    :param theta:
    :return:
    '''
    T = np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])
    return np.dot(coord,T)


def rotation_by_y(coord, theta):
    '''
    Make the coordinate shift with rotation to y-axis
    :param coord:
    :param theta:
    :return:
    '''
    T = np.array([[np.cos(theta), 0, -np.sin(theta)], [0, 1, 0], [np.sin(theta),0,np.cos(theta)]])
    return np.dot(coord, T)


def rotation_by_z(coord, theta):
    '''
    Make the coordinate shift with rotation to z-axis
    :param coord:
    :param theta:
    :return:
    '''
    T = np.array([[np.cos(theta), np.sin(theta),0], [-np.sin(theta),np.cos(theta),0],[0, 0, 1]])
    return np.dot(coord, T)


def transition_by_vector(coord, shift):
    '''

    :param coord:
    :param shift:
    :return:
    '''
    return coord+shift

def get_relative_coordinate_in_transform(coord, rotation, transition):
    '''

    :param coord:
    :param rotation:
    :param transition:
    :return:
    '''
    tar= transition_by_vector(coord, transition)
    tar= rotation_by_x(tar  , -rotation[0])
    tar = rotation_by_x(tar , -rotation[1])
    tar = rotation_by_x(tar , -rotation[2])
    return tar

def in_cubic_box(coords, center, BOXsize= 20):
    pass

class vector_generator:
    receptor= None

    def __init__(self,receptor_filename):
        '''
        :param receptor_filename: source pdb file or .gz file that can be opened and parsed
        '''

        try:
            parse = protein.parsePDB(receptor_filename)
        except:
            print 'Cannot parse file %s, please check your input' % receptor_filename.split('/')[-1]
            return False

        self.receptor= parse.select('protein')

        if parse.select('nucleic') is not None:
            print 'This program does not support nucleic cocrystal structure for now.'
            return False



