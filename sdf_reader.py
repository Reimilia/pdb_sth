__author__ = 'wy'

from rdkit import Chem
import json
import os, io
import time
from functools import wraps

'''
Basic example of how to use rdkit library
need to first install it follow the instruction on its website

'''





def query_specific_field(src , key, value=None):
    '''
    This program will search sdf file in the specifc location or search in group
    and pick up the specific field u want from the extension part in a single entry
    of a sdf file or a group

    From now on it only supports basic one-field search
    :param key: the search field, need exact match
    :param value: the target value, can be a part in the string (e.g. one pdb id in a string)
    :param src: two possible type :
                first : a filelocation must end in .sdf otherwise the result is nothing but '[]'
                second: a group of molecules (entries) , can be a list or bundle resource from rdkit library
                (i.e. the Chem.SDMolSupplier
    :return: the bundle result as a list, each is a mol in redit class-object format
             also, the result will be written in to a file in 'result' directory as a sdf file
    '''

    # write the result
    dir = 'result/filter_{}_{}'.format(key.split(' ')[0] + key.split(' ')[-1], value) + '.sdf'
    if not os.path.isfile(dir):
        if not os.path.exists('result'):
            os.mkdir('result')

    # sdf writer
    w = Chem.SDWriter(dir)

    # first possible input : filename and its location
    if isinstance(src,str):
        try:
            sdf = Chem.SDMolSupplier(src)
        except:
            print 'filename is wrong!'
            raise IOError
    else:
        sdf= src

    assert isinstance(key, str)
    # This list will return all molecules that satisfy the needs
    ans_list = []
    try:
        for mol in sdf:
            if mol is None:
                print 'Incomplete data?'
                continue
            if mol.HasProp(key):
                src_value = mol.GetProp(key)

                # write and store the corresponding files
                if value is None:
                    ans_list.append(mol)
                    w.write(mol)
                else:
                    if value == src_value or (
                                    isinstance(value, str) and isinstance(src_value, str) and value in src_value):
                        ans_list.append(mol)
                        w.write(mol)
    except:
        raise TypeError
    return ans_list


if __name__ == '__main__':

    #Both are correct usage

    sdf = Chem.SDMolSupplier('BindingDBSampleBig.sdf')
    filtered_result = query_specific_field(sdf, 'PDB ID(s) for Ligand-Target Complex', '2H6T')
    '''
    #Test



    sdf = Chem.SDMolSupplier('/media/wy/data/BindingDB_All_terse_3D.sdf')
    count =0
    for mol in sdf:
        if mol is not None:
            count+=1
    print count


    for mol in filtered_result:
        # This will show how to get protery names , basic info in a molecule class
        # For more info, see the website of rdkit

        # dict for property
        dict = {}
        for name in mol.GetPropNames():
            dict[name] = mol.GetProp(name)
        print json.dumps(dict, indent=4)

        # SMILE format
        print Chem.MolToSmiles(mol)

        # numofAtoms
        print mol.GetNumAtoms()

        # Atom and bonds
        print mol.GetAtomWithIdx(0).GetSymbol()
        print mol.GetAtomWithIdx(0).GetExplicitValence()
        print mol.GetBondWithIdx(0).GetBeginAtomIdx()
        print mol.GetBondWithIdx(0).GetEndAtomIdx()
        print mol.GetBondBetweenAtoms(0, 1).GetBondType()

        # This is useless for getting info
    '''