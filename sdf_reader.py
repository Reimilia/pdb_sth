__author__ = 'wy'

from rdkit import Chem
import json
import os,io

'''
Basic example of how to use rdkit library
need to first install it follow the instruction on its website

'''

'''
suppl = Chem.SDMolSupplier('BindingDBSampleSmall.sdf')

for mol in suppl:
    dict={}
    for name in mol.GetPropNames():
        dict[name]= mol.GetProp(name)
    print json.dumps(dict,indent=4)
    print Chem.MolToSmiles(mol)
    print '\n\n\n'
'''

def query_by_group(suppl, key, value=None):
    '''
    This program will search list suppl get all required part into one single bundle.
    :param suppl: a Chem.SDMolSupplier or a list of Mol
    :param key: the search field, need exact match
    :param value: the target value, can be a part in the string (e.g. one pdb id in a string)
    :return: the bundle and write into a single file
    '''
    ans_list=[]
    dir = 'result/filter_{}_{}'.format(key.split(' ')[0] + key.split(' ')[-1], value) + '.sdf'
    if not os.path.isfile(dir):
        if not os.path.exists('result'):
            os.mkdir('result')
    w = Chem.SDWriter(dir)
    try:
        for mol in suppl:
            assert isinstance(key, str)
            assert mol is not None
            if mol.HasProp(key):
                src_value = mol.GetProp(key)

            # write and store the corresponding files
            if value is None:
                ans_list.append(mol)
                w.write(mol)
            else:

                if value == src_value or (isinstance(value, str) and isinstance(src_value, str) and value in src_value):
                    ans_list.append(mol)
                    w.write(mol)
    except:
        raise TypeError

    return ans_list

def query_by_file(filedir, filename, key, value=None):
    '''
    This program will search sdf file in the specifc location and
    get all required part into one single bundle.

    From now on it only supports basic one field search
    :param key: the search field, need exact match
    :param value: the target value, can be a part in the string (e.g. one pdb id in a string)
    :param filename: the name of file , '*' means search all
    :param filedir: the location of a specific file
    :return: the bundle and write into a single file
    '''

    dir = 'result/filter_{}_{}'.format(key.split(' ')[0]+key.split(' ')[-1],value)+'.sdf'
    if not os.path.isfile(dir):
        if not os.path.exists('result'):
            os.mkdir('result')
    w = Chem.SDWriter(dir)
    if filename=='*':
        listedfile = os.listdir(filedir)
    else:
        listedfile = [os.path.join(filedir,filename)]
    #This list will return all molecules that satisfy the needs
    ans_list = []
    for filename in listedfile:
        if filename.split('.')[-1] == 'sdf':
            try:
                sdf = Chem.SDMolSupplier(os.path.join(filedir,filename))
            except:
                raise IOError
            for mol in sdf:
                assert isinstance(key,str)
                assert mol is not None
                if mol.HasProp(key):
                    src_value= mol.GetProp(key)

                #write and store the corresponding files
                if value is None:
                    ans_list.append(mol)
                    w.write(mol)
                else:

                    if value==src_value or (isinstance(value,str) and isinstance(src_value,str) and value in src_value):
                        ans_list.append(mol)
                        w.write(mol)
    return ans_list

if __name__=='__main__':



    filtered_result= query_by_file('.','BindingDBSampleBig.sdf','PDB ID(s) for Ligand-Target Complex','2H6T')
    for mol in filtered_result:
        #This will show how to get protery names , basic info in a molecule class
        #For more info, see the website of rdkit

        #dict for property
        dict={}
        for name in mol.GetPropNames():
            dict[name]= mol.GetProp(name)
        print json.dumps(dict,indent=4)

        #SMILE format
        print Chem.MolToSmiles(mol)

        #numofAtoms
        print mol.GetNumAtoms()

        #Atom and bonds
        print mol.GetAtomWithIdx(0).GetSymbol()
        print mol.GetAtomWithIdx(0).GetExplicitValence()
        print mol.GetBondWithIdx(0).GetBeginAtomIdx()
        print mol.GetBondWithIdx(0).GetEndAtomIdx()
        print mol.GetBondBetweenAtoms(0, 1).GetBondType()

        #This is useless for getting info
        print mol