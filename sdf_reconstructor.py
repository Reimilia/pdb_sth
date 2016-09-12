__author__ = 'wy'

from rdkit import Chem
from prody import *
import os
import csv
from vector_gen import pdb_container
import time
from functools import wraps
import logging
import pandas as pd

'''
Generate the data.
'''

fileHandler = logging.FileHandler('debug.log')
fileHandler.setLevel(logging.WARNING)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)


def fn_timer(function):
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

MUST = "PDB ID(s) for Ligand-Target Complex"
NAME = 'BindingDB Reactant_set_id'

@fn_timer
def query_specific_field(src , key):
    '''
    This program will search sdf file in the specifc location or search in group
    and pick up the specific field u want from the extension part in a single entry
    of a sdf file or a group

    From now on it only supports basic one-field search
    :param key: the search field, need exact match
    :param src: two possible type :
                first : a filelocation must end in .sdf otherwise the result is nothing but '[]'
                second: a group of molecules (entries) , can be a list or bundle resource from rdkit library
                (i.e. the Chem.SDMolSupplier
    :return: the bundle result as a list, each is a mol in redit class-object format
             also, the result will be written in to a file in 'result' directory as a sdf file
    '''

    #Wipe the data space:
    files = os.listdir('data/')
    for filename in files:
        loc = os.path.join('data/') + filename
        if os.path.exists(loc):
            os.remove(loc)

    # write the result
    filedir = 'result/filter_{}'.format(src.split('/')[-1].split('.')[0])+'.csv'
    if not os.path.isfile(filedir):
        if not os.path.exists('result'):
            os.mkdir('result')

    # csv writer
    writer = file(filedir, 'wb')
    w = csv.writer(writer)
    w.writerow([NAME,'Target PDB','ResIndex','Similarity']+key+['Vector'])
    # first possible input : filename and its location
    if isinstance(src,str):
        try:
            sdf = Chem.SDMolSupplier(src)
        except:
            print 'filename is wrong!'
            raise IOError
    else:
        sdf= src

    assert isinstance(key, list)
    # This list will return all molecules that satisfy the needs
    count =0
    active_count=0
    pdb_num=0
    bad_one=0

    PDBindex={}


    try:
        for mol in sdf:
            if mol is None:
                continue
            count +=1
            if count == 100:
                break
            if mol.HasProp(MUST):
                value = mol.GetProp(MUST)
                if value =='':
                    continue
                PDBs=value.split(',')
                if mol.HasProp(NAME):
                    one_line = [mol.GetProp(NAME),'','','']
                else:
                    one_line = ['Unknown','','','']
                for subkey in key:
                    if mol.HasProp(subkey):
                        one_line.append(mol.GetProp(subkey))
                    else:
                        one_line.append('')
                ww =Chem.SDWriter('data/a.sdf')
                ww.write(mol)
                ww.close()
                one_line.append('')
                for pdb in PDBs:
                    ##make index to the PDBs
                    if pdb not in PDBindex:
                        logging.info('Writing index to pdb {}'.format(pdb))
                        PDBindex[pdb] = pdb_container(pdb)
                        pdb_num +=1

                    #print '\n\n{}\n\n'.format(pdb)
                    ans_list = PDBindex[pdb].find_similar_target('data/a.sdf')

                    for eachone in ans_list:
                        assert 'id' in eachone
                        assert 'cp' in eachone
                        assert 'raw_vector' in eachone
                        one_line[1] = pdb
                        one_line[2] = eachone['id']
                        one_line[3] = eachone['cp']
                        one_line[-1] = eachone['raw_vector']
                        # print one_line
                        active_count += 1
                        w.writerow(one_line)

                    if len(ans_list)==0:
                        bad_one+=1
    except:
        logging.error('Unknown error here!')
    logging.warning('{} bad ligands found'.format(bad_one))

    logging.warning('{} molecules are detected, {} are recorded.'.format(count,active_count))

    writer.close()
    return filedir

@fn_timer
def sort_the_result(result_filename):
    try:
        logging.info('Begin Sorting by pdb')
        df = pd.read_csv(result_filename, sep=',')
        df = df.sort_values(by=['Target PDB'])
        df.to_csv(result_filename, sep=',', float_format='%.4f')
        logging.info('End Sorting')
    except:
        logging.error('Cannot sort, maybe it is too large.')

if __name__ == '__main__':

    #Test

    BIG_DATA_SRC= '/media/wy/data/BindingDB_All_terse_3D.sdf'

    result=query_specific_field(BIG_DATA_SRC,
                         ['Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)','kon (M-1-s-1)','koff (s-1)'])
    #result = 'result/filter_BindingDBSampleBig.csv'

    sort_the_result(result)