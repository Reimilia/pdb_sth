__author__ = 'wy'

from rdkit import Chem
from prody import *
import os
import csv
from vector_gen import pdb_container
import time
from functools import wraps
import logging

'''
Generate the data.
'''

fileHandler = logging.FileHandler('debug.log',mode='w')
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)



url_prefix = 'https://files.rcsb.org/download/'
filedir_PREFIX=  '/media/wy/data/raw_data/'
pdb_PREFIX = 'pdb_raw/'

MUST = "PDB ID(s) for Ligand-Target Complex"
NAME = 'BindingDB Reactant_set_id'
Total_columns= 12
key= ['Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)','kon (M-1-s-1)','koff (s-1)']


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" %
               (function.func_name, str(t1 - t0))
               )
        logging.warning ("Total time running %s: %s seconds" %
               (function.func_name, str(t1 - t0))
               )
        return result

    return function_timer

@fn_timer
def query_specific_field(src):
    '''
    This program will search sdf file in the specifc location or search in group
    and pick up the specific field u want from the extension part in a single entry
    of a sdf file or a group

    From now on it only supports basic one-field search

    :param src: pdb name
    :return: the bundle result as a list, each is a mol in redit class-object format
             also, the result will be written in to a file in 'result' directory as a sdf file
    '''

    #Wipe the data space:
    '''files = os.listdir('data/')
    for filename in files:
        loc = os.path.join('data/') + filename
        if os.path.exists(loc):
            os.remove(loc)'''

    # write the result
    filedir = 'result/filter_{}'.format(src.split('/')[-1].split('.')[0])+'.csv'
    if not os.path.isfile(filedir):
        if not os.path.exists('result'):
            os.mkdir('result')
    if os.path.exists(filedir):
        print '{} already done.'.format(src)
        logging.info('{} already done'.format(src))
        return True

    # csv writer
    writer = file(filedir, 'wb')
    w = csv.writer(writer)
    w.writerow(['Name',NAME,'Target PDB','ResIndex','Similarity']+key+['Vector'])
    # first possible input : filename and its location
    print 'here'
    sdfone = filedir_PREFIX + src.upper() + '.sdf'


    try:
        input_sdf = open(sdfone,'r')
    except:
        logging.error('PDB {} with ligands sdf not found!'.format(src))
        return False

    assert isinstance(key, list)
    # This list will return all molecules that satisfy the needs

    active_count=0
    count=0
    pdb_num=0
    bad_one=0

    PDBindex= pdb_container(src,filepos=pdb_PREFIX+src.lower()+'.pdb.gz')

    o =open( "a.sdf", "w")
    try:
        mol = ''
        LINE_BEGIN=True
        Wait_Signal= 0
        one_line=['']*Total_columns
        print 'here'
        for line in input_sdf:
            mol+=line
            if LINE_BEGIN:
                one_line[0] = line.lstrip(' ').rstrip('\n')
                LINE_BEGIN = False
            if Wait_Signal>0:
                if Wait_Signal==999:
                    one_line[1]=line.lstrip(' ').rstrip('\n')
                else:
                    one_line[4+Wait_Signal]= line.lstrip(' ').rstrip('\n')
                Wait_Signal= 0

            for i in range(len(key)):
                if key[i] in line:
                    Wait_Signal=i+1
                    break
            if NAME in line:
                Wait_Signal=999

            if '$$$$' in line:
                #end of a molecule
                o.write(mol)
                o.close()
                print 'here'
                ans_list =PDBindex.find_similar_target('a.sdf')
                print 'there'
                count +=1
                for eachone in ans_list:
                    assert 'id' in eachone
                    assert 'cp' in eachone
                    assert 'raw_vector' in eachone
                    one_line[2] = src
                    one_line[3] = eachone['id']
                    one_line[4] = eachone['cp']
                    one_line[-1] = eachone['raw_vector']
                    # print one_line
                    active_count += 1
                    w.writerow(one_line)

                if len(ans_list) == 0:
                    bad_one += 1
                    logging.info('not found ligand here: {}_{}.'.format(src,one_line[1]))

                mol = ''
                LINE_BEGIN=False
                o = open("a.sdf", "w")

    except:
        logging.error('Unknown error here!')
        return False
    logging.warning('{} bad ligands found'.format(bad_one))
    logging.warning('{} proteins are used'.format(pdb_num))
    logging.warning('{} molecules are detected, and {} pairs are recorded.'.format(count,active_count))

    writer.flush()
    writer.close()
    return True

def do_one(pdb):
    pdb = pdb.lower()
    filename = 'pdb_raw/{}.pdb.gz'.format(pdb)
    if os.path.exists(filename):
        print pdb + ' has downloaded'
        tag = query_specific_field(pdb)
        return tag
    else:
        urllib.urlretrieve(url_prefix + '{}.pdb.gz'.format(pdb.lower()), filename)
        time.sleep(1)
        o = open(filename, 'r')
        for l in o:
            if l.find('DOCTYPE') != -1:
                print 'download {} failed'.format(pdb)
                return False
                break
            else:
                print 'download {} successfully'.format(pdb)
                tag = query_specific_field(pdb)
                return tag
                break
        o.close()


if __name__ == '__main__':

    #Test

    from Config import PDB_tar
    import urllib
    import time

    DONE=[]
    FAIL=[]
    ct=0


    for pdb in PDB_tar:
        pdb=pdb.lower()
        filename ='pdb_raw/{}.pdb.gz'.format(pdb)
        if os.path.exists(filename):
            print pdb+ ' has downloaded'
            tag=query_specific_field(pdb)
            if tag:
                DONE.append(pdb)
            ct += 1
        else:
            urllib.urlretrieve(url_prefix+'{}.pdb.gz'.format(pdb.lower()),filename)
            time.sleep(1)
            o = open(filename, 'r')
            for l in o:
                if l.find('DOCTYPE') != -1:
                    print 'download {} failed'.format(pdb)
                    FAIL.append(pdb)
                    break
                else:
                    print 'download {} successfully'.format(pdb)
                    tag=query_specific_field(pdb)
                    if tag:
                        DONE.append(pdb)
                    break
            o.close()

    print ct
    logging.info('total: {}'.format(ct))__author__ = 'wy'

from rdkit import Chem
from prody import *
import os
import csv
from vector_gen import pdb_container
import time
from functools import wraps
import logging

'''
Generate the data.
'''

fileHandler = logging.FileHandler('debug.log',mode='w')
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)



url_prefix = 'https://files.rcsb.org/download/'
filedir_PREFIX=  '/media/wy/data/raw_data/'
pdb_PREFIX = 'pdb_raw/'

MUST = "PDB ID(s) for Ligand-Target Complex"
NAME = 'BindingDB Reactant_set_id'
Total_columns= 12
key= ['Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)','kon (M-1-s-1)','koff (s-1)']


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" %
               (function.func_name, str(t1 - t0))
               )
        logging.warning ("Total time running %s: %s seconds" %
               (function.func_name, str(t1 - t0))
               )
        return result

    return function_timer

@fn_timer
def query_specific_field(src):
    '''
    This program will search sdf file in the specifc location or search in group
    and pick up the specific field u want from the extension part in a single entry
    of a sdf file or a group

    From now on it only supports basic one-field search

    :param src: pdb name
    :return: the bundle result as a list, each is a mol in redit class-object format
             also, the result will be written in to a file in 'result' directory as a sdf file
    '''

    #Wipe the data space:
    '''files = os.listdir('data/')
    for filename in files:
        loc = os.path.join('data/') + filename
        if os.path.exists(loc):
            os.remove(loc)'''

    # write the result
    filedir = 'result/filter_{}'.format(src.split('/')[-1].split('.')[0])+'.csv'
    if not os.path.isfile(filedir):
        if not os.path.exists('result'):
            os.mkdir('result')
    if os.path.exists(filedir):
        print '{} already done.'.format(src)
        logging.info('{} already done'.format(src))
        return True

    # csv writer
    writer = file(filedir, 'wb')
    w = csv.writer(writer)
    w.writerow(['Name',NAME,'Target PDB','ResIndex','Similarity']+key+['Vector'])
    # first possible input : filename and its location
    print 'here'
    sdfone = filedir_PREFIX + src.upper() + '.sdf'


    try:
        input_sdf = open(sdfone,'r')
    except:
        logging.error('PDB {} with ligands sdf not found!'.format(src))
        return False

    assert isinstance(key, list)
    # This list will return all molecules that satisfy the needs

    active_count=0
    count=0
    pdb_num=0
    bad_one=0

    PDBindex= pdb_container(src,filepos=pdb_PREFIX+src.lower()+'.pdb.gz')

    o =open( "a.sdf", "w")
    try:
        mol = ''
        LINE_BEGIN=True
        Wait_Signal= 0
        one_line=['']*Total_columns
        print 'here'
        for line in input_sdf:
            mol+=line
            if LINE_BEGIN:
                one_line[0] = line.lstrip(' ').rstrip('\n')
                LINE_BEGIN = False
            if Wait_Signal>0:
                if Wait_Signal==999:
                    one_line[1]=line.lstrip(' ').rstrip('\n')
                else:
                    one_line[4+Wait_Signal]= line.lstrip(' ').rstrip('\n')
                Wait_Signal= 0

            for i in range(len(key)):
                if key[i] in line:
                    Wait_Signal=i+1
                    break
            if NAME in line:
                Wait_Signal=999

            if '$$$$' in line:
                #end of a molecule
                o.write(mol)
                o.close()
                print 'here'
                ans_list =PDBindex.find_similar_target('a.sdf')
                print 'there'
                count +=1
                for eachone in ans_list:
                    assert 'id' in eachone
                    assert 'cp' in eachone
                    assert 'raw_vector' in eachone
                    one_line[2] = src
                    one_line[3] = eachone['id']
                    one_line[4] = eachone['cp']
                    one_line[-1] = eachone['raw_vector']
                    # print one_line
                    active_count += 1
                    w.writerow(one_line)

                if len(ans_list) == 0:
                    bad_one += 1
                    logging.info('not found ligand here: {}_{}.'.format(src,one_line[1]))

                mol = ''
                LINE_BEGIN=False
                o = open("a.sdf", "w")

    except:
        logging.error('Unknown error here!')
        return False
    logging.warning('{} bad ligands found'.format(bad_one))
    logging.warning('{} proteins are used'.format(pdb_num))
    logging.warning('{} molecules are detected, and {} pairs are recorded.'.format(count,active_count))

    writer.flush()
    writer.close()
    return True

def do_one(pdb):
    pdb = pdb.lower()
    filename = 'pdb_raw/{}.pdb.gz'.format(pdb)
    if os.path.exists(filename):
        print pdb + ' has downloaded'
        tag = query_specific_field(pdb)
        return tag
    else:
        urllib.urlretrieve(url_prefix + '{}.pdb.gz'.format(pdb.lower()), filename)
        time.sleep(1)
        o = open(filename, 'r')
        for l in o:
            if l.find('DOCTYPE') != -1:
                print 'download {} failed'.format(pdb)
                return False
                break
            else:
                print 'download {} successfully'.format(pdb)
                tag = query_specific_field(pdb)
                return tag
                break
        o.close()


if __name__ == '__main__':

    #Test

    from Config import PDB_tar
    import urllib
    import time

    DONE=[]
    FAIL=[]
    ct=0


    for pdb in PDB_tar:
        pdb=pdb.lower()
        filename ='pdb_raw/{}.pdb.gz'.format(pdb)
        if os.path.exists(filename):
            print pdb+ ' has downloaded'
            tag=query_specific_field(pdb)
            if tag:
                DONE.append(pdb)
            ct += 1
        else:
            urllib.urlretrieve(url_prefix+'{}.pdb.gz'.format(pdb.lower()),filename)
            time.sleep(1)
            o = open(filename, 'r')
            for l in o:
                if l.find('DOCTYPE') != -1:
                    print 'download {} failed'.format(pdb)
                    FAIL.append(pdb)
                    break
                else:
                    print 'download {} successfully'.format(pdb)
                    tag=query_specific_field(pdb)
                    if tag:
                        DONE.append(pdb)
                    break
            o.close()

    print ct
    logging.info('total: {}'.format(ct))