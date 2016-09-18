__author__ = 'wy'

from prody import *
import os
import csv
from vector_gen import pdb_container
import time
from functools import wraps
import logging
from Config import *

'''
The main program to extract molecules in .sdf files and compare with ligands on PDB files.
Then accept all pairs with similarity >= 0.85 and generate the corresponding vectors.
'''

#This part is used to set debug log
fileHandler = logging.FileHandler('debug.log',mode='w')
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)




def fn_timer(function):
    '''
    This is the decorator used for time counting issue
    :param function:
    :return: no return. just print and record the time the devorated program ran.
    '''
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
def mol_ligand_tar_generator(src):
    '''
    This program will search sdf file in the specific location or search in group
    and pick up the specific field u want from the extension part in a single entry
    of a sdf file or a group

    From now on it only supports basic one-field search

    :param src: pdb name
    :return: no return, output will be in a .csv file with format : filter_[src].csv
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

    # in case for overwriting
    if os.path.exists(filedir):
        print '{} already done.'.format(src)
        logging.info('{} already done'.format(src))
        return True

    # csv writer
    writer = file(filedir, 'wb')
    w = csv.writer(writer)
    w.writerow(['Name',NAME,'Target PDB','ResIndex','Similarity']+key+['Vector'])

    # combine as file direction
    sdfone = filedir_PREFIX + src.upper() + '.sdf'

    try:
        input_sdf = open(sdfone,'r')
    except:
        logging.error('PDB {} with ligands sdf not found!'.format(src))
        return False

    assert isinstance(key, list)


    # This variables are used for counting and statistic issue.
    # TODO generate auto report using these numbers
    active_count=0
    count=0
    pdb_num=0
    bad_one=0

    # Combine as pdb file address
    # We generate a class to store each ligands as a dict, and use a method
    # to find the similar ones by tanimoto comparing scores to specific input files
    PDBindex= pdb_container(src,filepos=pdb_PREFIX+src.lower()+'.pdb.gz')

    # In each time
    # We write one single molecule in to a sdf file called a.sdf
    # Then use this file to compare with the one we extract from pdb files
    # Since we use scripts so I just extract them to make sure that is what
    # we really want to compare
    # (But there should be a smarter way to do so)
    o =open("a.sdf", "w")
    try:
        mol = ''
        LINE_BEGIN=True
        Wait_Signal= 0
        one_line=['']*Total_columns
        print 'here'
        for line in input_sdf:
            mol+=line

            #This part is finding columns need to be recorded in sdf files
            if LINE_BEGIN:
                one_line[0] = line.lstrip(' ').rstrip('\n')
                LINE_BEGIN = False

            # just lazy
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

                # Find pairs with at least 85% similarity scores
                ans_list =PDBindex.find_similar_target('a.sdf')

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


if __name__ == '__main__':

    #Test
    #ToDo Write a filemerger to merge all filter csv files into a desired bigone
    #ToDo Write a function or improve exising functions to combine the result from the same chemical compound.
    #ToDo Wrap as a function

    from Config import PDB_tar
    import urllib
    import time

    DONE=[]
    FAIL=[]
    ct=0


    for pdb in PDB_tar:

        #For each pdb name , there should be one corresponding molecule files.
        #This will generate one result file.
        pdb=pdb.lower()
        filename ='pdb_raw/{}.pdb.gz'.format(pdb)
        if os.path.exists(filename):
            #pdbfile exists
            print pdb+ ' has downloaded'
            tag=mol_ligand_tar_generator(pdb)
            if tag:
                DONE.append(pdb)
            else:
                FAIL.append(pdb)
            ct += 1
        else:
            #Not exists, download from the internet
            urllib.urlretrieve(url_prefix+'{}.pdb.gz'.format(pdb.lower()),filename)
            #Wait for 1 second from rejection on connection.
            time.sleep(1)

            # This is to check whether we download the file successfully
            o = open(filename, 'r')
            for l in o:
                if l.find('DOCTYPE') != -1:
                    print 'download {} failed'.format(pdb)
                    FAIL.append(pdb)
                    break
                else:
                    print 'download {} successfully'.format(pdb)
                    tag=mol_ligand_tar_generator(pdb)
                    if tag:
                        DONE.append(pdb)
                    else:
                        FAIL.append(pdb)
                    break
            o.close()

    print ct
    logging.info('total: {}'.format(ct))
    #Todo: Use DONE and FAIL wisely (generate a report)
