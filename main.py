from Config import *
import os
import csv
import logging
from vector_gen import pdb_container,fake_pdb_container
import time
from functools import wraps

#This part is used to set debug log
#This will generate a log that record every content logged with a specific security levels
fileHandler = logging.FileHandler('debug.log',mode='w')
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)


def initiate_report():
    csv_name = 'report.csv'
    writer = file(csv_name, 'wb')
    w = csv.writer(writer)
    w.writerow(['filename','pdb Name','molecules','paired','bad_one','pairtimes'])
    return csv_name
def fn_timer(function):
    '''
    This is the decorator used for time counting issue
    Need not understand this one. It has nothing to do with generating files
    :param function:
    :return: no return. just print and record the time the decorated program ran.
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


def bundle_result_mol2_file(source_mol_file ,experimentaldict, pdbdict):
    '''
    Generate the mol2 file with additional infomation in result folders
    :param ResId:
    :param remarkdict:
    :return:
    '''
    assert os.path.exists(source_mol_file)
    filename = source_mol_file.split('/')[-1]

    real_dir= os.path.join(result_PREFIX,filename)

    # first deal with experimental data
    comment = 'Remark: '
    if experimentaldict is not None:
        for k,v in experimentaldict.items():
            comment = comment + '|>'+ k + ':'+v
    if pdbdict is not None:
        for k,v in experimentaldict.items():
            comment = comment + '|>'+ k + ':'+v

    with open(real_dir,'wb') as w:
        w.write('# '+comment+'\n')
        with open(source_mol_file,'rb') as o:
            w.writelines(o.read())

    print '{} is moved into the result file with addtional info in first line.'.format(filename)






@fn_timer
def bindingDB_pdb_tar_generator(src,filepos,statistic_csv=None,CLEAN=False,fileforbabel='a.sdf'):
    '''

    :param src: pdb name
    :param statistic_csv: the report csv file's name
    :param CLEAN: Wipe temporary pdb files or not. Note I will not give options to wipe results. That's dangerous
    :return: True: If everything works fine
             False: Unexpected error happens. Note if there is no reuslt, it will return True because everything runs fine.
    '''

    # Wipe the pdb temporary files if you wish:
    if CLEAN:
        files = os.listdir('data/')
        for filename in files:
            loc = os.path.join('data/') + filename
            if os.path.exists(loc):
                os.remove(loc)

    # write the result

    result_file_name = 'filter_{}'.format(src.split('/')[-1].split('.')[0]) + '.csv'
    filedir = os.path.join(result_PREFIX, result_file_name)
    if not os.path.isfile(filedir):
        if not os.path.exists(result_PREFIX):
            os.mkdir(result_PREFIX)

    # in case for overwriting
    '''
    if os.path.exists(filedir):
        print '{} already done.'.format(src)
        logging.info('{} already done'.format(src))
        return True
    '''

    # combine as file direction
    sdfone = filedir_PREFIX + src.upper() + '.sdf'

    # open the source molecule files
    # Naming format [PDB name].sdf all lowercase
    try:
        input_sdf = open(sdfone, 'r')
    except:
        logging.error('PDB {} with ligands sdf not found!'.format(src))
        return False

    # This variables are used for counting and statistic issue.
    active_count = 0
    count = 0
    bad_one = 0

    # Combine as pdb file address
    # We generate a class to store each ligands as a dict, and use a method
    # to find the similar ones by tanimoto comparing scores to specific input files
    PDBindex = pdb_container(src, filepos=filepos, BOX=21, Size=0.35)

    if PDBindex.get_pdb_type() != 'Protein':
        return False
    # In each time
    # We write one single molecule in to a sdf file called a.sdf
    # Then use this file to compare with the one we extract from pdb files
    # Since we use scripts so I just extract them to make sure that is what
    # we really want to compare
    # (But there should be a smarter way to do so)

    try:
        mol = ''
        Wait_Signal = 0
        experimental_data = {}
        experiment_dict= {} # print 'here'
        for line in input_sdf:
            mol += line


            # just lazy
            if Wait_Signal > 0:
                experiment_dict[line_key]=line.lstrip(' ').rstrip(' ').rstrip('\n')
                Wait_Signal = 0

            for i in range(len(key)):
                if key[i] in line:
                    Wait_Signal = 1
                    line_key = key[i]
                    break


            if '$$$$' in line:
                # end of a molecule
                assert monomerID is not None
                if monomerID not in experimental_data:
                    fileforbabel = temp_pdb_PREFIX + '/{}/{}_{}.sdf'.format(src, src, monomerID)
                    o = open(fileforbabel, "a")
                    o.write(mol)
                    o.close()
                    experimental_data[monomerID] = one_line
                    one_line = [''] * len(experiment_part)
                else:
                    # combine experimental data together
                    for i in range(len(key)):
                        if experimental_data[monomerID][3 + i] == '':
                            experimental_data[monomerID][3 + i] = one_line[3 + i]
                        else:
                            if len(one_line[3 + i]) > 0:
                                experimental_data[monomerID][3 + i] += '|' + one_line[3 + i]

                mol = ''
                LINE_BEGIN = False
                monomerID = None

        for k, v in experimental_data.items():
            # print k,v
            fileforbabel = temp_pdb_PREFIX + '/{0}/{0}_{1}.sdf'.format(src, k)
            # Find pairs with at least 85% similarity scores
            ans_list = PDBindex.find_similar_target(fileforbabel)
            # print 'here'
            count += 1
            for eachone in ans_list:
                # Combine each part together
                v[2] = eachone['cp']
                active_count += 1
                w.writerow(v + PDBindex.bundle_result(eachone['id']))

            if len(ans_list) == 0:
                bad_one += 1
                logging.info('not found ligand here: {}_{}.'.format(src, one_line[1]))

    except:
        # raise TypeError
        logging.error('Unknown error here!')
        return False
    logging.warning('{} bad ligands found'.format(bad_one))
    logging.warning('{} molecules are detected, and {} pairs are recorded.'.format(count, active_count))

    # Discard unused one
    #PDBindex.clean_temp_data()

    # Do some record
    if statistic_csv is not None:
        writer = file(statistic_csv, 'ab')
        w = csv.writer(writer)
        w.writerow([src, count, count - bad_one, bad_one, active_count])
        writer.flush()
        writer.close()
    return True