__author__ = 'wy'

from prody import *
import os
import csv
from vector_gen import pdb_container,fake_pdb_container
import time
from functools import wraps
import logging
from Config import *
import urllib
from mapping import repair_pdbfile

'''
The main program to extract molecules in .sdf files and compare with ligands on PDB files.
Then accept all pairs with similarity >= 0.85 and generate the corresponding vectors.
'''

#This part is used to set debug log
#This will generate a log that record every content logged with a specific security levels
fileHandler = logging.FileHandler('debug.log',mode='w')
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('LINE %(lineno)-4d  %(levelname)-8s %(message)s', '%m-%d %H:%M')
fileHandler.setFormatter(formatter)
logging.getLogger('').addHandler(fileHandler)

def list_formatter(table):
    '''
    I don't know if there is a better solution to format a list into string
    :param table:
    :return:
    '''
    try:
        output='['+str(table[0])
        for i in range(len(table)-1):
            output+=(','+str(table[i+1]))
        output+=']'
    except:
        raise TypeError('This object is not iterable!')
    return output





def do_one_pdb(pdb,filename=None,REPORTCSV=None,index=0):
    '''
    For each target-complex pdb , this program check if .pdb file exists
    if not ,download first then call the function to match all possible target-ligands with
    molecules in sdf files in one single pdb
    :param pdb:name
    :parameter REPORTCSV: sometimes generate a list of report with the filename this one
    :return:
    '''

    # For each pdb name , there should be one corresponding molecule files.
    # This will generate one result file.
    pdb = pdb.lower()
    if filename is None:
        filename = os.path.join(pdb_PREFIX,'{}.pdb.gz'.format(pdb))
    if os.path.exists(filename):
        # pdbfile exists
        logging.info(pdb + ' has already exists')
        return mol_ligand_tar_generator(pdb,filename,statistic_csv=REPORTCSV,fileforbabel='{}.sdf'.format(index))

    else:
        # Not exists, download from the internet
        urllib.urlretrieve(url_prefix + '{}.pdb.gz'.format(pdb.lower()), filename)
        # Wait for 1 second from rejection on connection.
        time.sleep(1)

        # This is to check whether we download the file successfully
        o = open(filename, 'r')
        for l in o:
            if l.find('DOCTYPE') != -1:
                print 'download {} failed'.format(pdb)
                logging.error('download {} failed'.format(pdb))
                return False
            else:
                #If we download files successfully, then we will run the program
                print 'download {} successfully'.format(pdb)
                logging.info('download {} successfully'.format(pdb))
                return mol_ligand_tar_generator(pdb,filename,statistic_csv=REPORTCSV,fileforbabel='{}.sdf'.format(index))
        o.close()

def initiate_report():
    csv_name = 'report.csv'
    writer = file(csv_name, 'wb')
    w = csv.writer(writer)
    w.writerow(['filename','pdb Name','molecules','paired','bad_one','pairtimes'])
    return csv_name

def quick_split(pdb):
    pdb = pdb.lower()
    fake_pdb_container(pdb,filepos=os.path.join(pdb_PREFIX,pdb+'.pdb.gz'))


if __name__ == '__main__':

    DONE=[]
    FAIL=[]
    ct=0
    report = initiate_report()

    for pdb in PDB_tar[0:1]:
        #dirty way to do small scale tests
        #Use a count variable
        pdb =pdb.lower()
        #real_dir = repair_pdbfile(os.path.join(pdb_PREFIX,'{}.pdb.gz'.format(pdb)),pdb)

        if do_one_pdb(pdb, REPORTCSV=report):
            DONE.append(pdb)
        else:
            FAIL.append(pdb)
        ct+=1

    print ct
    logging.info('total: {}'.format(ct))
