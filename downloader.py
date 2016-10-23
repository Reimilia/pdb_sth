import sys,os
import time
import urllib
from Config import PDB_tar

url_prefix = 'https://files.rcsb.org/download/'
filepath = 'pdb_raw'

#This is a single script to download pdb source files


def down_one(line):
    filename=filepath+'/{}.pdb.gz'.format(line.lower())
    if not os.path.exists(filename):
        time.sleep(2)
        urllib.urlretrieve(url_prefix+'{}.pdb.gz'.format(line.lower()), filename)
        o=open(filename,'r')
        for l in o:
            if l.find('DOCTYPE')!=-1:
                print 'download {} failed'.format(line)
                break
            else:
                print 'download {} successfully'.format(line)
                break
        o.close()
    else:
        print '{} already exists'.format(line)

