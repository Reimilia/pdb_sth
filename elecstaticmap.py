'''
Generate autovina score and elecgrid files. Also differentiate
'''
from mapping import *
from Config import *
import os
from vector_gen import pdb_container
import csv

SUMMARY_COLUMN = ['PDB name','PDB type', 'ligand NAME', 'ligand index in PDB', 'vina score', 'box_scale', 'pure_protein_gridmap_filename',
                  'pure_ligand_gridmap_filename', 'ligand_receptor_complex_gridmap_filename']
PAIR_SUMMARY = 'lignad-receptor_pair.csv'
RESERVE_NAME = ['fake-ligand.pdb']


def generate_one_map(PDBname, BOX=20):
    # mapping files will be generated in this folder (for each pdb file)
    set_new_folder(PDBname,result_PREFIX)

    PDBIndex = pdb_container(PDBname,filepos=pdb_PREFIX+PDBname.lower()+'.pdb.gz')
    #PDBIndex.set_all_vina_benchmark(Box=BOX)
    PDBtype = PDBIndex.get_pdb_type()

    fake_ligand_filename = os.path.join(temp_pdb_PREFIX,'fake-ligand.pdb')

    writer = file('result/{}'.format(PAIR_SUMMARY), 'a')
    w = csv.writer(writer)

    for k,v in PDBIndex.heterodict.items():

        #Detect source file position
        ligand_filename = os.path.join(temp_pdb_PREFIX,v['naming']+'_ligand.pdb')
        receptor_filename = os.path.join(temp_pdb_PREFIX,v['naming']+'_receptor.pdb')
        complex_filename = os.path.join(temp_pdb_PREFIX, v['naming']+'_complex.pdb')
        #prepare auto grid map files
        do_auto_grid(receptor_filename, fake_ligand_filename, center=v['center'])
        do_auto_grid(ligand_filename, fake_ligand_filename, center=v['center'])
        do_auto_grid(complex_filename, fake_ligand_filename, center=v['center'])
        score=do_auto_vina_score(receptor_filename,ligand_filename,v['center'],Box=20)

        #Generate one data
        one_line=[PDBname,PDBtype, v['ligand'].getResname(), k, score, BOX, v['naming']+'_receptor', v['naming']+'_ligand', v['naming']+'_complex']
        w.writerow(one_line)
    writer.flush()
    writer.close()

    #Do this with the risk
    clean_temp_data()




def initialize_summary_file(filename):
    csv_name = filename
    writer = file('result/'+csv_name, 'wb')
    w = csv.writer(writer)
    w.writerow(SUMMARY_COLUMN)
    writer.close()

def clean_temp_data():
    '''
    Warning! This will wipe out all file in ./data except reserved names if not locked by root
    :return:
    '''
    files = os.listdir('data/')
    for filename in files:
        loc = os.path.join('data/') + filename
        if os.path.exists(loc) and filename not in RESERVE_NAME:
            os.remove(loc)


if __name__=='__main__':
    test =PDB_tar[0:5]
    print test

    initialize_summary_file(PAIR_SUMMARY)

    for each in test:
        each= each.lower()
        generate_one_map(each)
