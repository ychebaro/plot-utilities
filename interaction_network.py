#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import argparse
import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
import csv
import networkx as nx

from math     import sqrt
from sys      import argv,stdout,stdin 
from Bio.PDB import NeighborSearch, PDBParser, Selection
from collections import OrderedDict


np.set_printoptions(precision=3, threshold=np.nan)

__author__ = "Yasmine Chebaro"
__email__ = "yasmine.chebaro@gmail.com"

""" Extract interaction network around a specified residue of interest from a   """
""" contact probilibility matrix or from the average probabilities if several   """
""" matrices are available (from different simulations for example.             """
"""                                                                             """
""" For an illustration of the representation, you can have a look at the file  """
""" interaction_network_example.png                                             """


def get_nonzero_contacts(data):
    """ Returns a matrix where:
    - null values in row and column are removed
    - contacts between i, i+1 and i+2 are removed """
    
    nonzero_contacts = OrderedDict()
    for row in range(len(data)):
        for col in range(len(data[0])):
            if data[row,col] != 0.00 and (row not in {col,col-1,col+1,col-2,col+2}):
                nonzero_contacts[str(row+1)+','+str(col+1)] = data[row,col]
                
    return nonzero_contacts


def get_average_matrix(listdata, nb_resids):
    """ Returns average of matrices in list """

    filenames = []
    [filenames.append("data"+str(i)) for i in range(1,len(listdata)+1)]
    
    print filenames
    print len(filenames)
    
    matrix_all = np.zeros((nb_resids, nb_resids))
    for i in range(len(filenames)):
        print filenames[i]
        filenames[i] = np.genfromtxt(listdata[i],delimiter=' ')
        matrix_all += filenames[i]
    
    return np.array(matrix_all/len(filenames))       
            

def get_around_residue(listresids,pdb,chainid,cutoff):
    """ Returns list of residues within cutoff distance of the one specified, and the number of residues in the chain_id specified """
    
    structure = PDBParser(QUIET=True).get_structure('X',pdb)
    chain = structure[0][str(chainid)]
    center_residues = [chain[resi] for resi in listresids]

    center_atoms = Selection.unfold_entities(center_residues,str(chainid))
    atom_list = [atom for atom in structure[0][str(chainid)].get_atoms() if atom.name == 'CA']

    ns = NeighborSearch(atom_list)
    nearby_residues = {res for center_atom in center_atoms
                   for res in ns.search(center_atom.coord, float(cutoff), 'R')}
    
    return sorted(res.id[1] for res in nearby_residues), len(atom_list)



def get_around_contacts(around_residues,data):
    """ Extracts residues within cutoff from the contact probability matrix """
    around_contacts = []
    for key in data.keys():
        if (map(int,string.split(key,','))[0] or map(int,string.split(key,','))[1]) in around_residues:
            around_contacts.append([(map(int,string.split(key,',')))[0],(map(int,string.split(key,',')))[1],data[key]])
    
    return around_contacts


def parse_options():
    parser = argparse.ArgumentParser(description="Plots network of interaction around a specified residue using list of contact matrices (if you only have one just put one matrix in your input list). If several matrices are provided the average contact probability will be calculated.")
    parser.add_argument("-p", "--pdb", dest="pdb_file", required=True,
                        action="store", type=str,
                        help="pdb input file")
    parser.add_argument("-l", "--matlist", dest="matrices_list", required=True,
                        action="store", type=str,
                        help="list of matrix or matrices")
    parser.add_argument("-cid", "--chainid", dest="chainid", required=True,
                        action="store", type=str,
                        help="chain id of the residue considered")
    parser.add_argument("-c", "--cutoff", dest="cutoff", required=True,
                        action="store", type=str, 
                        help="cutoff distance used to find residues around the one considered")
    parser.add_argument("-r", "--resid", dest="resid", required=True,
                        action="store", type=str, 
                        help="number of the residue considered")

    return parser.parse_args()
        
       
def main():
    # get options
    options = parse_options()
    pdb = options.pdb_file
    matlist = options.matrices_list
    chainid = options.chainid
    cutoff = options.cutoff
    resid = options.resid
    
    
    listresids                 = map(int,string.split(resid,','))

    around_residues, nb_resids = get_around_residue(listresids,pdb,chainid,cutoff)
    
    listdata                   = np.genfromtxt(matlist,dtype='str')
    data                       = get_average_matrix(listdata, nb_resids)

    nonzero_data               = get_nonzero_contacts(data)
    around_contacts            = get_around_contacts(around_residues,nonzero_data)    
    
    G = nx.Graph()
    for elm in around_contacts:
        if elm[0] == int(resid) or elm[1] == int(resid):
            G.add_edge(elm[0],elm[1],color='deeppink',weight=4*elm[2])
        else:
            G.add_edge(elm[0],elm[1],color='deepskyblue',weight=4*elm[2])

          
    pos = nx.circular_layout(G)
    
    edges = G.edges()
    colors = [G[u][v]['color'] for u,v in edges]
    weights = [G[u][v]['weight'] for u,v in edges]

    nx.draw(G, pos, edges=edges, node_color='grey', edge_color=colors, width=weights, with_labels=True)    

    plt.show(G)




if __name__=="__main__":
    main()
