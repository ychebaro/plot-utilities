#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator 
from matplotlib.patches import Rectangle


__author__ = "Yasmine Chebaro"
__email__ = "yasmine.chebaro@gmail.com"

""" Plot some data with visualization of secondary structure """
""" uses dssp to get info on secondary structure             """
""" For info:                                                """
""" - H = alpha-helix                                        """ 
""" - B = residue in isolated beta-bridge                    """ 
""" - E = extended strand, participates in beta ladder       """ 
""" - G = 3-helix (310 helix)                                """ 
""" - I = 5 helix (pi-helix)                                 """
""" - T = hydrogen bonded turn                               """
""" - S = bend                                               """
"""                                                          """
""" User-defined selection of secondary structures to plot   """
""" DSSP reference:                                          """
""" W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637  """



def getsecstr(pdbfile,chainid):
    """ Run dssp on pdbfile 
    Extract secondary structure information
    of the pdb chain 
    """
    
    # get dssp output
    os.system("dssp %s > dssp.out" %(pdbfile))
    
    # create lists for sec str elements
    H = []
    B = []
    E = []
    G = []
    I = []
    T = []
    S = []
    
    # extract secondary structure info for the required chain
    with open("dssp.out") as myFile:
        lines = myFile.readlines()
        for line in lines:
            if line[11:12] == str(chainid):
                if line[16:18] == "H ":
                    H.append(int(line[0:5]))
                if line[16:18] == "B ":
                    B.append(int(line[0:5]))
                if line[16:18] == "E ":
                    E.append(int(line[0:5]))
                if line[16:18] == "G ":
                    G.append(int(line[0:5]))
                if line[16:18] == "I ":
                    I.append(int(line[0:5]))
                if line[16:18] == "T ":
                    T.append(int(line[0:5]))
                if line[16:18] == "S ":
                    S.append(int(line[0:5]))
                    
    return H,B,E,G,I,T,S

def consecutive(data, stepsize=1):
    """ split array or list into arrays of consecutive numbers """
    return np.split(data,np.where(np.diff(data) != stepsize)[0]+1)
      
        
def plotdatawithss(inputfile,listofss,H,B,E,G,I,T,S):
    fig = plt.figure(1, figsize=(10,4))
    ax1 = fig.add_subplot(111)
    
    # read and plot data
    x,y = np.genfromtxt(inputfile,unpack=True,dtype='float')
    ax1.plot(x,y,color='black',linewidth=1.5)  # adjust as see fit
    
    # labels and co
    ax1.set_xlabel('Residue number')
    ax1.set_ylabel('RMSF ($\AA$)')
    
    # add rectangles for the secondary elements selected
    yss = 0.8*min(y)
    rwidth = 0.1*min(y)
        
    for ss in listofss:
        if ss == "H":
            for elm in consecutive(H):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='firebrick'))
        if ss == "B":
            for elm in consecutive(B):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='mediumspringgreen'))
        if ss == "E":
            for elm in consecutive(E):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='forestgreen'))
        if ss == "G":
            for elm in consecutive(G):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='darkviolet'))
        if ss == "I":
            for elm in consecutive(I):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='deeppink'))
        if ss == "T":
            for elm in consecutive(T):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='deepskyblue'))
        if ss == "S":
            for elm in consecutive(S):
                ax1.add_patch(Rectangle((elm[0],yss),len(elm),rwidth,color='dodgerblue'))

    plt.show()
    

def parse_options():
    parser = argparse.ArgumentParser(description='Plot secondary structure elements as rectangles in a figure, runs dssp to get info. Argument -s (--sel) consists of a list of secondary structure elements chosen to be represented (H,B,E,G,I,T or S), format is ["E","G"] for example.')
    parser.add_argument("-f","--file", dest="datafile", required=True,
                        action='store', type=str,
                        help='data file, in (x,y) format')
    parser.add_argument("-p","--pdb", dest='pdbfile', required=True,
                        action='store', type=str,
                        help='pdb file')
    parser.add_argument("-c","--chain", dest='chainid', required=True,
                        action='store', type=str,
                        help='chain id for secondary structure representation')
    parser.add_argument("-s","--sel", dest='selss', required=True,
                        action='store', type=str,
                        help='selection of secondary structure elements chosen for representation')
    
    return parser.parse_args()
    
def main():
    # get options
    options   = parse_options()
    data      = options.datafile
    pdb       = options.pdbfile
    chain     = options.chainid
    selection = options.selss
    
    H,B,E,G,I,T,S = getsecstr(pdb, chain)
    plotdatawithss(data,selection,H,B,E,G,I,T,S)
    
    
if __name__ == "__main__":
    main()

    
    