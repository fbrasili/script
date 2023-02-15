import os
import sys
from numpy import *
from configurations import * 
#from init_file import *
#from utils import *

def compute_bonded_per_np_singleConf( data , threshold ):
    '''
    given the data of a single configuration compute the number of charged monomers in interaction with each np according to the criterium:
       one monomer is in interaction with one nanoparticle if their surfece-to-surface distance is lower than the threshold
    data is a nAtoms x 5 array, listing for each atom type, pos_x, pos_y, pos_z
    return a list of N integer values (number of interacgting monomers), where N is the number of nanoparticles
    '''
    cond_ions = data[:, 0].astype(int) == 3
    cond_nps  = data[:, 0].astype(int) == 4
    data_ions = data[cond_ions, 1:]
    data_nps  = data[cond_nps , 1:]

    return array([ sum( linalg.norm( data_ions - datum, axis=1) <= threshold ) for datum in data_nps ])
        

def compute_bonded_ions_np(confs_fname, threshold):
    '''
    read the configurations in the file confs_fname and for each compute:
        1. the distribution of the number of charged monomers bonded to each nanoparticle
        2. the number Nads of nanoparticles adsorbed to the microgel (nanopartilces bonded to at least one charged monomer)
    return the total distribution of the number of charged monomers bonded to each nanoparticle, the list of timesteps and the list of Nads
    '''
    nAtoms,first,interval,linesXtimestep,last,_ = readConfigurations(filename=confs_fname, retVals=True, check=False, verb=False)
    bonded_per_np,adsorbed_nps = [],[]
    for skiprows in linesXtimestep * arange( 1 + ( last - first ) / interval , dtype="int") + 9:
        data = loadConfiguration( confs_fname, skiprows=skiprows, max_rows=nAtoms, usecols=(1,2,3,4) )
        bonded_per_np_singleConf = array( compute_bonded_per_np_singleConf( data , threshold ) )
        bonded_per_np += [ bonded_per_np_singleConf ]
        adsorbed_nps  += [ (bonded_per_np_singleConf > 0).sum() ]
    return array(bonded_per_np).flatten(), arange(first, last+interval, interval), array(adsorbed_nps)
