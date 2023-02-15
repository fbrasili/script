import os
import sys
from numpy import *
from configurations import * 
from init_file import *
from utils import *

##########################################################################
### GENERATE A LIST OF CHAINS (WITH CROSSLINKERS OR DANGLING ENDS AS ENDS)

def gen_chains(init_fname, rem_loops=True, rem_dangling=False, saveChains=True, ret=True):
    fol = os.path.dirname( init_fname )

    nAtoms,nBonds,_,idxAtoms,idxBonds = read_init( init_fname )
    bonds = loadtxt( init_fname, skiprows=idxBonds, max_rows=nBonds, usecols=(2,3), dtype=int )
    atoms = loadtxt( init_fname, skiprows=idxAtoms, max_rows=nAtoms, usecols=(0,1), dtype=int )
    atomCounter = Counter( bonds.flatten() )
    ends = [ key  for key,value in atomCounter.items() if value != 2 ]

    chains = []
    freeEnds = {}
    for a1,a2 in bonds:
        isEnd1 = a1 in ends
        isEnd2 = a2 in ends
        isInChain1 = a1 in freeEnds.keys()
        isInChain2 = a2 in freeEnds.keys()
        if isEnd1 and isEnd2:
            chains += [ [a1,a2] ]
        elif isEnd1 or isEnd2:
            if isEnd2: a1,a2,isInChain2 = a2,a1,isInChain1
            if isInChain2:
                idxChain = freeEnds[a2]
                chain = chains[ idxChain ].copy()
                if chain[0] in ends:
                    chain = chain + [a1]
                else:
                    if chain.index(a2) == 0: chain = [a1] + chain
                    else                   : chain = [a1] + chain[::-1]
                chains[idxChain] = chain
                del freeEnds[a2]
            else:
                freeEnds[a2] = len( chains )
                chains += [ [a1,a2] ]
        else:
            if isInChain1 and isInChain2:
                idx1,idx2 = freeEnds[a1],freeEnds[a2]
                chain1,chain2 = chains[idx1].copy(),chains[idx2].copy()
                if chain1.index(a1) == 0: chain1 = chain1[::-1]
                if chain2.index(a2) != 0: chain2 = chain2[::-1]
                idx_min,idx_max = min(idx1,idx2),max(idx1,idx2)
                chains[ idx_min ] = chain1 + chain2
                del chains[idx_max]
                del freeEnds[a1]
                del freeEnds[a2]
                freeEnds.update( { key: idx_min   for key,value in freeEnds.items() if value == idx_max } )
                freeEnds.update( { key: value - 1 for key,value in freeEnds.items() if value >  idx_max } )
            elif isInChain1 or isInChain2:
                if isInChain2: a1,a2 = a2,a1
                idxChain = freeEnds[a1]
                chain = chains[idxChain].copy()
                if chain.index(a1) == 0: chain = chain[::-1]
                chains[idxChain] = chain + [a2]
                del freeEnds[a1]
                freeEnds[a2] = idxChain
            else:
                freeEnds[a1],freeEnds[a2] = len( chains ),len( chains )
                chains += [ [a1,a2] ]

    chains = array(chains, dtype=list)

    n0,n = len(chains),len(chains)
    print( 'total number of chains: %i'%n0 )

    if rem_loops: 
        chains = array( [chain for chain in chains if chain[0] != chain[-1] ] , dtype=list )
        print( 'discarded %i loops'%(n0-len(chains)) )
        n = len(chains)
    if rem_dangling:
        danglingAtoms = [ key for key,value in atomCounter.items() if value == 1 ]
        chains = array( [ chain for chain in chains if not ( chain[0] in danglingAtoms or chain[-1] in danglingAtoms ) ] , dtype=list )
        print( 'discarded %i dangling chains'%(n-len(chains)) )
        n = len(chains)

    if n < n0: print('final number of chains: %i'%n )

    chains_woCL,crosslinkers = [],[ key  for key,value in atomCounter.items() if value > 2 ]
    for chain in chains:
        if chain[ 0] in crosslinkers: chain = chain[1:  ]
        if chain[-1] in crosslinkers: chain = chain[ :-1]
        chains_woCL += [chain]
    chains = array(chains_woCL, dtype=list)

    N = array( [ len(chain) for chain in chains ], dtype=int )
    N,chains = N[argsort(N)][::-1],chains[argsort(N)][::-1]
    
    if saveChains:
        with open( os.path.join( fol , 'chains.dat' ) , "w" ) as of:
            for chain in chains:   of.write( '\t'.join( [ str(atom) for atom in chain ] ) + '\n' )
    if ret:
        return N,chains


#########################################################################
### MODIFY THE PARTICLE TYPES IN CONFIGURATION FILE TO DISTINGUISH CHAINS 

def modifyTypes(conf_file, mode='chains', chains_file=None, len_min=None, len_max=None, outfile='configuration_chains.dat'):
    '''
    conf_file is a file with only one configuration
    returna a file with the same configuration and atom types modified according to mode

    mode = 'chains'
        change the type of neutral monomers to identify the different chains
        (charged monomers, crosslinkers and other atoms are not considered)
    '''
    assert os.path.isfile(chains_file), "missing file " + chains_file + " with listed chains"
    if mode == 'chains':
        # read the chains from file chains_file
        with open(chains_file, 'r' ) as f:
            chains = array( [ [ int( atom ) for atom in chain.split() ] for chain in f] , dtype=list )
        N = array( [ len(chain) for chain in chains ], dtype=int )
        
        # removes chains with length outside len_range and excludes atoms at the ends of chains (crosslinkers)
        if len_min is None: len_min = N.min()
        if len_max is None: len_max = N.max()
        cond = (N >= len_min) * (N <= len_max)
        chains = array( [ chain[1:-1] for chain in chains[cond] ] , dtype=list )
        
        # creates the file outfile with modified types
        of = open( outfile , "w" )
        print("writing configuration file " + outfile + " ...")
        with open( conf_file, "r" ) as f:
            for i,line in enumerate(f):
                line = line.split()
                if i>8 and int( line[1] ) == 1:
                    atom = int( line[0] )
                    for ID,chain in enumerate(chains):
                        if atom in chain: line[1] = str(ID+10)
                line = " ".join(line) + '\n'
                of.write(line)
        of.close()


####################################################
### COMPARE THE CHAINS LISTED IN TWO DIFFERENT FILES

def compare_chains(chains_file1, chains_file2, verb=True, ret=False):
    if verb: print( "comparing chains in " + chains_file1 + " with chains in " + chains_file2 + " ..." )
    with open( chains_file1 , "r" ) as chf:
        chains1 = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
    with open( chains_file2 , "r" ) as chf:
        chains2 = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
    if verb:
        print( "number of chains in file " + chains_file1 + ": %i"%len(chains1) )
        print( "number of chains in file " + chains_file2 + ": %i"%len(chains2) )
    
    if all([ any([ chain1 == chain2 or chain1 == chain2[::-1] for chain1 in chains1 ]) for chain2 in chains2 ]) and all([ any([ chain2 == chain1 or chain2 == chain1[::-1] for chain2 in chains2 ]) for chain1 in chains1 ]):
        if verb: print( "the same chains are listed in the two files")
        if ret: return True
    else:
        print( "some bonds are different" )
        print()
        if verb:
            print( "the following chains from file " + chains_file2 + " are not present in " + chains_file1 )
            for chain2 in chains2:
                if all([ chain1 != chain2 and chain1 != chain2[::-1]  for chain1 in chains1 ]): print( chain2 )
            print()
            print( "the following chains from file " + chains_file1 + " are not present in " + chains_file2 )
            for chain1 in chains1:
                if all([ chain2 != chain1 and chain2 != chain1[::-1]  for chain2 in chains2 ]): print( chain1 )
        if ret: return False


###############################################
### FUNCTIONS FOR ANALYSIS OF CHAINS PROPERTIES

### functions for single chain
def calcCharge_singleChain(chain, configuration_data):
    '''
    chain is a list of int values, corresponding to the atom ID in che configuration_data
    configuration_data is an array with shape (N_monomers, 6), providing for each monomer ID, atom type, x, y, z and q 
    
    return the sum of the charges of the monomers in the chain
    '''
    idxs = array( [ where( configuration_data[:,0]==atom ) for atom in chain ] ).flatten()
    return configuration_data[idxs , 5].sum()
    
def calcRee_singleChain(chain, configuration_data):
    '''
    chain is a list of int values, corresponding to the atom ID in che configuration_data
    configuration_data is an array with shape (N_monomers, 6), providing for each monomer ID, atom type, x, y, z and q

    return the end-to-end distance of the chain
    '''
    if len(chain) == 1:
        return 0
    else:
        position1 = configuration_data[configuration_data[:,0]==chain[ 0], 2:5]
        position2 = configuration_data[configuration_data[:,0]==chain[-1], 2:5]
        return linalg.norm( position1 - position2 )

def calcRg_singleChain(chain, configuration_data):
    '''
    chain is a list of int values, corresponding to the atom ID in che configuration_data
    configuration_data is an array with shape (N_monomers, 6), providing for each monomer ID, atom type, x, y, z and q

    return the gyration radius of the chain
    '''
    if len(chain) == 1:
        return 0
    else:
        idxs = array( [ where( configuration_data[:,0]==atom ) for atom in chain ] ).flatten()
        positions = configuration_data[idxs ,2:5]
        return linalg.norm( positions - positions.mean(0) ) / sqrt( len(chain) )

def calcDistFromMgelCM_singleChain(chain, configuration_data):
    '''
    chain is a list of int values, corresponding to the atom ID in che configuration_data
    configuration_data is an array with shape (N_atoms, 6), providing for each atom ID, atom type, x, y, z and q

    return the distance between the center of mass of the chain and that of the microgel
    '''
    # position of the microgel center of mass
    cond_mgel = ( configuration_data[:,1] == 1 ) + ( configuration_data[:,1] == 3 ) 
    com_mgel  = configuration_data[cond_mgel, 2:5].mean(0)
    
    # positions of the monomers in the chain
    idxs_chain = array( [ where( configuration_data[:,0]==atom ) for atom in chain ] ).flatten()
    #chain_data = configuration_data[idxs_chain, 2:5]
    com_chain = configuration_data[idxs_chain, 2:5].mean(0)
    # chain distance from microgel center of mass
    #if   mode == 'com'    : return linalg.norm( chain_data.mean(0) - com_mgel )
    #elif mode == 'average':
    return linalg.norm( com_chain - com_mgel ) 

def calcMinDistFromNps_singleChain(chain, configuration_data, mode):
    '''
    chain is a list of int values, corresponding to the atom ID in che configuration_data
    configuration_data is an array with shape (N_atoms, 6), providing for each atom ID, atom type, x, y, z and q

    calculate the distance of the chain from each nanoparticle, where the chain position is determined according to mode:
      'charges' : uses the position of charged monomer
      'monomers': uses the position of all monomers

    return the minimum chain-np distance
    '''
    assert mode in ['charges','monomers'], 'allowed mode values are \'charges\' (positions of charged monomer) and \'monomers\' (positins of all monomers)'
 
    # configuration of the monomers in the chain
    idxs_chain = array( [ where( configuration_data[:,0]==atom ) for atom in chain ] ).flatten()
    chain_data = configuration_data[idxs_chain]
    
    # positions of the nanoparticles
    nps_positions = configuration_data[configuration_data[:,1]==4, 2:5]

    # chain charge
    q = chain_data[:,  -1 ].sum()
    
    if mode == 'charges' and q == 0:
        return nan
    else:
        if   mode == 'charges':
            chain_positions = chain_data[ chain_data[:,1]==3 , 2:5 ]
        elif mode == 'monomers':
            chain_positions = chain_data[ :                  , 2:5 ]
        return min( [ linalg.norm(nps_positions - pos, axis=1).min() for pos in chain_positions ] )
    
### functions for single configuration
def calcCharge_singleConf(chains_data, conf_data):
    '''
    chain is a list of N chains (each chain is a list of atoms)
    conf_data is an array with shape (N_monomers, 6), providing for each monomer ID, atom type, x, y, z and q
    
    return an array of N charge values, one for each chain
    '''
    return  array( [ calcCharge_singleChain(chain, conf_data) for chain in chains_data ] )

def calcRees_singleConf(chains_data, conf_data):
    '''
    chains_data is a list of N chains (each chain is a list of atoms)
    conf_data is an array with shape (N_monomers, 6), providing for each monomer the ID, atom type, x, y, z and q
    
    return an array of N end-to-end distances, one for each chain
    '''
    return array( [ calcRee_singleChain(chain, conf_data) for chain in chains_data ] )

def calcRgs_singleConf(chains_data, conf_data):
    '''
    chains_data is a list of N chains (each chain is a list of atoms)
    conf_data is an array with shape (N_monomers, 5), providing for each monomer the ID, atom type, x, y and z positions

    return an array of N gyration radii, one for each chain
    '''
    return array( [ calcRg_singleChain(chain, conf_data) for chain in chains_data ] )

def calcDistsFromMgelCM_singleConf(chains_data, conf_data):
    '''
    chains_data is a list of N chains (each chain is a list of atoms)
    conf_data is an array with shape (N_monomers, 5), providing for each monomer the ID, atom type, x, y and z positions

    return an array of N distances, one for each chain
    '''
    return array( [ calcDistFromMgelCM_singleChain(chain, conf_data) for chain in chains_data ] )

def calcDistsFromNps_singleConf(chains_data, conf_data, mode, verb=False):
    '''
    chains_data is a list of N chains (each chain is a list of atoms)
    conf_data is an array with shape (N_monomers, 5), providing for each monomer the ID, atom type, x, y and z positions

    return an array of N distances, one for each chain
    '''
    if verb:
        if   mode == 'charges' : print( "the chain-nanoparticle distance is calculated using the positions of charged monomers" )
        elif mode == 'com'     : print( "the chain-nanoparticle distance is calculated using the chain center-of-mass" )
        elif mode == 'monomers': print( "the chain-nanoparticle distance is calculated using the positions of all monomers" )
    return array( [ calcMinDistFromNps_singleChain(chain, conf_data, mode) for chain in chains_data ] )


### functions for all configurations
def compute_ChainsProperties(chains_data, confs_fname, lMin=None, computeDistsFromNps=False):
    '''
    chains_data is a list of N chains (each chain is a list of atoms)
    confs_fname is the name of the configurations file
    
    analyze only chain with lengh (number of monomers) equal or larger than lMin (if lMin is None, lMin = 2)
    
    for each configuration, calculate:
     1. length of chains
     2. charge of chains
     3. distances of chains from the center of mass of the microgel according to mode:
        'com': uses the center of mass of the chain
        'average': uses the position of each monomer and compute the average distance

     4. gyration radii and end-to-end distances of chains
     5. chains distance from nanoparticles (only if computeDistsFromNps is True): shorter monomer-nanoparticle distance
     6. chains distance from nanoparticles (only if computeDistsFromNps is True): shorter charge-nanoparticle distance

    return the computed quantities
    '''
    #  number of atoms and configurations
    nAtoms,first,step,linesXstep,last,_ = readConfigurations( filename=confs_fname, retVals=True, check=False, verb=False )
    nConfs = int( (last - first) / step ) + 1
    
    # chains lenghts and minimum number of monomers
    L = array( [ len(chain) for chain in chains_data ] ) 
    
    if lMin is None: lMin = 1
    print( "selecting chains with at least %i monomers..."%lMin )
    cond = L >= lMin
    print( "total number of chains: %i"%len(L[cond]) )
    #L = tile(L[cond] , nConfs ) 
    L = tile(L[cond] , (nConfs,1) ) 
    chains_data = chains_data[cond]
    
    # chains charge, radial position, Ree and Rg 
    print() ; print("computing chains properties ...")
    qs, rs, Rees, Rgs = [], [], [], []
    if computeDistsFromNps: ds,dcs = [], []
    for i in arange(nConfs):
        print('reading configuration %i'%(first + i*step) + ' of %i ...'%last)
        conf_data = loadConfiguration( confs_fname, skiprows=(linesXstep*(i+1) - nAtoms), max_rows=nAtoms, usecols=(0,1,2,3,4,5) )
        qs   += [ calcCharge_singleConf(          chains_data , conf_data ) ]
        rs   += [ calcDistsFromMgelCM_singleConf( chains_data , conf_data ) ]
        Rees += [ calcRees_singleConf(            chains_data , conf_data ) ]
        Rgs  += [ calcRgs_singleConf(             chains_data , conf_data ) ]
        if computeDistsFromNps:
            ds  += [ calcDistsFromNps_singleConf( chains_data , conf_data, 'monomers' ) ]
            dcs += [ calcDistsFromNps_singleConf( chains_data , conf_data, 'charges'  ) ]
        del conf_data
    print( )
    
    qs, rs, Rees, Rgs = array(qs), array(rs), array(Rees), array(Rgs)
    
    if computeDistsFromNps:
        ds,dcs = array(ds), array(dcs)
        return L,qs,rs,Rgs,Rees,ds,dcs
    else:
        return L,qs,rs,Rgs,Rees
    
##############################################
### FUNCTIONS FOR ANALYSIS OF CHAINS SHRINKING

def add_columnRg0(path_data, path_data_ref, ret=False):
    '''
    read reference data of gyration radius and radial position from file path_data_ref
    add them as new columns in the file path_data
    '''
    # read data from path_data_ref
    L0,q0,r0,Rgs0 = loadtxt(path_data_ref, skiprows=1, unpack=True, usecols=(0,1,2,3))
    # read data from path_data
    data    = loadtxt(path_data    , skiprows=1, unpack=True)
    # discard any columns with ref data
    if len(data) > 7: data = data[:7]
    
    # check that the number of chains in path_data is a multople of those in path_data_ref (repeated configurations)
    assert shape(data)[1]%len(L0) == 0, "the number of chain in " + path_data + " is not consistent with the number of chains in " + path_data_ref
    nConfs = int( shape(data)[1]/len(L0) )
    
    # the number of columns in the file path_data must be higher or equal to 7
    #  7: chains properties (included chains-np distance)
    # >7: chains properties (included chains-np distance) + other variables from ref data (remove the columns and replace them with new ones)
    assert len(data) == 7, "wrong number of columns in file " + p + "\n run again compute_chains_properties"
    
    # heading of the new file
    heading = "# L, q, r (distance from mgel COM), Rg, Ree, d (shorter monomer-NP distance), dc (shorter charge-NP distance), <r0>, <Rg0> \n"

    # append the new column with values of Rgs0 to data 
    data =  append( data , [ tile(r0, nConfs) , tile(Rgs0, nConfs) ] , axis=0 )
    
    # write the file path_data with the new column
    temp_file = os.path.join( os.path.dirname( path_data ) , ".temp_data" )
    os.system( "rm " + path_data )
    with open( path_data , "w" ) as f: f.write( heading )
    savetxt( temp_file , data.transpose() )
    os.system( "cat " + temp_file + " >> " + path_data )
    os.system( "rm " + temp_file )
    print( "file " + path_data + " updated" )

    if ret: return L0,q0,r0,Rgs0

def compute_shrinkingFromChains(chains_properties_fname, usecols=(0,1,2,3,7,8)):
    L,q,r,Rg,r0,Rg0 = loadtxt(chains_properties_fname, skiprows=1, usecols=usecols, unpack=True)
    Nchains = len(L)
    
    neutral = q == 0
    charged = q >  0
    
    Rm,Rm0 = average( Rg , weights=L ) , average( Rg0 , weights=L )
    rm,rm0 = average(  r , weights=L ) , average(  r0 , weights=L )
    sR,sR0 = sqrt( ( average( Rg**2 , weights=L ) - Rm**2 ) / Nchains ) , sqrt( ( average( Rg0**2 , weights=L ) - Rm0**2 ) / Nchains )
    sr,sr0 = sqrt( ( average(  r**2 , weights=L ) - rm**2 ) / Nchains ) , sqrt( ( average(  r0**2 , weights=L ) - rm0**2 ) / Nchains )
    shr,dis = Rm/Rm0 , rm/rm0
    errShr =  shr * hypot( sR/Rm , sR0/Rm0 )
    errDis =  dis * hypot( sr/rm , sr0/rm0 )
    
    shr_n = average( Rg[neutral] , weights=L[neutral] ) / average( Rg0[neutral] , weights=L[neutral] )
    shr_c = average( Rg[charged] , weights=L[charged] ) / average( Rg0[charged] , weights=L[charged] )
    dis   = average(  r          , weights=L          ) / average(  r0          , weights=L          )
    dis_n = average(  r[neutral] , weights=L[neutral] ) / average(  r0[neutral] , weights=L[neutral] )
    dis_c = average(  r[charged] , weights=L[charged] ) / average(  r0[charged] , weights=L[charged] )

    A,A0 = average(  r**2 , weights=L ),average(  r0**2 , weights=L )
    B,B0 = average( Rg**2 , weights=L ),average( Rg0**2 , weights=L )
    sqshr = sqrt( A/A0 )
    sqdis = sqrt( B/B0 )
    sqtot = sqrt( (A+B) / (A0+B0) )
    
    return shr, errShr, shr_c, shr_n, dis, errDis, dis_c, dis_n, sqshr, sqdis, sqtot
'''
def calcReeRg_Ave_allConfs(chains_data, chains_charged_data, chains_neutral_data, confs_fname, verb=False):
'''
#    chains_data is a list of N chains (each chain is a list of atoms)
#    confs_fname is the name of the configurations file
#
#    for each configuration calculate the Rees and Rgs of chains
#
#    return the average Rees and Rgs and the RMS average of Rees of chains
'''
    nAtoms,first,step,linesXstep,last,_ = checkConfigurations(filename=confs_fname, ret=True, printInfo=False)
    nConfs = int( (last - first) / step ) + 1
    Rees,Rgs,Rgs_charged,Rgs_neutral = [],[],[],[]
    for i in arange(nConfs):
        if verb: print('reading configuration %i ...'%(first + i*step) + ' of %i'%last)
        conf_data = loadConfiguration(confs_fname, skiprows=(linesXstep*(i+1) - nAtoms), max_rows=nAtoms, usecols=(0,1,2,3,4))
        Rees  += [ calcRees_singleConf( chains_data, conf_data ) ]
        Rgs   += [ calcRgs_singleConf(  chains_data, conf_data ) ]
        Rgs_charged += [ calcRgs_singleConf( chains_charged_data, conf_data ) ]
        Rgs_neutral += [ calcRgs_singleConf( chains_neutral_data, conf_data ) ]
    Rees,Rgs,Rgs_charged,Rgs_neutral = array( Rees ),array( Rgs ),array( Rgs_charged ),array( Rgs_neutral )
    mRees,mRgs,mRgs_charged,mRgs_neutral = Rees.mean(0),Rgs.mean(0),Rgs_charged.mean(0),Rgs_neutral.mean(0)
    rmsRees = linalg.norm(Rees, axis=0) / sqrt( nConfs )
    return rmsRees,mRgs,mRees,mRgs_charged,mRgs_neutral
'''
'''
def compute_DistsFromMgelCM_allConfs(chains_data, chains_charged_data, chains_neutral_data, confs_fname, mode='com', verb=False):
'''
#    chains_data is a list of N chains (each chain is a list of atoms)
#    confs_fname is the name of the configurations file
#    
#    for each configuration, calculate the distance of chains from the center of mass of the microgel according to mode:
#        'com': uses the center of mass of the chain
#        'average': uses the position of each monomer and compute the average distance
#    
#    return the average distance from of chains from the center of mass of the microgel
'''
    nAtoms,first,step,linesXstep,last,_ = checkConfigurations(filename=confs_fname, ret=True, printInfo=False)
    nConfs = int( (last - first) / step ) + 1

    dists,dists_charged,dists_neutral = [],[],[]
    for i in arange(nConfs):
        if verb: print('reading configuration %i'%(first + i*step) + ' of %i ...'%last)
        conf_data = loadConfiguration(confs_fname, skiprows=(linesXstep*(i+1) - nAtoms), max_rows=nAtoms, usecols=(0,1,2,3,4))
        dists += [ calcDistsFromMgelCM_singleConf(chains_data, conf_data, 'com') ]
        dists_charged += [ calcDistsFromMgelCM_singleConf(chains_charged_data, conf_data, 'com') ]
        dists_neutral += [ calcDistsFromMgelCM_singleConf(chains_neutral_data, conf_data, 'com') ]
        del conf_data
    if verb: print( )

    return array( dists ).mean(0),array( dists_charged ).mean(0),array( dists_neutral ).mean(0)
'''
 
def compute_DistsFromNps(chains_data, chains_charged_data, chains_neutral_data, confs_fname, lMin=None):
    '''
    chains_data, chains_charged_data and chains_neutral_data are lists of chains (each chain is a list of atoms)
    confs_fname is the name of the configurations file

    select only chain with lengh (number of monomers) equal or larger than lMin (if lMin is None, lMin = 2)
    
    for each configuration, calculate the distances of chains frome the nearest nanoparticle, where the chain position is determined according to mode:
    'charges' : uses the positions of charged monomers (selecting the shorter distance from nanoparticles)
                only chains with charges are selected in this case
    'monomers': uses the positions of all monomers (selecting the shorter distance from nanoparticles)
    'com'     : uses the center of mass of the chain

    return the computed quantities
    '''
    nAtoms,first,step,linesXstep,last,_ = checkConfigurations(filename=confs_fname, ret=True, printInfo=False)
    nConfs = int( (last - first) / step ) + 1

    L   = array([ len(chain) for chain in chains_data         ])
    L_c = array([ len(chain) for chain in chains_charged_data ])
    L_n = array([ len(chain) for chain in chains_neutral_data ])

    if lMin is None: lMin = 2
    print( "selecting chains with at least %i monomers"%lMin )
    cond_c = L_c >= lMin ; print( "number of charged chains: %i"%len(L_c[cond_c]) ) ; L_c = tile(L_c[cond_c], (nConfs,1)) ; chains_charged_data = chains_charged_data[cond_c]
    cond_n = L_n >= lMin ; print( "number of neutral chains: %i"%len(L_n[cond_n]) ) ; L_n = tile(L_n[cond_n], (nConfs,1)) ; chains_neutral_data = chains_neutral_data[cond_n]
    cond   = L   >= lMin ; print( "total number of chains: %i"%len(L[cond]  )     ) ; L   = tile(L[cond]    , (nConfs,1)) ; chains_data         = chains_data[cond]
    

    print() ; print("computing chains-nanoparticles distances ...")
    dists_charges,dists_monomers,dists_monomers_c,dists_monomers_n = [],[],[],[]
    for i in arange(nConfs):
        print('reading configuration %i'%(first + i*step) + ' of %i ...'%last)
        conf_data = loadConfiguration(confs_fname, skiprows=(linesXstep*(i+1) - nAtoms), max_rows=nAtoms, usecols=(0,1,2,3,4))
        dists_charges    += [ calcDistsFromNps_singleConf(chains_charged_data, conf_data, mode='charges' , verb=False) ]
        dists_monomers   += [ calcDistsFromNps_singleConf(chains_data        , conf_data, mode='monomers', verb=False) ]
        dists_monomers_c += [ calcDistsFromNps_singleConf(chains_charged_data, conf_data, mode='monomers', verb=False) ]
        dists_monomers_n += [ calcDistsFromNps_singleConf(chains_neutral_data, conf_data, mode='monomers', verb=False) ]
        del conf_data
    print( )

    return L,L_c,L_n,array(dists_charges), array(dists_monomers), array(dists_monomers_c), array(dists_monomers_n)

'''
def chains_analysis(chains_file, confs_file, dist_mode, Ree0_file, Rg0_file, chains_charged_file=None, chains_neutral_file=None, makePlots=False, verb=False):
'''
#    chains_file is the name of the file with the list of all the chains in a microgel
#    confs_file is the name of the file with the configurations
#    makePlots is True, plots all the analyses in the directory ./plots
#    
#    ### scrivere cosa fa ###
'''
    fol = os.path.dirname( chains_file )
    if not chains_charged_file is None: calcCharged = True
    else: calcCharged = False
    if not chains_neutral_file is None: calcNeutral = True
    else: calcNeutral = False

    ### read chains from chains_file
    with open( chains_file , "r" ) as chf:
        chains = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
    N = array([ len(chain) for chain in chains ])

    # excludes chains outside range [nMin, nMax]
    if nMin is None: nMin = N.min()
    if nMax is None: nMax = N.max()
    cond = ( N >= nMin ) * ( N <= nMax )
    chains,N = chains[cond],N[cond]


    ### read chains from chains_file
    print("reading chains from file " + chains_file + " ..." )
    with open( chains_file , "r" ) as chf:
        chains = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
    if calcCharged:
        print("reading charged chains from file " + chains_charged_file + " ..." )
        with open( chains_charged_file , "r" ) as chf:
            chargedChains = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
        Rgs_charged = []
    if calcNeutral:
        print("reading neutral chains from file " + chains_neutral_file + " ..." )
        with open( chains_neutral_file , "r" ) as chf:
            neutralChains = array( [ [ int(atom) for atom in line.split() ] for line in chf ], dtype=list )
        Rgs_neutral = []

    ### read single configurations from confs_file and performs the chain analysid for each 
    nAtoms,first,step,linesXtimestep,last,_ = checkConfigurations(filename=confs_file, ret=True, printInfo=False)
    nConfs = int( (last - first) / step ) + 1
    Rees_all,Rgs_all,dists = [],[],[]
    for i in arange(nConfs):
        if verb: print('reading configuration %i'%(first + i*step) + ' of %i'%last)
        conf_data = loadtxt(confs_file, skiprows=(linesXtimestep*(i+1) - nAtoms), max_rows=nAtoms, usecols=(0,1,2,3,4) )
        Rees_all += [ calcRees_singleConf( chains, conf_data ) ]
        Rgs_all  += [ calcRgs_singleConf(  chains, conf_data ) ]
        if calcCharged: Rgs_charged += [ calcRgs_singleConf(  chargedChains, conf_data ) ]
        if calcNeutral: Rgs_neutral += [ calcRgs_singleConf(  neutralChains, conf_data ) ]
        dists_all += [ calcDistsFromNps_singleConf( chains, conf_data, dist_mode ) ]
    Rees_all,Rgs_all,dists_all = array( Rees_all ),array( Rgs_all ),array( dists_all )
    if calcCharged: Rgs_charged = array( Rgs_charged )
    if calcNeutral: Rgs_neutral = array( Rgs_neutral )

    Rees0,Rgs0 = loadtxt(Ree0_file).flatten(),loadtxt(Rg0_file).flatten()
    
    mRees,mRgs = Rees_all.mean(0),Rgs_all.mean(0)
    rmsRees = linalg.norm(Rees_all, axis=0) / sqrt( nConfs )
'''

#####################################################################
# COMPUTE MICROGEL GYRATION RADIUS-RELATED QUANTITIES BASED ON CHAINS



###############################################################
# COMPUTE g(r) BETWEEN CHAINS CENTERS-OF-MASS AND NANOPARTICLES

#funzione per calcolare la gr tra due tipi di particle
def compute_gr_twoTypes(positions1, positions2, bins=100):
    '''
    compute the g(r) between particles of two distinct groups
    the particles positions are given separately for the two groups 
    the positions arrays, positions1 and positions2, have shape (N1,3) and (N2,3)
    '''
    dists = concatenate([linalg.norm( positions1 - point, axis=1 ) for point in positions2])
    g,r = histogram(dists, bins=bins, density=True)
    return array([r[:-1],g/2.])

'''
funzioni specifiche per calcolare la gr tra catene e nanoparticelle
def compute_gr_chainNP_singleChain(chain, atoms, positions_monomers, positions_nps):
    idxs = array( [ where( atoms==atom ) for atom in chain ] ).flatten()
    positions = positions_monomers[idxs]
    cm = positions.mean(1)

def compute_gr_chainNP_allChains( chains, data_monomers, data_nps ):


def compute_gr_chainNP_allConfs(filename, chains, verb=True):
    nAtoms,first,step,linesXtimestep,last,_ = checkConfigurations(filename=filename, ret=True, printInfo=False)
    headingLines = linesXtimestep - nAtoms
    nConfs = int( (last - first) / step ) + 1
    coms_all,grs = [],[]
    for i in arange(nconfs):
        if verb: print('reading configuration %i'%(first + i*step) + ' of %i'%last, end='\r')
        skiprows = headingLines + i * linesXtimestep
        data = loadtxt(filename, skiprows=skiprows, max_rows=nAtoms, usecols=(0,2,3,4), unpack=True )
        nps = data[0] == 4 ; momomers = (data[0] == 1) or (data[0] == 3)
        data_nps      = data[1: , nps     ]
        data_monomers = data[ : , monomers]
        coms,gr = compute_gr_chainNP_allChains( chains, data_monomers, data_nps )
        coms_all + = [coms]
        grs += [gr]
    com_trajectories = array(coms_all).transpose(1,0,2) # questo èda verificare
    gr_ave = array(grs).mean(0) # questo èda verificare
    return com_trajectories, gr_ave
'''



