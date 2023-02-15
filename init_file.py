#!/usr/bin/python3

from numpy import *
from collections import Counter
import sys
import os
from utils import *

##########################################################################
### READ INFO FROM THE HEADINGS OF THE INIT FILE AND CHECK ITS CONSISTENCY

def read_init(filename):
    box= []
    with open( filename, "r" ) as infile:
        for i,line in enumerate(infile):
            if 'atoms' in line:
                nAtoms = int( line.split()[0] )
            elif 'bonds'in line:
                nBonds = int( line.split()[0] )
            elif array( [ string in line for string in ['xlo','ylo','zlo'] ] ).any():
                line = line.split()
                box += [ float( line[1] ) - float( line[0] ) ]
            elif 'Atoms' in line:
                idxAtoms = i + 2
            elif 'Bonds' in line:
                idxBonds = i+2
                break
    return nAtoms,nBonds,box,idxAtoms,idxBonds

def check_init(filename, retPartType=4, ret=False, verb=True):
    if verb:
        print()
        print("init file: " + filename)
    lut_partType = {
            1: 'neutral microgel monomers',
            2: 'negative microgel counterions',
            3: 'positive microgel monomers',
            4: 'negative nanoparticles',
            5: 'positive nanoparticles counterions'
            }
    
    nAtoms,nBonds,box,skip,_ = read_init(filename)
    
    types,charges = loadtxt(filename, skiprows=skip, max_rows=nAtoms, usecols=(1,5), dtype=float).transpose()
    types.astype(int)

    isTypes = sort( list( set( types ) ) )
    typeCounter = { lut_partType[key]: Counter(types)[key] for key in isTypes }

    #charges = genfromtxt(filename, skip_header=skip, max_rows=nAtoms, usecols=(5), dtype=float).transpose()
    totalCharge = charges.sum()
    particleCharges = { lut_partType[key]: charges[types==key] for key in isTypes }
    
    #for key in particleCharges.keys(): print(key)
    for partType,charges in particleCharges.items(): assert all(charges == charges[0]) , 'not all the particles of type %i have the same charge' %partType
    for key,charges in particleCharges.items() : particleCharges[key] = int(charges[0])
    
    if ret:
        if retPartType is None or retPartType == 'all':
            return { int(key): particleCharges[lut_partType[key]] for key in isTypes }, { int(key): typeCounter[lut_partType[key]] for key in isTypes }
            #return particleCharges.values(), typeCounter.values()
        else:
            if retPartType in lut_partType.keys(): retPartType = lut_partType[retPartType]
            Q = particleCharges[retPartType]
            N = typeCounter[retPartType]
            return Q,N#,nAtoms,box,typeCounter,particleCharges,totalCharge
    else:
        print("")
        print("Number of atoms: %i" %nAtoms)
        print("Number of bonds: %i" %nBonds)
        print("Box size: ", box)
        print("")
        for item in typeCounter.items(): print("number of " + item[0] + ": %i" %item[1])
        print("")
        for item in particleCharges.items(): print("charge of " + item[0] + ": %i" %item[1])
        print("")
        print("total charge of the system: %i" %totalCharge)
        print("")


#############################
### REORDER THE LIST OF ATOMS

def reorder_atoms(filename, outfile='out.init', replace=False):
    atomsRange,velsRange = [],[]
    doVels = False
    with open( filename, "r" ) as infile:
        for i,line in enumerate(infile):
            if 'atoms' in line:
                nAtoms = int( line.split()[0] )
            elif 'Atoms' in line:
                idxAtoms = i + 2
                atomsRange = range(idxAtoms, idxAtoms + nAtoms)
            elif 'Velocities' in line:
                doVels = True
                idxVels = i + 2
                velsRange = range(idxVels, idxVels + nAtoms)
                break
            elif 'Bonds' in line:
                break
    
    linesAtoms = loadtxt( filename, skiprows=idxAtoms, max_rows=nAtoms, dtype=str )
    linesAtoms = linesAtoms[ argsort( linesAtoms[:,0].flatten().astype(int) ) ]
    if doVels:
        linesVels = loadtxt( filename, skiprows=idxVels, max_rows=nAtoms, dtype=str )
        linesVels = linesVels[ argsort( linesVels[:,0].flatten().astype(int) ) ]
    
    of = open( outfile , "w" )
    with open( filename, "r" ) as infile:
        for i,line in enumerate(infile):
            if   i in atomsRange: line = ' '.join( linesAtoms[i-idxAtoms] ) + '\n'
            elif doVels and i in velsRange : line = ' '.join( linesVels[ i-idxVels ] ) + '\n'
            of.write( line )
    
    of.close()
    
    if replace: os.system( 'mv ' + outfile + ' '+ filename )


#########################################################################
### GENERATE A LIST OF CHAINS (WITH CROSSLINKERS OR DANGLING ENDS AS ENDS)

def gen_chains(filename, rem_loops=True, sep_dangling=False, saveChains=True, ret=True):
    
    nAtoms,nBonds,_,idxAtoms,idxBonds = read_init(filename)
    bonds = loadtxt( filename, skiprows=idxBonds, max_rows=nBonds, usecols=(2,3), dtype=int )
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
                freeEnds.update( { key: idx_min for key,value in freeEnds.items() if value == idx_max } )
                freeEnds.update( { key: value - 1 for key,value in freeEnds.items() if value > idx_max } )
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

    n0 = len(chains)
    print( 'total number of chains: %i'%n0 )

    if rem_loops: 
        chains = array( [chain for chain in chains if chain[0] != chain[-1] ] , dtype=list)
        print( 'discarded %i loops'%(n0-len(chains)) )
    
    N = array( [ len(chain) for chain in chains ], dtype=int )
    N,chains = N[argsort(N)][::-1],chains[argsort(N)][::-1]

    if sep_dangling: 
        danglingAtoms = [ key for key,value in atomCounter.items() if value == 1 ]
        danglingChains = array( [ chain for chain in chains if chain[0] in danglingAtoms or chain[-1] in danglingAtoms ] , dtype=list )
        innerChains = array( [ chain for chain in chains if chain not in danglingChains ] , dtype=list )
        if saveChains:
            with open( 'dangling_chains.dat' , "w" ) as of:
                for chain in danglingChains: of.write( '\t'.join(chain) + '\n' )
            with open( 'inner_chains.dat' , "w" ) as of:
                for chain in innerChains: of.write( '\t'.join(chain) + '\n' )
        if ret:
            return innerChains,danglingChains
    else:
        if saveChains:
            with open( 'chains.dat' , "w" ) as of:
                for chain in chains: of.write( '\t'.join( [ str(atom) for atom in chain ] ) + '\n' )
        if ret:
            return N,chains


####################################################
### REMOVE ATOMS OF SPECIFIC TYPE FROM THE INIT FILE

def remove_atoms(filename, atomType, toRemove, outfile='out.init', replace=False):
    '''
    filename: string; name of the init file to modify
    atomType: string or int; type of the atoms to remove
    toRemove: int or "all"; number of atoms to remove
    outfile: string; name of the output init file
    replace: bool; if True the function replaces the init file
                   if False the function creates a new init file, using outfile as name
    '''
    lut_partType = {
            1: 'neutral microgel monomers',
            2: 'negative microgel counterions',
            3: 'positive microgel monomers',
            4: 'negative nanoparticles',
            5: 'positive nanoparticles counterions'
            }
    
    assert filename != outfile, "Choose a name for the output file different from " + filename
    

    for key,val in lut_partType.items():
        if atomType == val: atomType = key
    atomType = int( atomType )
    
    _,nType = check_init(filename, retPartType=atomType, ret=True)
    if toRemove=='all' or int(toRemove)==nType: 
        toRemove = nType
        remType = True
    else:
        toRemove = int(toRemove)
        remType = False
    
    # reorders the Atoms and Velocities lists
    reorder_atoms(filename, outfile='.temp_init.dat', replace=False)

    # open temp output file for writing
    of = open( outfile , "w" )
    
    # open input file for reading
    with open( '.temp_init.dat', "r" ) as infile:
        removedBonds = 0
        atomsRange,massesRange,velsRange,bondsRange,removedAtoms = [],[],[],[],[]
        for i,line in enumerate(infile):
            if 'atoms' in line:
                line = line.split()
                nAtoms = int( line[0] )
                line[0] = str( int( nAtoms - toRemove ) )
                line = ' '.join(line) + '\n'
            if 'atom types' in line:
                line = line.split()
                nTypes = int( line[0] )
                if remType: line[0] = str( nTypes - 1 )
                line = ' '.join(line) + '\n'
            elif 'bonds' in line:
                nBonds = int( line.split()[0] )
            elif 'Masses' in line:
                massesRange = range(i+2, i+2+nTypes)
            elif i in massesRange and remType:
                line = line.split()
                if int( line[0] ) == atomType:
                    continue
                else:
                    line = ' '.join(line) + '\n'
            elif 'Atoms' in line:
                idxLine = 0
                atomsRange = range(i+2, i+2+nAtoms)
                print("Removing %i "%toRemove + lut_partType[atomType] + "...")
            elif i in atomsRange:
                line = line.split()
                if toRemove > 0 and int( line[1] ) == atomType:
                    removedAtoms += [ int( line[0] ) ]
                    toRemove -= 1
                    continue
                else:
                    idxLine += 1
                    line[0] = str( int(idxLine) )
                    line = '\t'.join(line) + '\n'
            elif 'Velocities' in line:
                idxLine = 0
                velsRange = range(i+2, i+2+nAtoms)
            elif i in velsRange:
                line = line.split()
                if int( line[0] ) in removedAtoms:
                    continue
                else:
                    idxLine += 1
                    line[0] = str( int(idxLine) )
                    line = '\t'.join(line) + '\n'
            elif 'Bonds' in line:
                idxLine = 0
                bondsRange = range(i+2, i+2+nBonds)
            elif i in bondsRange:
                line = line.split()
                if int( line[2] ) in removedAtoms or int( line[3] ) in removedAtoms:
                    removedBonds += 1
                    continue
                else:
                    idxLine += 1
                    line[0] = str( int(idxLine) )
                    line = '\t'.join(line) + '\n'
            # write the line in the output file
            of.write( line )
    # close the output file
    of.close()

    os.system( 'rm .temp_init.dat' )
    
    if removedBonds > 0:    ### rivedere come funzione "r+" è necessario riscrivere anche tutte le altre righe? o agisce solo su quella selezionata?
        nBonds -= removedBonds
        with open( outfile , "r+") as of:
            for i,line in enumerate(of):

                if 'bonds' in line:
                    line = line.split()
                    line[0] = str( int(nBonds) )
                    line = ' '.join(line) + '\n'
                    of.write( line )
                    break

    # replace the input init file
    if replace: os.system( 'mv ' + outfile + ' '+ filename )


#######################################################
### ADD A BOND TYPE BETWEEN ATOMS OF SAME SELECTED TYPE

def add_bondType(filename, atomType=1, outfile='out.init', replace=False):
    
    atomType = str( atomType )
    # open temp output file for writing
    of = open( outfile , "w" )
    
    # open input file for reading
    with open( filename, "r" ) as infile:
        idxHead = None
        skipline = False
        atomsRange,bondsRange,neutralMonomers = [],[],[]
        # read info from init files
        for i,line in enumerate(infile):
            if 'atoms' in line:
                nAtoms = int( line.strip().split()[0] )
            elif 'bonds' in line:
                nBonds = int( line.strip().split()[0] )
            # add 1 bond type
            elif 'bond types' in line:
                line = line.split()
                newBondType = str( int( line[0] ) + 1 )
                line[0] = newBondType
                line = ' '.join(line) + '\n'
            # skip the 'Bond Coeffs' lines 
            if 'Bond Coeffs' in line:
                skipline = True
            # identify the 'Atoms' lines
            elif 'Atoms' in line:
                line = 'Atoms\n'
                skipline = False
                atomsRange = range(i+2, i+2+nAtoms)
            # identify neutral (selected) atoms
            elif i in atomsRange:
                words = line.strip().split()
                if words[1] == atomType: neutralMonomers += [ words[0] ]
            # identify 'Bonds' lines
            elif 'Bonds' in line:
                bondsRange = range(i+2, i+2+nBonds)
            # change the type of bonds between neutral (selected) atoms
            elif i in bondsRange:
                line = line.strip().split()
                if line[2] in neutralMonomers and line[3] in neutralMonomers: line[1] = newBondType
                line = '\t'.join(line) + '\n'

            # write the line in the output file
            if skipline: continue
            else: of.write( line )
    
    # close the output file
    of.close()

    # replace the input init file
    if replace: os.system( 'mv ' + outfile + ' '+ filename )



##########################################################
### REMOVE A BOND TYPE BETWEEN ATOMS OF SAME SELECTED TYPE

def rem_bondType(filename, bondTypeToRemove=2, newBondType=1, outfile='out.init', replace=False):
    newBondType = str( newBondType )
    # open temp output file for writing
    of = open( outfile , "w" )

    # open input file for reading
    with open( filename, "r" ) as infile:
        idxHead = None
        skipline = False
        atomsRange,bondsRange,neutralMonomers = [],[],[]
        # read info from init files
        for i,line in enumerate(infile):
            if 'bonds' in line:
                nBonds = int( line.strip().split()[0] )
            # remove 1 bond type
            elif 'bond types' in line:
                line = line.split()
                #newBondType = str( int( line[0] ) - 1 )
                line[0] = str( int( line[0] ) - 1 )
                line = ' '.join(line) + '\n'
            # skip the 'Bond Coeffs' lines 
            elif 'Bond Coeffs' in line:
                skipline = True
            # identify the 'Atoms' lines
            elif 'Atoms' in line:
                skipline = False
            # identify 'Bonds' lines
            elif 'Bonds' in line:
                bondsRange = range(i+2, i+2+nBonds)
            # change the type of bonds between neutral (selected) atoms
            elif i in bondsRange:
                line = line.strip().split()
                if int( line[1] ) == bondTypeToRemove: line[1] = newBondType
                line = '\t'.join(line) + '\n'
            
            # write the line in the output file
            if skipline: continue
            else: of.write( line )

    # close the output file
    of.close()

    # replace the input init file
    if replace: os.system( 'mv ' + outfile + ' '+ filename )


####################################################################################################
### CHANGE THE SIZE OF THE BOX (RECENTER THE MICROGEL IN THE CENTER OF THE BOX AND REWRAP ALL ATOMS)
def change_box(filename, newBox, outfile=None, replace=False):
    assert newBox%2 == 0, "use an even integer for the box size"
    newBox = [ -newBox/2 , newBox/2 ]

    # open input file for reading
    with open( filename, "r" ) as infile:
        atomsRange=[]
        for i,line in enumerate(infile):
            if 'atoms' in line:
                nAtoms = int( line.split()[0] )
            elif 'Atoms' in line:
                skiprows = i+2
                atomsRange = range(i+2, i+2+nAtoms)
                break
    datas = loadtxt(filename, skiprows=skiprows, max_rows=nAtoms)
    datas = recenter(datas)
    datas = rewrap( datas, tile(newBox,(3,1)) )
    #nvals = range( shape(datas)[1] )
    
    if outfile is None: outfile = 'out.init'
    if replace: outfile = '.temp_init'
    # open temp output file for writing
    of = open( outfile , "w" )
    
    # open input file for reading
    with open( filename, "r" ) as infile:
        for i,line in enumerate(infile):
            if array( [ string in line for string in ['xlo','ylo','zlo'] ] ).any():
                lo,hi = newBox
                line = line.split()
                line[0],line[1] = str(lo),str(hi)
                line = ' '.join(line) + '\n'
            elif i in atomsRange:
                data = datas[i-skiprows]
                line = []
                for j,val in enumerate(data):
                    if j in [0,1,6]:
                        line += [ str( int( val ) ) ]
                    else:
                        line += [ str( val ) ]
                del data
                line = '\t'.join(line) + '\n'
            of.write(line)
    del datas
    # close the output file
    of.close()
    # replace the input init file
    if replace: os.system( 'mv ' + outfile + ' '+ filename )


##################################
### OTHER CHANGES IN THE INIT FILE

def modify_init(filename, partType=None, newCharge=None, newType=None, convert=False, outfile=None, replace=False):
    if outfile is None: outfile = 'out.init'
    if replace: outfile = '.temp_init.dat'
    #if not replace and outfile is None: outfile = 'out.init'
    #elif replace and outfile is None: outfile = '.temp_init'
    if newCharge is not None or newType is not None:
        assert partType is not None, 'Specify the type of particle to modify' 
        partType = int( partType )

    # open temp output file for writing
    of = open( outfile , "w" )
                                                                                                                                        
    # open input file for reading
    with open( filename, "r" ) as infile:
        atomsRange=[]
        for i,line in enumerate(infile):
            if 'atoms' in line:
                nAtoms = int( line.split()[0] )
            # change the number of atom types
            elif 'atom types' in line and newType is not None:
                line = line.split()
                line[0] = str( int(line[0]) - 1 )
                line = ' '.join(line) + '\n'
            # identify the 'Atoms' lines
            elif 'Atoms' in line:
                atomsRange = range(i+2, i+2+nAtoms)
            # change particle properties
            elif i in atomsRange:
                line = line.split()
                if int( line[1] ) == partType:
                    # change charge
                    if newCharge is not None:
                        if newCharge == 'invert': line[5] = '-' + line[5]
                        else                    : line[5] = newCharge
                    # change type
                    if newType is not None:
                        line[1] = str( newType )
                line = '\t'.join(line) + '\n'
            
            # write the line in the output file
            of.write( line )

    # close the output file
    of.close()
    # replace the input init file
    if replace: os.system( 'mv ' + outfile + ' '+ filename )
