import os
import sys
from numpy import *
from utils import *

###########################################
# CHECK CONFIGURATIONS FILE
# - read number of atoms
# - read first and last timestep
# - read the interval between timesteps
# - check that timesteps are consecutive without repetitions or jumps

def readConfigurations(filename='configurations.dat', retVals=False, check=True, verb=True):
    if check: retVals=False
    if verb: print( 'analyzing configurations file ' + filename + ' ...' )
    timesteps = []
    repeatedTimesteps = []
    with open(filename, 'r') as f:
        nAtoms,linesXtimestep = 1,-2
        for i,line in enumerate(f):
            if i == 1:
                firstTimestep = int( line.strip() )
                timesteps += [ firstTimestep ]
            elif i == 3:
                nAtoms = int( line.strip() )
                if verb:
                    print( 'number of atoms: %i'%nAtoms )
                    print( 'first timestep %i'%firstTimestep )
            elif linesXtimestep < 0 and line.split()[0] == '1':
                linesXtimestep = i + nAtoms
            elif i == linesXtimestep + 1:
                secondTimestep = int( line.strip() )
                assert secondTimestep >= firstTimestep, 'the second step preceds the first!!'
                timesteps += [ secondTimestep ]
                interval = secondTimestep - firstTimestep
                if verb: print( 'time interval between configurations: %i'%interval )
            elif i%(linesXtimestep) == 1:
                timestep = int( line.strip() )
                assert timestep > firstTimestep, 'there are steps preceding the first!!'
                if timestep in timesteps: repeatedTimesteps += [ timestep ]
                elif timestep > timesteps[-1] + interval: print( 'missing configurations between timestep %i'%timesteps[-1] + ' and timestep %i'%timestep )
                elif timestep < timesteps[-1] + interval: print( 'timstep %i'%timestep + ' follows timestep %i'%timesteps[-1] )
                timesteps += [timestep]
    
    lastTimestep = timesteps[-1]
    if verb: 
        print( 'last timestep: %i'%lastTimestep )
        if len(repeatedTimesteps) > 0:
            print( 'the following timesteps are repeated:' )
            for repTs in repeatedTimesteps: print( repTs )
    if verb: print( 'check completed.' ) ; print()
    if check:
        return check
    elif retVals:
        return nAtoms, firstTimestep, interval, linesXtimestep, lastTimestep, array( repeatedTimesteps )

################################################
# EXTRACT CONSECUTIVE CONFIGURATIONS FROM A FILE
def extractConfigurations(path, n, first=0, out_fname="selected_configurations.dat", verb=True ):
    assert os.path.isfile(path) , "the file " + path + " does not exist"
    assert first >= 0
    assert n > 0 , "select at least one configuration"
    if verb: print( "extracting the selected configurations from file " + path )
    _,f,interval,lines,l,_ = readConfigurations( filename=path, retVals=True, check=False, verb=False )
    assert n <= 1 + ( l - f ) / interval - first , "not enough timesteps in the file " + path
    
    os.system( "head -n %i "%(lines*(n+first)) + path + " | tail -n %i > "%(lines*n) + out_fname )
    if verb: print( "selected configurations saved to file " + out_fname )
    

################################
# MERGE TWO CONFIGURATIONS FILES
def mergeConfigurations(path1, path2, replace=True, out_fname="mergedConfs.dat"):
    assert all([ os.path.isfile(p) for p in [path1,path2] ]) , "one of the two files to merge does not exist"
    n1,first1,int1,l1,last1,_ = readConfigurations( filename=path1, retVals=True, check=False, verb=False )
    n2,first2,int2,l2,last2,_ = readConfigurations( filename=path2, retVals=True, check=False, verb=False )
    assert n1 == n2 , "the number of atoms in the two files is different, cannot proceed merging"
    if int1 != int2: print( "the timestep interval of the two files is different!!" )
    nLines1 = l1 * int( (first2 - first1)/int1 )
    os.system( "head -n %i "%nLines1 + path1 + " > .temp_data " )
    os.system( "cat " + path2 + " >> .temp_data " )
    if replace:
        os.system( "mv .temp_data " + path1 )
        os.system( "rm " + path2 )
    else:
        os.system( "mv .temp_data " + out_fname )

#######################################################
# LOAD DATA FROM SINGLE CONFIGURATION, RECENTER, REWRAP

def loadConfiguration(conf_file, skiprows=None, max_rows=None, usecols=None):
    data = loadtxt(conf_file, skiprows=skiprows, max_rows=max_rows )
    data = recenter(data)
    box = loadtxt(conf_file, skiprows=5, max_rows=3)
    data = rewrap(data, box)
    return data[:, array(usecols)]
'''
def loop_on_configs(filename='configurations.dat', function):
    print( 'configurations file: ' + filename)
    #timesteps = []
    #repeatedTimesteps = []
    with open(filename, 'r') as f:
        nAtoms,linesXtimestep = 1,-2
        for i,line in enumerate(f):
            if i == 1:
                firstTimestep = int( line.strip() )
            elif i == 3:
                nAtoms = int( line.strip() )
                if printInfo:
            elif linesXtimestep < 0 and line.split()[0] == '1':
                linesXtimestep = i + nAtoms
            elif i == linesXtimestep + 1:
                interval = secondTimestep - firstTimestep
           elif i%(linesXtimestep) == 1:

'''
###########################################
# REMOVE TIMESTEPS IN CONFIGURATIONS FILE

def removeTimesteps(whatRemove = 'head', filename='configurations.dat', firstTimstep2keep = None, outfile='new_confs.dat', replace=False):
    nAtoms, first, interval, lines4timestep, last, repeated = checkConfigurations(filename=filename, ret=True, printInfo=False)
    if   whatRemove == 'doubled':
        timesteps2remove = repeated
    elif whatRemove == 'head':
        assert firstTimstep2keep is not None, 'insert the number of the first timesteps to keep'
        assert firstTimstep2keep > first, 'the first timestep is already %i'%first
        timesteps2remove = arange(first, firstTimstep2keep, interval)
    #lines2remove = concatenate( [ ( n * lines4timestep ) + arange( lines4timestep ) for n in (timesteps2remove - first)/interval ] ).astype(int)
    
    removeFrom,removeTo = timesteps2remove.min(),timesteps2remove.max()
    print( 'removing the the configurations between timestep %i'%removeFrom + ' and timestep %i ...'%removeTo )
    
    headLines2keep = int( lines4timestep * ( removeFrom - first    )/interval )
    tailLines2keep = int( lines4timestep * ( last       - removeTo )/interval )
    
    if os.path.isfile(outfile): os.system( 'rm ' + outfile )
    os.system( 'head -n %i '%headLines2keep + filename + ' >> ' + outfile )
    os.system( 'tail -n %i '%tailLines2keep + filename + ' >> ' + outfile )

    if replace: os.system( 'mv ' + outfile + ' ' + filename )
    else:       print( 'output file: ' + outfile )
    print('Done.')
