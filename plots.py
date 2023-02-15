import matplotlib
matplotlib.use('Agg')
from numpy import loadtxt,diff#,savetxt,shape,argsort
from pylab import *
#from glob import glob
import os
import sys
from scipy import integrate 
sys.path.append('/home/fbrasili/.python')
from init_file import *


def compute_pot(g):
    #g = where( g==0 , 1.e-18 , g )
    return -log(g)

def compute_Fr(r,g):
    #g = where( g==0 , 1.e-18 , g )
    return gradient(g)/gradient(r)/g

def compute_threshold_gr_pot(fol=None):
    if fol is None:
        os.system( "pwd >> path.dat" )
        fol = str( loadtxt( "path.dat" , dtype=str ) )
        os.system( "rm path.dat" )
    gr_path = os.path.join( fol , 'gr_ions_npart_nopbc.dat' )
    r,g = loadtxt( gr_path, skiprows=1 , unpack=True )
    cond = g != 0
    r,g = r[cond],g[cond]
    U = -log(g[1:])
    F = diff(g)/diff(r)/g[1:]
    r,g = r[1:],g[1:]
    savetxt( os.path.join( fol , 'Fr.dat' ) , c_[r,F] )
    cond = r<20
    r,F,U = r[cond],F[cond],U[cond]
    cond = ( r>max(r[F==F.min()]) ) #* ( U+1>0 )
    return r[cond].min()

def compute_threshold_gr_fact(frac=0.2 , fol=None):
    if fol is None:
        os.system( "pwd >> path.dat" )
        fol = str( loadtxt( "path.dat" , dtype=str ) )
        os.system( "rm path.dat" )
    gr_path = os.path.join( fol , 'gr_ions_npart_nopbc.dat' )

    init_fname = os.path.join( fol , "mgel_nanoparts_eq10e6.dat" )
    _,_,box,_,_ = read_init(init_fname)
    V = box[0]**3
    _,N = check_init(init_fname, retPartType=3, ret=True, verb=False)
    
    r,g = loadtxt( gr_path, skiprows=1 , unpack=True )
    cond = ( r > max(r[ g == g.max() ]) ) * ( g <= g.max() * frac )
    return r[cond].min()


def compute_threshold_gr_Coul(fol=None):
    gr_paths = {
            2 : 'gr_neg-cions_npart.dat',
            3 : 'gr_ions_npart_nopbc.dat',
            4 : 'gr_npart_npart_nopbc.dat',
            5 : 'gr_pos-cions_npart.dat'
            }
    if fol is None:
        os.system( "pwd >> path.dat" )
        fol = str( loadtxt( "path.dat" , dtype=str ) )
        os.system( "rm path.dat" )
    init_fname = os.path.join( fol , "mgel_nanoparts_eq10e6.dat" )
    gr_paths = { i: os.path.join( fol , p ) for i,p in gr_paths.items() }
    
    _,_,box,_,_ = read_init(init_fname)
    V = box[0]**3
    Zs,Ns = check_init(init_fname, retPartType="all", ret=True, verb=False)
    Zr = 0
    for key in gr_paths.keys():
        Z,N = Zs[key],Ns[key]
        if key==4: N=N-1
        r,g = loadtxt( gr_paths[key] , skiprows=1 , unpack=True )
        Zr += 4*pi * Z*N/V * integrate.cumtrapz( g * r**2 , r , initial=0 ) 
    Ir = integrate.cumtrapz( Zr/r**2 , r , initial=0 )
    Ur = Zs[4]/r - Ir + Ir[-1]
    
    figure()
    plot(r,Ur)
    savefig("fol"); colse()
    return r[Ur+1<0].max()

def myPlot(observable , system , fol=None , norm=False , ret=False):
    filename = {
            'profile': {
                'nanoparticles': 'nanoparts_profile.dat',
                'microgel': 'monomers_profile.dat',
                'charge': 'charge_profile.dat',
                'positive counterions': 'pos_cions_profile.dat',
                'negative counterions': 'neg_cions_profile.dat'
                },
            'g(R)': {
                'nanoparticle - nanoparticle': 'gr_npart_npart_nopbc.dat', #nr_npart_npart_pbc.dat
                'charged monomer - charged monomer': 'gr_ions_ions_nopbc.dat', #nr_ions_ions_pbc.dat
                'nanoparticle - charged monomer': 'gr_ions_npart_nopbc.dat', #nr_ions_npart_pbc.dat
                'nanoparticle - positive counterion': 'gr_pos-cions_npart.dat',
                'nanoparticle - negative counterion': 'gr_neg-cions_npart.dat',
                'charged monomer - negative counterion': 'gr_ions_neg-cions.dat',
                'charged monomer - positive counterion': 'gr_ions_pos-cions.dat',
                'positive counterion - negative counterion': 'gr_neg-cions_pos-cions.dat'
                },
            'gyration radius': {
                'chull': os.path.join( "chull_rh" , "rh.dat" ),
                'zeno': 'swelling'
                },
            'hydrodynamic radius': {
                'chull': os.path.join( "chull_rh" , "rh.dat" ),
                'zeno': 'swelling'
                }
            }
    args = {
            'profile': [4,(1,3)],
            'g(R)': [1,(0,1)],
            'gyration radius': [0,(0,1)],
            'hydrodynamic radius': [0,(0,2)]
            }
    
    if   observable == 'prof'    : observable = 'profile'
    elif observable == 'gr'      : observable = 'g(R)'
    elif observable == 'Rg'      : observable = 'gyration radius'
    elif observable == 'Rh'      : observable = 'hydrodynamic radius'
    
    if observable == 'profile':
        if system in [ 'nanoparticle' , 'nanoparts' , 'nanopart' , 'nps' , 'np' ]: system,retPartType = 'nanoparticles',4
        elif system in [ 'mgel' , 'monomer' , 'monomers' ]: system,retPartType = 'microgel',None
        elif system == 'pos_cions': system,retPartType = 'positive counterions',5
        elif system == 'neg_cions': system,retPartType = 'negative counterions',2
    elif observable == 'g(R)':
        if system in [ 'nanoparticle' , 'nanoparts' , 'nanopart' , 'nps' , 'np' ]: 
            system = 'nanoparticle - nanoparticle'
            retPartType = 4
        elif system in [ 'ion' , 'ions' , 'monomer' , 'monomers' ]: 
            system = 'charged monomer - charged monomer'
            retPartType = 3
        elif system == 'np-ion': 
            system = 'nanoparticle - charged monomer'
            retPartType = 3
        elif system == 'np-pos_cion': 
            system = 'nanoparticle - positive counterion'
            retPartType = 5
        elif system == 'np-neg_cion': 
            system = 'nanoparticle - negative counterion'
            retPartType = 2
        elif system == 'ion-neg_cion': 
            system = 'charged monomer - negative counterion'
            retPartType = 2
        elif system == 'ion-pos_cion': 
            system = 'charged monomer - positive counterion'
            retPartType = 5
        elif system == 'pos_cion-neg_cion': 
            system = 'positive counterion - negative counterion'
            retPartType = 2
    if system == 'swelling': system = 'chull'

    skiprows,usecols = args[observable]
    if observable == 'profile' and system == 'charge': skiprows = 1
    
    if fol is None:
        os.system( "pwd >> path.dat" )
        fol = str( loadtxt( "path.dat" , dtype=str ) )
        os.system( "rm path.dat" )
    if 'alpha' not in fol:
        if system in ['chull', 'zeno']: fol = os.path.join( fol , 'swelling' )
        else                          : fol = os.path.join( fol , 'alpha_0.00' )
    
    path = os.path.join( fol , filename[observable][system] )
    x,y = loadtxt( path , skiprows=skiprows , usecols=usecols , unpack=True )
    
    if system in ['chull', 'zeno']: label = observable
    else                          : label = system
    
    if type(norm) is float:
        y /= norm
    elif norm:
        init_fname = os.path.join( fol , "mgel_nanoparts_eq10e6.dat")
        if not os.path.isfile( init_fname ): init_fname = os.path.join( fol , "mgel_cions_eq10e6.dat")
        if system == 'microgel':
            _,Nions = check_init(init_fname, retPartType=3, ret=True, verb=False)
            _,Nneut = check_init(init_fname, retPartType=1, ret=True, verb=False)
            Nparts = Nions + Nneut
        else: 
            _,Nparts = check_init(init_fname, retPartType=retPartType, ret=True, verb=False)
        #_,_,box,_,_ = read_init(init_fname)
        y /= Nparts

    if ret:
        return x, y, label
    else:
        plot(x, y, label=label)

'''
non ha senso scrivere una funzione per i due/tre sublplot
def plot_swellingCurve(fol='./', withZeno=False, ret=False):
    Xrg,Yrg,label = myPlot('Rg', 'swelilng' , fol=fol , ret=True)
    plot(x,y)
'''
