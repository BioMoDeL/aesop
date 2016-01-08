
import os as os
import sys as sys
import datetime as dt
import numpy as np
import prody as pd
from modeller import environ, model, alignment, selection

######################################################################################################################################################
# Container for performing an Alanine Scan with AESOP
#   pdb     -   PDB file for performing Alascan. Must contain all chain selections with standard aminoacid nomenclature
#   selstr  -   List of selection strings for performing mutations
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant   
#   sdie    -   Solvent dielectric constant
######################################################################################################################################################
class Alascan:
    """Summary
    
    Attributes
    ----------
    apbs : TYPE
        Description
    dirs : int
        Description
    ion : TYPE
        Description
    pdb : TYPE
        Description
    pdb2pqr : TYPE
        Description
    pdie : TYPE
        Description
    prefix : TYPE
        Description
    sdie : TYPE
        Description
    selstr : TYPE
        Description
    """
    def __init__(self, pdb, pdb2pqr_exe, apbs_exe, selstr=['protein'], region=None, ion=0.150, pdie=20.0, sdie=78.54):
        self.pdb = pdb
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        self.selstr = selstr
        if region is None:
            region = selstr
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        # Insert code to instantiate dirs and prefix
        self.dirs = 0
        self.prefix = '%4d%02d%02d'%(dt.date.today().year, dt.date.today().month, dt.date.today().day)

    def getPDB(self):
        return self.pdb
    def getSel(self):
        return self.selstr
    def getDirs(self):
        return self.dirs
    def getPrefix(self):
        return self.prefix
    def getEnergies():
        return 0
    def getMutID():
        return 0

    def genPQR():
        0
    def genMut():
        0
    def batchAPBS():
        0

    def run():
        0

######################################################################################################################################################
# Function to mutate a single residue in a PDB structure
######################################################################################################################################################
def mutatePDB(pdb, mutid, resnum, chain=None, resid='ALA'):
    """Summary
    
    Parameters
    ----------
    pdb : TYPE
        Description
    mutid : TYPE
        Description
    resnum : TYPE
        Description
    chain : TYPE, optional
        Description
    resid : str, optional
        Description
    
    Returns
    -------
    name : TYPE
        Description
    """
    # pdb is the pdb file
    # resnum is the residue number to be mutated
    # chain (optional) can specify what chain the residue to be mutated is located on
    # mutid is the prefix for the written mutated structure to be written
    # resid is the residue to mutate to

    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    aln = alignment(env)
    mdl = model(env, file=pdb)
    aln.append_model(mdl, atom_files=pdb, align_codes='parent')

    if chain is None:
        sel = selection(mdl.residue_range(resnum-1, resnum-1))
    else:
        sel = selection(mdl.residue_range(str(resnum)+':'+chain, str(resnum)+':'+chain))

    sel.mutate(residue_type=resid)

    aln.append_model(mdl, align_codes='mutant')
    mdl.clear_topology()
    mdl.generate_topology(aln['mutant'])
    mdl.transfer_xyz(aln)

    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.write(file=mutid+'.pdb')

######################################################################################################################################################
# Function to run PDB2PQR.exe - should work on any supported OS
######################################################################################################################################################
def execPDB2PQR(path_pdb2pqr_exe, pdbfile, optargs='--ff charmm --chain', outfile=None):
    """Summary
    
    Parameters
    ----------
    path_pdb2pqr_exe : TYPE
        Description
    pdbfile : TYPE
        Description
    optargs : str, optional
        Description
    outfile : TYPE, optional
        Description
    
    Returns
    -------
    name : TYPE
        Description
    """
    if outfile is None:
        outfile = os.path.splitext(pdbfile)[0]+'.pqr'
    os.system('{0} {1} {2} {3}'.format(path_pdb2pqr_exe, optargs, pdbfile, outfile))
    return outfile

######################################################################################################################################################
# Function to run APBS.exe - should work on any supported OS
######################################################################################################################################################
def execAPBS(path_apbs_exe, pqr_chain, pqr_complex, prefix=None, grid=1.0, ion=0.150, pdie=20.0, sdie=78.54):
    """Summary
    
    Parameters
    ----------
    path_apbs_exe : TYPE
        Description
    pqr_chain : TYPE
        Description
    pqr_complex : TYPE
        Description
    prefix : TYPE, optional
        Description
    grid : float, optional
        Description
    ion : float, optional
        Description
    pdie : float, optional
        Description
    sdie : float, optional
        Description
    
    Returns
    -------
    name : TYPE
        Description
    """
    # path_apbs_exe -   full path to apbs executable ('C:\\APBS\\apbs.exe')
    # pqr_chain     -   path to file with single chain pqr (mutant or parent)
    # pqr_complex   -   path to file with complex pqr, contains mutant or parent chain that is used in pqr_chain
    # prefix        -   string to pre-pend to output files (log file, dx file, energy file)
    # grid          -   grid spacing using in APBS calculation
    # ion           -   ionic strength for calculation
    # pdie          -   protein dielectric constant
    # sdie          -   solvent dielectric constant

    if prefix is None:
        prefix = os.path.splitext(pqr_chain)[0]

    cfac = 1.5 # hard-coded scaling factor for mesh dimension, for now

    pqr = pd.parsePQR(pqr_complex)
    coords = pqr.getCoords()
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]

    # Determine mesh dimensions according to Ron's AESOP protocol in the R source file
    fg = np.array((np.ceil(np.max(x)-np.min(x)), np.ceil(np.max(y)-np.min(y)), np.ceil(np.max(z)-np.min(z))))
    fg = np.ceil((fg + 5) * cfac)
    dime_list = (32*np.linspace(1,100,100))+1   # list of possible dime values
    dime_ind = np.ceil(fg/(32*grid))    # index of dime to use from list

    glen = fg
    dime = np.array((dime_list[dime_ind[0]], dime_list[dime_ind[1]], dime_list[dime_ind[2]]))

    # Format APBS input file
    cmd_read = ['read\n',
                '   mol pqr %s\n'%(pqr_chain),
                '   mol pqr %s\n'%(pqr_complex),
                'end\n']
    cmd_solv = ['elec name solv\n',
                '   mg-manual\n',
                '   dime %d %d %d\n'%(dime[0], dime[1], dime[2]),
                '   glen %d %d %d\n'%(glen[0], glen[1], glen[2]),
                '   gcent mol 2\n',
                '   mol 1\n',
                '   lpbe\n',
                '   bcfl sdh\n',
                '   srfm smol\n',
                '   chgm spl2\n',
                '   ion 1 %.2f 2.0\n'%(ion),
                '   ion -1 %.2f 2.0\n'%(ion),
                '   pdie %.2f\n'%(pdie),
                '   sdie %.2f\n'%(sdie),
                '   sdens 10.0\n',
                '   srad 0.0\n',
                '   swin 0.3\n',
                '   temp 298.15\n',
                '   calcenergy total\n',
                '   write pot dx %s\n'%(prefix),
                'end\n']
    cmd_ref = ['elec name ref\n',
                '   mg-manual\n',
                '   dime %d %d %d\n'%(dime[0], dime[1], dime[2]),
                '   glen %d %d %d\n'%(glen[0], glen[1], glen[2]),
                '   gcent mol 2\n',
                '   mol 1\n',
                '   lpbe\n',
                '   bcfl sdh\n',
                '   srfm smol\n',
                '   chgm spl2\n',
                '   ion 1 %.2f 2.0\n'%(ion),
                '   ion -1 %.2f 2.0\n'%(ion),
                '   pdie %.2f\n'%(pdie),
                '   sdie %.2f\n'%(pdie),
                '   sdens 10.0\n',
                '   srad 0.0\n',
                '   swin 0.3\n',
                '   temp 298.15\n',
                '   calcenergy total\n',
                'end\n']
    cmd_write = ['print elecEnergy solv end\n',
                'print elecEnergy ref end\n',
                'quit\n']
    apbs_in = cmd_read + cmd_solv + cmd_ref + cmd_write

    # Write APBS input file
    file_apbs_in = prefix+'.in'
    file_apbs_log = prefix+'.log'
    with open(file_apbs_in, 'w') as f:
        for line in apbs_in:
            f.write(line)

    # Execute APBS
    os.system('{0} {1} {2}'.format(path_apbs_exe, '--output-file=%s.log --output-format=flat'%(prefix), file_apbs_in))
    # os.system('{0} {1}'.format(path_apbs_exe, file_apbs_in))

    return prefix



