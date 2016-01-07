
import os as os
import sys as sys
import datetime as dt
import numpy as np
import prody as pd

######################################################################################################################################################
# Container for performing an Alanine Scan with AESOP
#   pdb     -   PDB file for performing Alascan. Must contain all chain selections with standard aminoacid nomenclature
#   selstr  -   List of selection strings for performing mutations
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant   
#   sdie    -   Solvent dielectric constant
######################################################################################################################################################
class Alascan:
    def __init__(self, pdb, selstr=['Protein'], ion=0.150, pdie=20.0, sdie=78.54):
        self.pdb = pdb
        self.selstr = selstr
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
    # pdb is the pdb file
    # residue is the residue number
    # mutid is the prefix for the written mutated structure
    # resid is the residue to mutate to

    from modeller import environ, model, alignment, selection

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
def execPDB2PQR(path_pdb2pqr_exe, pdbfile, optargs='--ff charmm --chain --apbs-input', outfile=None):
    import os as os
    if outfile is None:
        outfile = os.path.splitext(pdbfile)[0]+'.pqr'
    os.system('{0} {1} {2} {3}'.format(path_pdb2pqr_exe, optargs, pdbfile, outfile))
    return outfile