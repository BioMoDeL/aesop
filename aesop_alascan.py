
import os as os
import sys as sys
import subprocess as sp
import datetime as dt
import re as re
import numpy as np
import prody as pd
import matplotlib.pyplot as plt
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
    def __init__(self, pdb, pdb2pqr_exe, apbs_exe, coulomb_exe=None, selstr=['protein'], jobname=None, region=None, grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse'):
        self.pdb = pdb
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        if coulomb_exe is None:
            self.coulomb = os.path.split(apbs_exe)[0]
        else:
            self.coulomb = coulomb_exe
        self.selstr = selstr
        if region is None:
            self.region = selstr
        else:
            self.region = region
        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        if jobname is None:
            self.jobname = '%4d%02d%02d_%02d%02d%02d'%(dt.date.today().year, dt.date.today().month, dt.date.today().day, dt.datetime.now().hour, dt.datetime.now().minute, dt.datetime.now().second)
        else:
            self.jobname = jobname
        self.jobdir = jobname
        if not os.path.exists(os.path.join(self.jobdir)):
            os.makedirs(os.path.join(self.jobdir))
        self.E_ref = np.zeros(0)
        self.E_solv = np.zeros(0)
        self.mutid = []
        self.ff = ff

    def getPDB(self):
        return self.pdb
    def getSel(self):
        return self.selstr
    def getDirs(self):
        return self.dirs
    def getPrefix(self):
        return self.prefix
    # def getEnergies(self):
    #     return 0
    def getMutids(self):
        l = self.mutid
        return [item for sublist in l for item in sublist]

    def genDirs(self):
        # Create necessary directories for PDB files
        pdb_complex_dir = 'complex_pdb'
        if not os.path.exists(os.path.join(self.jobdir, pdb_complex_dir)):
            os.makedirs(os.path.join(self.jobdir, pdb_complex_dir))
        self.pdb_complex_dir = pdb_complex_dir

        # Create necessary directories for PQR files
        pqr_complex_dir = 'complex_pqr'
        if not os.path.exists(os.path.join(self.jobdir, pqr_complex_dir)):
            os.makedirs(os.path.join(self.jobdir, pqr_complex_dir))
        pqr_sel_dir=[]
        for i in xrange(0,len(self.selstr)):
            pqr_sel_dir.append('seg%d_pqr'%(i+1))
            if not os.path.exists(os.path.join(self.jobdir, 'seg%d_pqr'%(i+1))):
                os.makedirs(os.path.join(self.jobdir, 'seg%d_pqr'%(i+1)))
        self.pqr_complex_dir = pqr_complex_dir
        self.pqr_sel_dir = pqr_sel_dir

        # Create necessary directories for APBS files
        logs_apbs_dir = 'apbs_logs'
        if not os.path.exists(os.path.join(self.jobdir, logs_apbs_dir)):
            os.makedirs(os.path.join(self.jobdir, logs_apbs_dir))
        self.logs_apbs_dir = 'apbs_logs'


    def genMutid(self):
        selstr = self.selstr
        region = self.region

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        list_mutids = [[] for x in xrange(len(selstr)+1)]
        list_chids = [[] for x in xrange(len(selstr)+1)]
        list_resnums = [[] for x in xrange(len(selstr)+1)]
        list_resnames = [[] for x in xrange(len(selstr)+1)]
        
        list_mutids[0] = [parent_file_prefix]
        list_chids[0] = ['']
        list_resnums[0] = [np.zeros(0)]
        list_resnames[0] = ['']

        index = np.linspace(1,len(selstr),len(selstr)).astype(int)
        for i, sel, reg in zip(index, selstr, region):
            combined_selection = parent_pdb.select(''.join(['(', ') and ('.join((sel, reg, 'charged', 'calpha')), ')']))
            list_chids[i] = combined_selection.getChids().tolist()
            list_resnums[i] = combined_selection.getResnums().tolist()
            list_resnames[i] = combined_selection.getResnames().tolist()
            code = ['seg%d'%(i)+'_'+AA_dict[res_id]+res_no+'A' for ch_id, res_no, res_id in zip(list_chids[i], map(str, list_resnums[i]), list_resnames[i])]
            list_mutids[i] = code

        dim_sel = len(list_mutids)
        dim_mut = np.sum([len(x) for x in list_mutids])
        mask_by_sel = np.zeros((dim_mut, dim_sel)).astype(bool)
        counter = 0
        for i in xrange(dim_sel):
            for j in xrange(counter, counter+len(list_mutids[i]), 1):
                mask_by_sel[j, i] = True
            counter += len(list_mutids[i])

        self.mutid = list_mutids
        self.list_chids = list_chids
        self.list_resnums = list_resnums
        self.list_resnames = list_resnames
        self.mask_by_sel = mask_by_sel

    def genParent(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        # list_mutids = [item for sublist in self.mutid for item in sublist]
        # list_chids = [item for sublist in self.list_chids for item in sublist]
        # list_resnums = [item for sublist in self.list_resnums for item in sublist]

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix+'.pdb')
        system = parent_pdb.select('(('+') or ('.join(selstr)+'))')
        print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)

    def genPDB(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        list_mutids = [item for sublist in self.mutid for item in sublist]
        list_chids = [item for sublist in self.list_chids for item in sublist]
        list_resnums = [item for sublist in self.list_resnums for item in sublist]

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix+'.pdb')
        system = parent_pdb.select('(('+') or ('.join(selstr)+'))')
        print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)

        for mutid, chain, resnum in zip(list_mutids[1:], list_chids[1:], list_resnums[1:]):
            outpath = os.path.join(jobdir, pdb_complex_dir, mutid)
            print '\nGenerating PDB for mutant: %s'%(mutid)
            mutatePDB(pdb=infile, mutid=outpath, resnum=resnum, chain=chain, resid='ALA')

    def genTruncatedPQR(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_pdb2pqr = self.pdb2pqr
        ff = self.ff

        list_mutids = [item for sublist in self.mutid for item in sublist]
        list_chids = [item for sublist in self.list_chids for item in sublist]
        list_resnums = [item for sublist in self.list_resnums for item in sublist]

        # Calculate PQR for parent
        infile = os.path.join(jobdir, pdb_complex_dir, list_mutids[0]+'.pdb')
        outfile = os.path.join(jobdir, pqr_complex_dir, list_mutids[0]+'.pqr')
        print '\nGenerating PQR for mutant: %s'%(list_mutids[0])
        execPDB2PQR(path_pdb2pqr, infile, outfile=outfile, ff=ff)
        complex_pqr = pd.parsePQR(outfile)
        for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, list_mutids[0]+'.pqr')
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)

        for mutid, chain, resnum in zip(list_mutids[1:], list_chids[1:], list_resnums[1:]):
            outpath = os.path.join(jobdir, pqr_complex_dir, mutid)
            print '\nGenerating PQR for mutant: %s'%(mutid)
            # print 'mutid %s, chain %s, resnum %d'%(mutid, chain, resnum)
            # print outpath+'.pqr'
            mutatePQR(outfile, mutid=outpath, resnum=resnum, chain=chain)
            complex_pqr = pd.parsePQR(outpath+'.pqr')
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid+'.pqr')
                # print selfile
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)

    def genPQR(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_pdb2pqr = self.pdb2pqr
        ff = self.ff

        list_mutids = [item for sublist in self.mutid for item in sublist]

        for mutid in list_mutids:
            print '\nGenerating PQR for mutant: %s'%(mutid)
            infile = os.path.join(jobdir, pdb_complex_dir, mutid+'.pdb')
            outfile = os.path.join(jobdir, pqr_complex_dir, mutid+'.pqr')
            execPDB2PQR(path_pdb2pqr, infile, outfile=outfile, ff=ff)
            complex_pqr = pd.parsePQR(outfile)
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid+'.pqr')
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)

    def calcAPBS(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr)+1

        Gsolv = np.zeros((dim_mutid, dim_sel))
        Gref = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\nCalculating solvation and reference energies for mutant: %s'%(mutid)
            complex_pqr = os.path.join(jobdir, pqr_complex_dir, mutid+'.pqr')
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir]+pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid+'.pqr')
                path_prefix_log = os.path.join(jobdir, logs_apbs_dir, mutid)
                energies = execAPBS(path_apbs, subunit_pqr, complex_pqr, prefix=path_prefix_log, grid=1.0, ion=0.150, pdie=20.0, sdie=78.54)
                print energies[0][0]
                print energies[0][1]
                Gsolv[i,j] = energies[0][0]
                Gref[i,j] = energies[0][1]

        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcCoulomb(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_coulomb = self.coulomb
        pdie = self.pdie

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr)+1

        Gcoul = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\nCalculating coulombic energies for mutant: %s'%(mutid)
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir]+pqr_sel_dir):
            	subunit_pqr = os.path.join(jobdir, seldir, mutid+'.pqr')
            	energies = execCoulomb(path_coulomb, subunit_pqr)
            	Gcoul[i,j] = energies/pdie

    	self.Gcoul = Gcoul

    def ddGbind_rel(self):
        Gsolv = self.Gsolv
        Gref = self.Gref
        Gcoul = self.Gcoul

        dGsolv = Gsolv-Gref
        dGsolu = Gsolv[:,0]-Gsolv[:,1:].sum(axis=1)
        dGcoul = Gcoul[:,0]-Gcoul[:,1:].sum(axis=1)
        ddGsolv = dGsolv[:,0]-dGsolv[:,1:].sum(axis=1)

        dGbind = ddGsolv + dGcoul
        dGbind_rel = dGbind - dGbind[0]
        return dGbind_rel

    def run(self):
        self.genDirs()
        self.genMutid()
        self.genPDB()
        self.genPQR()
        self.calcAPBS()
        self.calcCoulomb()

    def run_truncated(self):
        self.genDirs()
        self.genMutid()
        self.genParent()
        self.genTruncatedPQR()
        self.calcAPBS()
        self.calcCoulomb()

    def summary(self):
        plotResults(self)
 

######################################################################################################################################################
# Function to run commands, recording output
######################################################################################################################################################
def runProcess(command):
    proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    # proc = sp.Popen(command, stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    # print "program output:", out
    return (out, err)

######################################################################################################################################################
# Function to mutate a single residue in a PDB structure, mutates with side-chain truncation
######################################################################################################################################################
def mutatePQR(pqrfile, mutid, resnum, chain=None): # Only use this function with PARSE for now ...
    parent = pd.parsePQR(pqrfile)
    if chain is None:
        residue = parent.select('resnum %d'%(int(resnum)))
        preceed = parent.select('resnum < %d'%(int(resnum)))
        follow = parent.select('resnum > %d'%(int(resnum)))
    elif chain is not None:
        residue = parent.select('chain %s and resnum %d'%(str(chain), int(resnum)))
        preceed = parent.select('chain %s and resnum < %d'%(str(chain), int(resnum)))
        follow = parent.select('chain %s and resnum > %d'%(str(chain), int(resnum)))
        otherchains = parent.select('not chain %s'%(str(chain)))
    bb = residue.select('not sidechain')
    sc = residue.sidechain
    cg = residue.select('name CG')
    cb = residue.select('name CB')

    # Set charge and radii of side chain to 0, change CG atom to HB1
    residue.setResnames('ALA')
    sc.setCharges(0)
    sc.setRadii(0)
    cg.setNames('HB1')

    # Shorten the HB1-CB bond
    pos_hb1 = (0.7105*(cg.getCoords() - cb.getCoords())) + cb.getCoords()
    cg.setCoords(pos_hb1)

    # Compile mutated pdb
    ala_atoms = ['N','H','H2','H3','CA','HA','CB','HB1','HB2','HB3','C','O','OXT']
    if chain is None:
        if preceed is None:
            mutant = residue.select('name '+' '.join(ala_atoms))+follow
        if follow is None:
            mutant = preceed+residue.select('name '+' '.join(ala_atoms))
        if (preceed is None) and (follow is None):
            mutant = residue.select('name '+' '.join(ala_atoms))
        if (preceed is not None) and (follow is not None):
            mutant = preceed+residue.select('name '+' '.join(ala_atoms))+follow
    else:
        if otherchains is None:
            if preceed is None:
                mutant = residue.select('name '+' '.join(ala_atoms))+follow
            if follow is None:
                mutant = preceed+residue.select('name '+' '.join(ala_atoms))
            if (preceed is None) and (follow is None):
                mutant = residue.select('name '+' '.join(ala_atoms))
            if (preceed is not None) and (follow is not None):
                mutant = preceed+residue.select('name '+' '.join(ala_atoms))+follow
        if otherchains is not None:
            if preceed is None:
                mutant = residue.select('name '+' '.join(ala_atoms))+follow+otherchains
            if follow is None:
                mutant = preceed+residue.select('name '+' '.join(ala_atoms))+otherchains
            if (preceed is None) and (follow is None):
                mutant = residue.select('name '+' '.join(ala_atoms))+otherchains
            if (preceed is not None) and (follow is not None):
                mutant = preceed+residue.select('name '+' '.join(ala_atoms))+follow+otherchains

    # Write mutant pqr
    pd.writePQR(mutid+'.pqr', mutant)

######################################################################################################################################################
# Function to mutate a single residue in a PDB structure, mutates with modeller by building internal coordinates of residue
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

    if((chain is None) or (chain.isspace())):
        sel = selection(mdl.residue_range(int(resnum)-1, int(resnum)-1))
    else:
        sel = selection(mdl.residue_range(str(resnum)+':'+chain, str(resnum)+':'+chain))

    sel.mutate(residue_type=resid)

    aln.append_model(mdl, align_codes='mutant')
    mdl.clear_topology()
    mdl.generate_topology(aln['mutant'])
    mdl.transfer_xyz(aln)

    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.write(file=mutid+'.pdb')

    h = model(env, file=mutid+'.pdb') # Without this section, chainids and resnums from parent won't be retained!
    m = model(env, file=pdb)
    aln = alignment(env)
    aln.append_model(m, atom_files=pdb, align_codes='parent')
    aln.append_model(h, atom_files=mutid+'.pdb', align_codes='mutant')
    h.res_num_from(m, aln)  # Restore old residue numbering and chain indexing
    h.write(file=mutid+'.pdb')

######################################################################################################################################################
# Function to run PDB2PQR.exe - should work on any supported OS
######################################################################################################################################################
def execPDB2PQR(path_pdb2pqr_exe, pdbfile, outfile=None, ff='parse'):
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
    # os.system('"{0}" {1} {2} {3}'.format(path_pdb2pqr_exe, optargs, pdbfile, outfile))
    (log, err) = runProcess([path_pdb2pqr_exe, '--ff=%s'%(ff), '--chain', pdbfile, outfile])
    return (log, err)

######################################################################################################################################################
# Function to run APBS.exe - should work on any supported OS
######################################################################################################################################################
def execAPBS(path_apbs_exe, pqr_chain, pqr_complex, prefix=None, grid=1.0, ion=0.150, pdie=20.0, sdie=78.54):
    """Summary
    
    Parameters
    ----------
    path_apbs_exe : STRING
        Full path to APBS executable, EX: 'C:\\APBS\\apbs.exe'
    pqr_chain : STRING
        PQR file name containing the segment that will undergo electrostatic calculations
    pqr_complex : STRING
        PQR file name containing the complex that AESOP is analyzing, must contain pqr_chain
    prefix : STRING, optional
        Phrase to prepend before any file that is generated before writing
    grid : float, optional
        Grid spacing for the mesh grid based electrostatic calculations. Suggested value of 1 or below
    ion : float, optional
        Ionic strength for APBS calculation
    pdie : float, optional
        Protein dielectric constant for APBS calculation
    sdie : float, optional
        Solvent dielectric constant for APBS calculation
    
    Returns
    -------
    file_apbs_log : STRING
        File name for the log file that APBS generates. This file contains results from calculations performed and must be parsed
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
    dime_ind = np.ceil(fg/(32*grid))-1    # index of dime to use from list, subtract one to be consistent with python indexing!

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
    # os.system('"{0}" {1} {2}'.format(path_apbs_exe, '--output-file=%s --output-format=flat'%(file_apbs_log), file_apbs_in))
    # os.system('{0} {1}'.format(path_apbs_exe, file_apbs_in))
    # (log, err) = runProcess([path_apbs_exe, file_apbs_in])
    (log, err) = runProcess([path_apbs_exe, '--output-file=%s'%(file_apbs_log), '--output-format=flat', file_apbs_in])
    pattern = re.compile('(?<=Global net ELEC energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?')
    elec = np.asarray([x.split() for x in re.findall(pattern, log)]).astype(np.float)
    elec = elec.reshape((1,elec.size))

    # return file_apbs_log
    return elec

######################################################################################################################################################
# Function to run coulomb.exe - should work on any supported OS
######################################################################################################################################################
def execCoulomb(path_coulomb_exe, pqr):
	(log, err) = runProcess([path_coulomb_exe, pqr])
	pattern = re.compile('(?<=Total energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?') #May need to update regex
	coul = np.asarray(re.findall(pattern, log)).astype(np.float)
	return coul

######################################################################################################################################################
# Function to plot results of Alascan
######################################################################################################################################################
def plotResults(Alascan, filename=None):
    plt.style.use('seaborn-talk') 
    figure, axarr = plt.subplots(len(Alascan.mutid) - 1, sharey=True)
    for i in xrange(1,len(Alascan.mutid)):
        axarr[i-1].set_title(np.unique(np.array([w.split('_') for w in Alascan.mutid[i]])[:,0])[0]+' ddGbind relative to WT')
        axarr[i-1].set_ylabel('kJ/mol')
        axarr[i-1].set_xticks(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]])))
        axarr[i-1].set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[i]])[:,1], rotation='vertical', ha='left')
        axarr[i-1].bar(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]))[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]] > 0], Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]  > 0], color = 'red' )
        axarr[i-1].bar(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]))[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]] < 0], Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]  < 0], color = 'blue' )
        axarr[i-1].xaxis.set_ticks_position('bottom')
        axarr[i-1].yaxis.set_ticks_position('left')
    plt.tight_layout()
    if filename is not None:
        figure.savefig(filename)

######################################################################################################################################################
# Function to parse APBS log file - REMOVED as it is not required!
######################################################################################################################################################
# def parseAPBS_totEnergy(path_log):
#     """Searches for a 'totEnergy' calculation result in the APBS log file
    
#     Parameters
#     ----------
#     path_log : STRING
#         Full path to APBS log file, EX: 'C:\\Users\\User\\Documents\\AESOP\\apbs.log'
    
#     Returns
#     -------
#     data : NDARRAY
#         Array that contains results of calculations in the log file, units should be kJ/mol
#     """
#     data = []
#     with open(path_log, 'r') as f:
#         lines = f.read()
#         # The following pattern extracts a scientific notation number only if preceded by totEnergy
#         # RegEx for scientific notation is: "[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?"
#         pattern = re.compile('(?<=totEnergy)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?')
#         matches = re.findall(pattern, lines)
#         data = np.asarray([x.split() for x in matches]).astype(np.float)
#         data = data.reshape((1,data.size))
#     return data

######################################################################################################################################################
# Dictionary to convert between 3 letter and 1 letter amino acid codes
######################################################################################################################################################
AA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

