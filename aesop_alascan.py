
import os as os
import sys as sys
# import glob as glob
import datetime as dt
import re as re
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
    def __init__(self, pdb, pdb2pqr_exe, apbs_exe, selstr=['protein'], jobname=None, region=None, grid=1, ion=0.150, pdie=20.0, sdie=78.54):
        self.pdb = pdb
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        self.selstr = selstr
        if region is None:
            # self.region = ''.join(['(',') or ('.join(selstr),')'])
            self.region = selstr
        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        # Insert code to instantiate dirs and prefix
        # self.dirs = 0
        #jobname must be alphanumeric, no spaces
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

    def getPDB(self):
        return self.pdb
    def getSel(self):
        return self.selstr
    def getDirs(self):
        return self.dirs
    def getPrefix(self):
        return self.prefix
    def getEnergies(self):
        return 0
    def getMutids(self):
        return self.mutids

    # def genMut(self):
    #     #Need to enforce that users either include a pdb with no chain IDs or all chain IDs, cannot have a comparison with a region that has no chain ID and one that does
    #     pdb = pd.parsePDB(self.pdb)
    #     selstr = self.selstr
    #     jobdir = self.jobdir
    #     jobname = self.jobname
    #     region = self.region

    #     list_of_pdb_dirs = []
    #     chainid_check = np.unique(pdb.select(''.join(['(',') or ('.join(selstr),')'])).getChids())[0].isspace()
    #     if chainid_check is True:
    #         complex_pdb_dir = jobname+'_pdb'
    #         if not os.path.exists(os.path.join(jobdir, complex_pdb_dir)):
    #             os.makedirs(os.path.join(jobdir, complex_pdb_dir))
    #         list_of_pdb_dirs.append(os.path.join(jobdir, complex_pdb_dir))
    #         complex_parent_pdb_file = os.path.join(jobdir, complex_pdb_dir, jobname+'.pdb')
    #         pd.writePDB(complex_parent_pdb_file, pdb.select(''.join(['(',') or ('.join(selstr),')'])))
    #     else:
    #         complex_pdb_dir = ''.join(('chain', '_chain'.join(np.unique(pdb.select(''.join(['(',') or ('.join(selstr),')'])).getChids())),'_complex_pdb'))
    #         if not os.path.exists(os.path.join(jobdir, complex_pdb_dir)):
    #             os.makedirs(os.path.join(jobdir, complex_pdb_dir))
    #         list_of_pdb_dirs.append(os.path.join(jobdir, complex_pdb_dir))
    #         complex_parent_pdb_file = os.path.join(jobdir, complex_pdb_dir, complex_pdb_dir.replace('_pdb','.pdb'))
    #         pd.writePDB(complex_parent_pdb_file, pdb.select(''.join(['(',') or ('.join(selstr),')'])))  # changed selections to selstr
    #         for i in selstr:
    #             indiv_pdb_dir = ''.join(('chain', '_chain'.join(np.unique(pdb.select(i).getChids())), '_pdb'))
    #             if not os.path.exists(os.path.join(jobdir, indiv_pdb_dir)):
    #                 os.makedirs(os.path.join(jobdir, indiv_pdb_dir))
    #             list_of_pdb_dirs.append(os.path.join(jobdir, indiv_pdb_dir))
    #             indiv_pdb_file = os.path.join(jobdir, indiv_pdb_dir, indiv_pdb_dir.replace('_pdb','.pdb'))
    #             pd.writePDB(indiv_pdb_file, pdb.select(i))

    #     for i,j in zip(selstr, region):
    #         combined_selection = pdb.select(''.join(['(',') and ('.join((i, j, 'charged', 'calpha')),')']))
    #         # combined_selection = pdb.select('(%s) and (%s)'%(i,j)) # Methods do not work, it can be much simpler though
    #         list_of_res_chainids = combined_selection.getChids().tolist()
    #         list_of_res_nums = map(str, combined_selection.getResnums().tolist())
    #         list_of_res_names = combined_selection.getResnames().tolist()
    #         if chainid_check is True:
    #             for res_no, res_id in zip(list_of_res_nums, list_of_res_names):
    #                 mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_'+AA_dict[res_id]+res_no+'A'
    #                 mutatePDB(pdb=complex_parent_pdb_file, mutid=mut_file_prefix, resnum=res_no, chain=None, resid='ALA')
    #         else:
    #             for ch_id, res_no, res_id in zip(list_of_res_chainids, list_of_res_nums, list_of_res_names):
    #                 mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
    #                 print mut_file_prefix
    #                 mutatePDB(pdb=complex_parent_pdb_file, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')
    #                 mut_indiv_file_dir = ''.join(('chain', '_chain'.join(np.unique(pdb.select(i).getChids())), '_pdb'))
    #                 mut_file_parent_pdb = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb','.pdb'))
    #                 mut_file_prefix = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb',''))+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
    #                 mutatePDB(pdb=mut_file_parent_pdb, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')


    def genPDB(self):
        parent_pdb = pd.parsePDB(self.pdb)

        # Create necessary directories for PDB files
        pdb_complex_dir = 'complex_pdb'
        if not os.path.exists(os.path.join(self.jobdir, pdb_complex_dir)):
            os.makedirs(os.path.join(self.jobdir, pdb_complex_dir))
        pdb_sel_dir=[]
        for i in xrange(0,len(self.selstr)):
            pdb_sel_dir.append('sel_%d_pdb'%(i))
            if not os.path.exists(os.path.join(self.jobdir, 'sel_%d_pdb'%(i))):
                os.makedirs(os.path.join(self.jobdir, 'sel_%d_pdb'%(i)))
        self.pdb_complex_dir = pdb_complex_dir
        self.pdb_sel_dir = pdb_sel_dir

        # Create necessary directories for PQR files
        pqr_complex_dir = 'complex_pqr'
        if not os.path.exists(os.path.join(self.jobdir, pqr_complex_dir)):
            os.makedirs(os.path.join(self.jobdir, pqr_complex_dir))
        pqr_sel_dir=[]
        for i in xrange(0,len(self.selstr)):
            pqr_sel_dir.append('sel_%d_pqr'%(i))
            if not os.path.exists(os.path.join(self.jobdir, 'sel_%d_pqr'%(i))):
                os.makedirs(os.path.join(self.jobdir, 'sel_%d_pqr'%(i)))
        self.pqr_complex_dir = pqr_complex_dir
        self.pqr_sel_dir = pqr_sel_dir

        # Create necessary directories for APBS files
        logs_apbs_dir = 'apbs_logs'
        if not os.path.exists(os.path.join(self.jobdir, logs_apbs_dir)):
            os.makedirs(os.path.join(self.jobdir, logs_apbs_dir))
        self.logs_apbs_dir = 'apbs_logs'

        # Set reference parent structure file names
        parent_file_prefix = 'wt'
        complex_parent_pdb_file = os.path.join(self.jobdir, pdb_complex_dir, parent_file_prefix+'.pdb')

        # Find all residues that must be mutated, generating mutids
        chainid_check = np.unique(pd.parsePDB(self.pdb).select(''.join(['(',') or ('.join(self.selstr),')'])).getChids())[0].isspace() # Check to see if there are no Chids
        list_mutids = [[] for x in xrange(len(self.selstr)+1)]
        counter = 0
        for id_by_sel in list_mutids:
            if counter == 0: # if counter is zero, then this is the wild type
                pdb = parent_pdb.select(' or '.join(self.selstr))
                pd.writePDB(complex_parent_pdb_file, pdb)
                id_by_sel.append(parent_file_prefix)
                for sel, sel_dir in zip(self.selstr, pdb_sel_dir):
                    pdb = parent_pdb.select(sel)
                    pd.writePDB(os.path.join(self.jobdir, sel_dir, parent_file_prefix+'_'+sel_dir.replace('_pdb','.pdb')), pdb)
                    id_by_sel.append(sel_dir.replace('_pdb','')+'_'+parent_file_prefix)

            else:
            # for sel, reg, sel_dir in zip(self.selstr, self.region, pdb_sel_dir):
                sel = self.selstr[counter-1]
                reg = self.region[counter-1]
                sel_dir = pdb_sel_dir[counter-1]
                combined_selection = parent_pdb.select(''.join(['(', ') and ('.join((sel, reg, 'charged', 'calpha')), ')']))
                list_of_res_chainids = combined_selection.getChids().tolist()
                list_of_res_nums = map(str, combined_selection.getResnums().tolist())
                list_of_res_names = combined_selection.getResnames().tolist()
                if chainid_check is True:
                    for res_no, res_id in zip(list_of_res_nums, list_of_res_names):
                        mut_file_prefix = sel_dir.replace('_pdb','')+'_'+AA_dict[res_id]+res_no+'A'
                        id_by_sel.append(mut_file_prefix)
                        path_mut_complex = os.path.join(self.jobdir, pdb_complex_dir, mut_file_prefix)
                        path_mut_chain = os.path.join(self.jobdir, sel_dir, mut_file_prefix)
                        mutatePDB(pdb=complex_parent_pdb_file, mutid=path_mut_complex, resnum=res_no, chain=None, resid='ALA')
                        pd.writePDB(path_mut_chain+'.pdb', pd.parsePDB(path_mut_complex+'.pdb').select(sel))
                else:
                    for ch_id, res_no, res_id in zip(list_of_res_chainids, list_of_res_nums, list_of_res_names):
                        mut_file_prefix = sel_dir.replace('_pdb','')+'_'+AA_dict[res_id]+res_no+'A'
                        id_by_sel.append(mut_file_prefix)
                        path_mut_complex = os.path.join(self.jobdir, pdb_complex_dir, mut_file_prefix)
                        path_mut_chain = os.path.join(self.jobdir, sel_dir, mut_file_prefix)
                        mutatePDB(pdb=complex_parent_pdb_file, mutid=path_mut_complex, resnum=res_no, chain=ch_id, resid='ALA')
                        pd.writePDB(path_mut_chain+'.pdb', pd.parsePDB(path_mut_complex+'.pdb').select(sel))
            counter += 1
        self.mutids = list_mutids

    def genPQR(self):
        # Convert all PDB files to PQR files
        path_in = os.path.join(self.jobdir, self.pdb_complex_dir)
        path_out = os.path.join(self.jobdir, self.pqr_complex_dir)
        for name in os.listdir(path_in):
            name = os.path.basename(name)
            execPDB2PQR(self.pdb2pqr, os.path.join(path_in, name), outfile=os.path.join(path_out, name.replace('.pdb', '.pqr')))

        # Write individual chains for mutants
        mutids = self.mutids
        for i in xrange(len(mutids)): # Loop for each selection: wt, sel0, sel1, etc.
            if i == 0:
                j = 1
                for sel in self.selstr:
                    path_pqr_complex = os.path.join(self.jobdir, self.pqr_complex_dir, mutids[i][0]+'.pqr')
                    path_pqr_sel = os.path.join(self.jobdir, self.pqr_sel_dir[j-1], mutids[i][j]+'.pqr')
                    pqr = pd.parsePQR(path_pqr_complex)
                    pd.writePQR(path_pqr_sel, pqr.select(sel))
                    j += 1
            else:
                for name in mutids[i]:
                    path_pqr_complex = os.path.join(self.jobdir, self.pqr_complex_dir, name+'.pqr')
                    path_pqr_sel = os.path.join(self.jobdir, self.pqr_sel_dir[i-1], name+'.pqr')
                    pqr = pd.parsePQR(path_pqr_complex)
                    pd.writePQR(path_pqr_sel, pqr.select(self.selstr[i-1]))

    def calcAPBS(self):
        mutids = self.mutids
        for i in xrange(len(mutids)): # Loop for each selection: wt, sel0, sel1, etc.
            if i == 0:
                j = 1
                for sel in self.selstr:
                    path_pqr_complex = os.path.join(self.jobdir, self.pqr_complex_dir, mutids[i][0]+'.pqr')
                    path_pqr_sel = os.path.join(self.jobdir, self.pqr_sel_dir[j-1], mutids[i][j]+'.pqr')
                    path_prefix = os.path.join(self.jobdir, self.logs_apbs_dir, mutids[i][j])
                    path_log = os.path.join(self.jobdir, self.logs_apbs_dir, mutids[i][j]+'.log')
                    execAPBS(self.apbs, path_pqr_sel, path_pqr_complex, prefix=path_prefix, grid=self.grid, ion=self.ion, pdie=self.pdie, sdie=self.sdie)
                    data = parseAPBS_totEnergy(path_log)
                    self.E_solv = np.append(self.E_solv, data[0][0])
                    self.E_ref = np.append(self.E_ref, data[0][1])
                    j += 1
            else:
                for name in mutids[i]:
                    path_pqr_complex = os.path.join(self.jobdir, self.pqr_complex_dir, name+'.pqr')
                    path_pqr_sel = os.path.join(self.jobdir, self.pqr_sel_dir[i-1], name+'.pqr')
                    path_prefix = os.path.join(self.jobdir, self.logs_apbs_dir, name)
                    path_log = os.path.join(self.jobdir, self.logs_apbs_dir, name+'.log')
                    execAPBS(self.apbs, path_pqr_sel, path_pqr_complex, prefix=path_prefix, grid=self.grid, ion=self.ion, pdie=self.pdie, sdie=self.sdie)
                    data = parseAPBS_totEnergy(path_log)
                    self.E_solv = np.append(self.E_solv, data[0][0])
                    self.E_ref = np.append(self.E_ref, data[0][1])

        # for ids, i in zip(mutids, xrange(len(mutids))):
        #     for name, j in zip(ids, xrange(len(mutids))):
        #         if i == 0:
        #             path_pqr_chain = os.path.join(self.jobdir, self.pqr_sel_dir, mutids[i][i])

        #         else:
        #             path_pqr_chain = os.path.join(self.jobdir, self.pqr_sel_dir[i])
        #             path_pqr_complex =
        #             path_prefix =
        #             execAPBS(path_apbs_exe, pqr_chain, pqr_complex, prefix=name, grid=self.grid, ion=self.ion, pdie=self.pdie, sdie=self.sdie)

    def run(self):
        self.genPDB()
        self.genPQR()
        self.calcAPBS()






        #         combined_selection = pdb.select(''.join(['(',') and ('.join((i, j, 'charged', 'calpha')),')']))
        #         list_of_res_chainids = combined_selection.getChids().tolist()
        #         list_of_res_nums = map(str, combined_selection.getResnums().tolist())
        #         list_of_res_names = combined_selection.getResnames().tolist()
        #         if chainid_check is True:
        #             for res_no, res_id in zip(list_of_res_nums, list_of_res_names):
        #                 mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_'+AA_dict[res_id]+res_no+'A'
        #                 mutatePDB(pdb=complex_parent_pdb_file, mutid=mut_file_prefix, resnum=res_no, chain=None, resid='ALA')
        #         else:
        #             for ch_id, res_no, res_id in zip(list_of_res_chainids, list_of_res_nums, list_of_res_names):
        #                 mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
        #                 mutatePDB(pdb=self.pdb, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')
        #                 mut_indiv_file_dir = ''.join(('chain', '_chain'.join(np.unique(pdb.select(i).getChids())), '_pdb'))
        #                 mut_file_parent_pdb = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb','.pdb'))
        #                 mut_file_prefix = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb',''))+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
        #                 mutatePDB(pdb=self.pdb, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')





        # for i,j, list_id in zip(selstr, region):

        #     list_mutfiles

        #     combined_selection = pdb.select(''.join(['(',') and ('.join((i, j, 'charged', 'calpha')),')']))
        #     list_of_res_chainids = combined_selection.getChids().tolist()
        #     list_of_res_nums = map(str, combined_selection.getResnums().tolist())
        #     list_of_res_names = combined_selection.getResnames().tolist()
        #     if chainid_check is True:
        #         for res_no, res_id in zip(list_of_res_nums, list_of_res_names):
        #             mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_'+AA_dict[res_id]+res_no+'A'
        #             mutatePDB(pdb=complex_parent_pdb_file, mutid=mut_file_prefix, resnum=res_no, chain=None, resid='ALA')
        #     else:
        #         for ch_id, res_no, res_id in zip(list_of_res_chainids, list_of_res_nums, list_of_res_names):
        #             mut_file_prefix = complex_parent_pdb_file.replace('.pdb','')+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
        #             mutatePDB(pdb=complex_parent_pdb_file, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')
        #             mut_indiv_file_dir = ''.join(('chain', '_chain'.join(np.unique(pdb.select(i).getChids())), '_pdb'))
        #             mut_file_parent_pdb = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb','.pdb'))
        #             mut_file_prefix = os.path.join(jobdir, mut_indiv_file_dir, mut_indiv_file_dir.replace('_pdb',''))+'_mutchain'+ch_id+'_'+AA_dict[res_id]+res_no+'A'
        #             mutatePDB(pdb=mut_file_parent_pdb, mutid=mut_file_prefix, resnum=res_no, chain=ch_id, resid='ALA')


        # Generate all PDB files, including mutants
        # Run APBS on all structures, storing results


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
    os.system('"{0}" {1} {2} {3}'.format(path_pdb2pqr_exe, optargs, pdbfile, outfile))
    return outfile

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
    os.system('"{0}" {1} {2}'.format(path_apbs_exe, '--output-file=%s --output-format=flat'%(file_apbs_log), file_apbs_in))
    # os.system('{0} {1}'.format(path_apbs_exe, file_apbs_in))

    return file_apbs_log

######################################################################################################################################################
# Function to parse APBS log file
######################################################################################################################################################
def parseAPBS_totEnergy(path_log):
    """Searches for a 'totEnergy' calculation result in the APBS log file
    
    Parameters
    ----------
    path_log : STRING
        Full path to APBS log file, EX: 'C:\\Users\\User\\Documents\\AESOP\\apbs.log'
    
    Returns
    -------
    data : NDARRAY
        Array that contains results of calculations in the log file, units should be kJ/mol
    """
    data = []
    with open(path_log, 'r') as f:
        lines = f.read()
        # The following pattern extracts a scientific notation number only if preceded by totEnergy
        # RegEx for scientific notation is: "[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?"
        pattern = re.compile('(?<=totEnergy)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?')
        matches = re.findall(pattern, lines)
        data = np.asarray([x.split() for x in matches]).astype(np.float)
        data = data.reshape((1,data.size))
    return data

######################################################################################################################################################
# Dictionary to convert between 3 letter and 1 letter amino acid codes
######################################################################################################################################################
AA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
