import os as os
# import sys as sys
import subprocess as sp
import datetime as dt
import timeit as ti
import re as re
import numpy as np
import prody as pd
import scipy.spatial as spatial
# import scipy.interpolate as interp
import scipy.cluster.hierarchy as cluster
import matplotlib.pyplot as plt
from modeller import environ, model, alignment, selection
from multiprocessing import Pool#, freeze_support
import gridData as gd
import itertools as it
# import plotly.plotly as py
# import plotly.graph_objs as go
# from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
# from plotly.tools import FigureFactory as FF

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

    def __init__(self, pdb, pdb2pqr_exe, apbs_exe, coulomb_exe=None, selstr=['protein'], jobname=None, region=None,
                 grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse', cfac=1.5, dx=False):
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
            self.jobname = '%4d%02d%02d_%02d%02d%02d' % (
            dt.date.today().year, dt.date.today().month, dt.date.today().day, dt.datetime.now().hour,
            dt.datetime.now().minute, dt.datetime.now().second)
        else:
            self.jobname = jobname
        self.jobdir = jobname
        if not os.path.exists(os.path.join(self.jobdir)):
            os.makedirs(os.path.join(self.jobdir))
        self.E_ref = np.zeros(0)
        self.E_solv = np.zeros(0)
        self.mutid = []
        self.ff = ff
        self.cfac = cfac
        self.dx = dx

    def getPDB(self):
        return self.pdb

    def getSel(self):
        return self.selstr

    def getDirs(self):
        return self.dirs

    def getPrefix(self):
        return self.prefix

    def getMutids(self):
        l = self.mutid
        return [item for sublist in l for item in sublist]

    def getDX(self):
        return self.dx_files

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
        pqr_sel_dir = []
        for i in xrange(0, len(self.selstr)):
            pqr_sel_dir.append('seg%d_pqr' % (i + 1))
            if not os.path.exists(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1))):
                os.makedirs(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1)))
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

        list_mutids = [[] for x in xrange(len(selstr) + 1)]
        list_chids = [[] for x in xrange(len(selstr) + 1)]
        list_resnums = [[] for x in xrange(len(selstr) + 1)]
        list_resnames = [[] for x in xrange(len(selstr) + 1)]

        list_mutids[0] = [parent_file_prefix]
        list_chids[0] = ['']
        list_resnums[0] = [np.zeros(0)]
        list_resnames[0] = ['']

        index = np.linspace(1, len(selstr), len(selstr)).astype(int)
        for i, sel, reg in zip(index, selstr, region):
            # print ' and '.join([sel, reg, 'charged', 'calpha'])
            combined_selection = parent_pdb.select(' and '.join([sel, reg, 'charged', 'calpha']))
            # if sel is not reg:
            #     combined_selection = parent_pdb.select(''.join(['(', ') and ('.join((sel, reg, 'charged', 'calpha')), ')']))
            # elif sel is reg:
            #     combined_selection = parent_pdb.select(''.join(['(', ') and ('.join((sel, 'charged', 'calpha')), ')']))
            list_chids[i] = combined_selection.getChids().tolist()
            list_resnums[i] = combined_selection.getResnums().tolist()
            list_resnames[i] = combined_selection.getResnames().tolist()
            code = ['seg%d' % (i) + '_' + AA_dict[res_id] + res_no + 'A' for ch_id, res_no, res_id in
                    zip(list_chids[i], map(str, list_resnums[i]), list_resnames[i])]
            list_mutids[i] = code

        dim_sel = len(list_mutids)
        dim_mut = np.sum([len(x) for x in list_mutids])
        mask_by_sel = np.zeros((dim_mut, dim_sel)).astype(bool)
        counter = 0
        for i in xrange(dim_sel):
            for j in xrange(counter, counter + len(list_mutids[i]), 1):
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

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix + '.pdb')
        system = parent_pdb.select('((' + ') or ('.join(selstr) + '))')
        # print '(('+') or ('.join(selstr)+'))'
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

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix + '.pdb')
        system = parent_pdb.select('((' + ') or ('.join(selstr) + '))')
        # print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)

        for mutid, chain, resnum in zip(list_mutids[1:], list_chids[1:], list_resnums[1:]):
            outpath = os.path.join(jobdir, pdb_complex_dir, mutid)
            print '\n%s:\tgenerating PDB for mutant: %s' % (self.jobname, mutid)
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
        infile = os.path.join(jobdir, pdb_complex_dir, list_mutids[0] + '.pdb')
        outfile = os.path.join(jobdir, pqr_complex_dir, list_mutids[0] + '.pqr')
        print '\n%s:\tgenerating PQR for parent: %s' % (self.jobname, list_mutids[0])
        execPDB2PQR(path_pdb2pqr, infile, outfile=outfile, ff=ff)
        complex_pqr = pd.parsePQR(outfile)
        for sel, seldir in zip(selstr, pqr_sel_dir):
            selfile = os.path.join(jobdir, seldir, list_mutids[0] + '.pqr')
            pqr = complex_pqr.select(sel)
            pd.writePQR(selfile, pqr)

        for mutid, chain, resnum in zip(list_mutids[1:], list_chids[1:], list_resnums[1:]):
            outpath = os.path.join(jobdir, pqr_complex_dir, mutid)
            print '\n%s:\tgenerating PQR for mutant: %s' % (self.jobname, mutid)
            # print 'mutid %s, chain %s, resnum %d'%(mutid, chain, resnum)
            # print outpath+'.pqr'
            mutatePQR(outfile, mutid=outpath, resnum=resnum, chain=chain)
            complex_pqr = pd.parsePQR(outpath + '.pqr')
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid + '.pqr')
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
            print '\n%s:\tgenerating PQR for mutant: %s' % (self.jobname, mutid)
            infile = os.path.join(jobdir, pdb_complex_dir, mutid + '.pdb')
            outfile = os.path.join(jobdir, pqr_complex_dir, mutid + '.pqr')
            execPDB2PQR(path_pdb2pqr, infile, outfile=outfile, ff=ff)
            complex_pqr = pd.parsePQR(outfile)
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid + '.pqr')
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
        dim_sel = len(selstr) + 1

        mask_by_sel = np.copy(self.mask_by_sel)  # Mask parts that are true will be run with APBS
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel))
        Gref = np.zeros((dim_mutid, dim_sel))

        complex_pqr = os.path.join(jobdir, pqr_complex_dir, list_mutids[0] + '.pqr')
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating solvation and reference energies for mutant: %s' % (self.jobname, mutid)
            # complex_pqr = os.path.join(jobdir, pqr_complex_dir, mutid+'.pqr')
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_prefix_log = os.path.join(jobdir, logs_apbs_dir, mutid)
                if mask_by_sel[i, j]:
                    energies = execAPBS(path_apbs, subunit_pqr, complex_pqr, prefix=path_prefix_log, grid=self.grid,
                                        ion=self.ion, pdie=self.pdie, sdie=self.sdie, cfac=self.cfac)
                    # print energies[0][0]
                    # print energies[0][1]
                    Gsolv[i, j] = energies[0][0]
                    Gref[i, j] = energies[0][1]
                if not mask_by_sel[i, j]:
                    Gsolv[i, j] = Gsolv[0, j]
                    Gref[i, j] = Gref[0, j]
        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcAPBS_parallel(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        mask_by_sel = np.copy(self.mask_by_sel)  # Mask parts that are true will be run with APBS
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel)).astype(float)
        Gref = np.zeros((dim_mutid, dim_sel)).astype(float)

        path_list = []
        pqr_chain_list = []
        pqr_complex_list = []
        prefix_list = []
        grid_list = []
        ion_list = []
        pdie_list = []
        sdie_list = []
        cfac_list = []
        dx_list = []
        i_list = []
        j_list = []

        # Find all calculations to be done
        complex_pqr = os.path.join(jobdir, pqr_complex_dir, list_mutids[0] + '.pqr')
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                if mask_by_sel[i, j]:
                    path_list.append(path_apbs)
                    pqr_chain_list.append(subunit_pqr)
                    pqr_complex_list.append(complex_pqr)
                    prefix_list.append(os.path.join(jobdir, logs_apbs_dir, '%d_%d_' % (i, j)+mutid))  # added to make sure apbs.in file is unique!
                    grid_list.append(self.grid)
                    ion_list.append(self.ion)
                    pdie_list.append(self.pdie)
                    sdie_list.append(self.sdie)
                    cfac_list.append(self.cfac)
                    dx_list.append(self.dx)
                    i_list.append(i)
                    j_list.append(j)

        # Organize kernel and run batch process
        kernel = zip(path_list, pqr_chain_list, pqr_complex_list, prefix_list, grid_list, ion_list, pdie_list,
                     sdie_list, cfac_list, dx_list, i_list, j_list)
        apbs_results = []
        p = Pool()
        print '%s:\trunning batchAPBS ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchAPBS, kernel):
            counter += 1
            print '.... %s:\tbatch APBS %d percent complete ....' % (self.jobname, int(counter * 100 / max_count))
            i = int(result[0])
            j = int(result[1])
            solv = result[2]
            ref = result[3]
            # print '%d, %d, %f, %f'%(i, j, solv, ref)
            apbs_results.append([i, j, solv, ref])
            Gsolv[i, j] = solv
            Gref[i, j] = ref
        apbs_results = np.asarray(apbs_results)
        self.apbs_results = apbs_results
        if self.dx == True:
            self.dx_files = [x+'.dx' for x in prefix_list]

        # Fill in results that are duplicates
        for i in xrange(dim_mutid):
            for j in xrange(dim_sel):
                if not mask_by_sel[i, j]:
                    Gsolv[i, j] = Gsolv[0, j]
                    Gref[i, j] = Gref[0, j]

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
        dim_sel = len(selstr) + 1

        Gcoul = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating coulombic energies for mutant: %s' % (self.jobname, mutid)
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                energies = execCoulomb(path_coulomb, subunit_pqr)
                Gcoul[i, j] = energies / pdie

        self.Gcoul = Gcoul

    def calcCoulomb_parallel(self):
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_coulomb = self.coulomb

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        mask_by_sel = np.copy(self.mask_by_sel)  # Mask parts that are true will be run with APBS
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gcoul = np.zeros((dim_mutid, dim_sel)).astype(float)

        path_list = []
        pqr_chain_list = []
        pdie_list = []
        i_list = []
        j_list = []

        # Find all calculations to be done
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                if mask_by_sel[i, j]:
                    path_list.append(path_coulomb)
                    pqr_chain_list.append(subunit_pqr)
                    pdie_list.append(self.pdie)
                    i_list.append(i)
                    j_list.append(j)

        # Organize kernel and run batch process
        kernel = zip(path_list, pqr_chain_list, pdie_list, i_list, j_list)
        coulomb_results = []
        p = Pool()
        print '%s:\trunning batchCoulomb ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchCoulomb, kernel):
            counter += 1
            print '.... %s:\tbatch coulomb %d percent complete ....' % (self.jobname, int(counter * 100 / max_count))
            i = int(result[0])
            j = int(result[1])
            coul = result[2]
            coulomb_results.append([i, j, coul])
            Gcoul[i, j] = coul
        coulomb_results = np.asarray(coulomb_results)
        self.coulomb_results = coulomb_results

        # Fill in results that are duplicates
        for i in xrange(dim_mutid):
            for j in xrange(dim_sel):
                if not mask_by_sel[i, j]:
                    Gcoul[i, j] = Gcoul[0, j]
        self.Gcoul = Gcoul

    def ddGbind_rel(self):
        Gsolv = self.Gsolv
        Gref = self.Gref
        Gcoul = self.Gcoul

        dGsolv = Gsolv - Gref

        dGsolu = Gsolv[:, 0] - Gsolv[:, 1:].sum(axis=1)
        dGcoul = Gcoul[:, 0] - Gcoul[:, 1:].sum(axis=1)
        ddGsolv = dGsolv[:, 0] - dGsolv[:, 1:].sum(axis=1)

        dGbind = ddGsolv + dGcoul
        dGbind_rel = dGbind - dGbind[0]
        return dGbind_rel

    def dGsolv_rel(self):
        Gsolv = self.Gsolv
        Gref = self.Gref
        dGsolv = Gsolv - Gref
        dGsolv = dGsolv - dGsolv[0, 0]
        return dGsolv[:, 0]

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

    def run_parallel(self):
        start = ti.default_timer()
        self.genDirs()
        self.genMutid()  # contains warning: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison if tokens[0] == 'and' or tokens[-1] == 'and':
        self.genParent()
        self.genTruncatedPQR()
        self.calcAPBS_parallel()
        self.calcCoulomb_parallel()
        stop = ti.default_timer()
        print '%s:\tAESOP alanine scan completed in %.2f seconds' % (self.jobname, stop - start)

    def summary(self, filename=None):
        plotResults(self, filename=None)

######################################################################################################################################################
# Container for performing an Directed Mutagenesis Scan with AESOP
#   pdb     -   PDB file for performing Alascan. Must contain all chain selections with standard amino acid nomenclature
#   selstr  -   List of selection strings for chain selection. This does not change what is mutated.
#   target  -   List of selection strings for mutants, each will element of list will be mutated simultaneously
#   mutation -  Mutate each target residue to the corresponding element of this list. (3 letter code)
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant
#   sdie    -   Solvent dielectric constant
######################################################################################################################################################
class DirectedMutagenesis:
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

    def __init__(self, pdb, target, mutation, pdb2pqr_exe, apbs_exe, coulomb_exe=None, selstr=['protein'], jobname=None,
                 grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse', cfac=1.5, dx=False):
        self.pdb = pdb
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        if coulomb_exe is None:
            self.coulomb = os.path.split(apbs_exe)[0]
        else:
            self.coulomb = coulomb_exe
        self.selstr = selstr
        self.target = target
        self.mutation = mutation
        # if region is None:
        #     self.region = selstr
        # else:
        #     self.region = region
        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        if jobname is None:
            self.jobname = '%4d%02d%02d_%02d%02d%02d' % (
            dt.date.today().year, dt.date.today().month, dt.date.today().day, dt.datetime.now().hour,
            dt.datetime.now().minute, dt.datetime.now().second)
        else:
            self.jobname = jobname
        self.jobdir = jobname
        if not os.path.exists(os.path.join(self.jobdir)):
            os.makedirs(os.path.join(self.jobdir))
        self.E_ref = np.zeros(0)
        self.E_solv = np.zeros(0)
        self.mutid = []
        self.ff = ff
        self.cfac = cfac
        self.dx = dx

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
        pqr_sel_dir = []
        for i in xrange(0, len(self.selstr)):
            pqr_sel_dir.append('seg%d_pqr' % (i + 1))
            if not os.path.exists(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1))):
                os.makedirs(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1)))
        self.pqr_complex_dir = pqr_complex_dir
        self.pqr_sel_dir = pqr_sel_dir

        # Create necessary directories for APBS files
        logs_apbs_dir = 'apbs_logs'
        if not os.path.exists(os.path.join(self.jobdir, logs_apbs_dir)):
            os.makedirs(os.path.join(self.jobdir, logs_apbs_dir))
        self.logs_apbs_dir = 'apbs_logs'

    def genMutid(self):
        selstr = self.selstr
        target = self.target
        mutation = self.mutation

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        list_mutids = [[] for x in xrange(len(selstr) + 1)]
        list_chids = [[] for x in xrange(len(selstr) + 1)]
        list_resnums = [[] for x in xrange(len(selstr) + 1)]
        list_resnames = [[] for x in xrange(len(selstr) + 1)]

        dim_sel = len(selstr) + 1
        dim_mut = len(target) + 1
        mask_by_sel = np.zeros((dim_mut, dim_sel)).astype(bool)
        mask_by_sel[0, 0] = True

        list_mutids[0] = [parent_file_prefix]
        list_chids[0] = ['']
        list_resnums[0] = [np.zeros(0).tolist()]
        list_resnames[0] = ['']

        for i, region in zip(xrange(len(selstr)), selstr):
            for j, sel, mut in zip(xrange(len(target)), target, mutation):
                combined_selection = parent_pdb.select('('+' and '.join([sel, region, 'calpha'])+')')

                if combined_selection is not None:
                    chids = combined_selection.getChids().tolist()
                    resnums = combined_selection.getResnums().tolist()
                    resnames = combined_selection.getResnames().tolist()

                    label = 'sel%d' % (i+1)
                    for resnum, resname in zip(resnums, resnames):
                        label = '_'.join([label, AA_dict[resname]+str(resnum)+AA_dict[mut]])

                    list_mutids[i+1].append(label)
                    list_chids[i+1].append(chids)
                    list_resnums[i+1].append(resnums)
                    list_resnames[i+1].append(resnames)
                    mask_by_sel[j+1, i+1] = True

        self.mutid = list_mutids
        self.list_chids = list_chids
        self.list_resnums = list_resnums
        self.list_resnames = list_resnames
        self.mask_by_sel = mask_by_sel

    def genParent(self):
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        # list_mutids = [item for sublist in self.mutid for item in sublist]
        # list_chids = [item for sublist in self.list_chids for item in sublist]
        # list_resnums = [item for sublist in self.list_resnums for item in sublist]

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix + '.pdb')
        system = parent_pdb.select('((' + ') or ('.join(selstr) + '))')
        # print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)

    def genPDB(self):
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        mutation = self.mutation

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        list_mutids = [item for sublist in self.mutid for item in sublist]
        list_chids = [item for sublist in self.list_chids for item in sublist]
        list_resnums = [item for sublist in self.list_resnums for item in sublist]

        infile = os.path.join(jobdir, pdb_complex_dir, parent_file_prefix + '.pdb')
        system = parent_pdb.select('(' + ') or ('.join(selstr) + ')')
        pd.writePDB(infile, system)

        for mutid, chain, resnum, mut in zip(list_mutids[1:], list_chids[1:], list_resnums[1:], mutation):
            outpath = os.path.join(jobdir, pdb_complex_dir, mutid)
            print '\n%s:\tgenerating PDB for mutant: %s' % (self.jobname, mutid)
            mutatePDB(pdb=infile, mutid=outpath, resnum=resnum, chain=chain, resid=mut)

    def genPQR(self):
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_pdb2pqr = self.pdb2pqr
        ff = self.ff

        list_mutids = [item for sublist in self.mutid for item in sublist]

        # Calculate PQR of complexes
        for mutid in list_mutids:
            print '\n%s:\tgenerating PQR for mutant: %s' % (self.jobname, mutid)
            infile = os.path.join(jobdir, pdb_complex_dir, mutid + '.pdb')
            outfile = os.path.join(jobdir, pqr_complex_dir, mutid + '.pqr')
            execPDB2PQR(path_pdb2pqr, infile, outfile=outfile, ff=ff)
            complex_pqr = pd.parsePQR(outfile)
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid + '.pqr')
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)

    def calcAPBS(self):
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        mask_by_sel = np.copy(self.mask_by_sel)  # Mask parts that are true will be run with APBS
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel))
        Gref = np.zeros((dim_mutid, dim_sel))

        complex_pqr = os.path.join(jobdir, pqr_complex_dir, list_mutids[0] + '.pqr')
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating solvation and reference energies for mutant: %s' % (self.jobname, mutid)
            # complex_pqr = os.path.join(jobdir, pqr_complex_dir, mutid+'.pqr')
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_prefix_log = os.path.join(jobdir, logs_apbs_dir, mutid)
                if mask_by_sel[i, j]:
                    energies = execAPBS(path_apbs, subunit_pqr, complex_pqr, prefix=path_prefix_log, grid=self.grid,
                                        ion=self.ion, pdie=self.pdie, sdie=self.sdie, cfac=self.cfac)
                    # print energies[0][0]
                    # print energies[0][1]
                    Gsolv[i, j] = energies[0][0]
                    Gref[i, j] = energies[0][1]
                if not mask_by_sel[i, j]:
                    Gsolv[i, j] = Gsolv[0, j]
                    Gref[i, j] = Gref[0, j]
        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcCoulomb(self):
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_coulomb = self.coulomb
        pdie = self.pdie

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        Gcoul = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating coulombic energies for mutant: %s' % (self.jobname, mutid)
            for j, seldir in zip(xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                energies = execCoulomb(path_coulomb, subunit_pqr)
                Gcoul[i, j] = energies / pdie

        self.Gcoul = Gcoul

    def ddGbind_rel(self):
        Gsolv = self.Gsolv
        Gref = self.Gref
        Gcoul = self.Gcoul

        dGsolv = Gsolv - Gref

        dGsolu = Gsolv[:, 0] - Gsolv[:, 1:].sum(axis=1)
        dGcoul = Gcoul[:, 0] - Gcoul[:, 1:].sum(axis=1)
        ddGsolv = dGsolv[:, 0] - dGsolv[:, 1:].sum(axis=1)

        dGbind = ddGsolv + dGcoul
        dGbind_rel = dGbind - dGbind[0]
        return dGbind_rel

    def dGsolv_rel(self):
        Gsolv = self.Gsolv
        Gref = self.Gref
        dGsolv = Gsolv - Gref
        dGsolv = dGsolv - dGsolv[0, 0]
        return dGsolv[:, 0]

    def run(self):
        self.genDirs()
        self.genMutid()
        self.genParent()
        self.genPDB()
        self.genPQR()
        self.calcAPBS()
        self.calcCoulomb()

######################################################################################################################################################
# Container for performing ESD analysis on set of PDB files
#   alascan     -   Alascan class with certain class functions.
######################################################################################################################################################
class ElecSimilarity: # PLEASE SUPERPOSE SYSTEM BEFORE USING THIS METHOD! Coordinates must be consistent!
    def __init__(self, pdbfiles, pdb2pqr_exe, apbs_exe, selstr=None, jobname=None,
                 grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse', cfac=1.5):
        self.pdbfiles = pdbfiles
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        self.dx = True
        if jobname is None:
            self.jobname = '%4d%02d%02d_%02d%02d%02d' % (dt.date.today().year, dt.date.today().month,
                                                         dt.date.today().day, dt.datetime.now().hour,
                                                         dt.datetime.now().minute, dt.datetime.now().second)
        else:
            self.jobname = jobname
        self.jobdir = jobname
        if not os.path.exists(os.path.join(self.jobdir)):
            os.makedirs(os.path.join(self.jobdir))
        self.pdbdir = os.path.join(self.jobdir, 'pdb_files')
        self.pqrdir = os.path.join(self.jobdir, 'pqr_files')
        if not os.path.exists(self.pdbdir):
            os.makedirs(self.pdbdir)
        if not os.path.exists(self.pqrdir):
            os.makedirs(self.pqrdir)

        for pdbfile in pdbfiles:
            pdb = pd.parsePDB(pdbfile)
            if selstr is None:
                pd.writePDB(os.path.join(self.pdbdir, os.path.basename(pdbfile)), pdb)
            elif selstr is not None:
                pd.writePDB(os.path.join(self.pdbdir, os.path.basename(pdbfile)), pdb)
        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        self.ff = ff
        self.cfac = cfac

    def findGLEN(self):
        0
        # Determine mesh dimensions according to Ron's AESOP protocol in the R source file
        # pdbfiles = self.pdbfiles
        # grid = self.grid
        # cfac = self.cfac
        # glen = np.zeros(1, 3)
        # for pdbfile in pdbfiles
        #     pdb = pd.parsePDB(pdbfile)
        #     coords = pdb.getCoords()
        #     x = coords[:,0]
        #     y = coords[:,1]
        #     z = coords[:,2]
        #     fg = np.array((np.ceil(np.max(x) - np.min(x)), np.ceil(np.max(y) - np.min(y)), np.ceil(np.max(z) - np.min(z))))
        #     fg = np.ceil((fg + 5) * cfac)
        #     glen = np.vstack((glen, fg)).max(axis=0)
        # dime_list = (32 * np.linspace(1, 100, 100)) + 1  # list of possible dime values
        # dime_ind = np.ceil(glen / (32 * grid)) - 1  # index of dime to use from list, subtract one to be consistent with python indexing!
        # dime = np.array((dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])], dime_list[int(dime_ind[2])]))
        #
        # self.dime = dime
        # self.glen = glen.reshape((1, 3))





######################################################################################################################################################
# Container for performing ESD analysis
#   alascan     -   Alascan class with certain class functions.
######################################################################################################################################################
class ESD:
    def __init__(self, alascan):
        dx_files = alascan.getDX() # Get only files with index [0-9]+.
        pattern_id = re.compile('\d+_0_(wt|seg\d+_[A-Z]\d+[A-Z])[.]dx')
        pattern_file = re.compile('(.*[0-9]+_0_(wt|seg[0-9]+_[A-Z]\d+[A-Z])[.]dx)')
        self.ids = [f.group(1) for dx_file in dx_files for f in [re.search(pattern_id, os.path.basename(dx_file))] if f]
        self.files = [f.group(1) for dx_file in dx_files for f in [re.search(pattern_file, dx_file)] if f]

        self.ion = alascan.ion
        self.pdie = alascan.pdie
        self.sdie = alascan.sdie
        self.jobname = alascan.jobname
        self.pdbfile = alascan.pdb

        grid = gd.Grid(self.files[0])
        self.midpoints = grid.midpoints
        self.edges = grid.edges
        self.dim_dx = grid.grid.shape

        self.mask = np.ones((self.dim_dx[0]*self.dim_dx[1]*self.dim_dx[2])).astype(bool)
        # num_files = len(self.dx_files)
        # dim_coords = grid.grid.shape[0] * grid.grid.shape[1] * grid.grid.shape[2]
        # self.coords = np.zeros((num_files, dim_coords)) # currently can't allocate enough memory
        # for i in xrange(num_files):
        #     grid = gd.Grid(self.dx_files[i])
        #     self.coords[i,:] = grid.grid.reshape((1, dim_coords))

    def findSurfaceGridPts(self, path_dssp):
        pdbfile = self.pdbfile
        pdb = pd.parsePDB(pdbfile)
        xyz = pdb.getCoords()
        hull = spatial.Delaunay(xyz)
        # resid, rsa, exposed = calcRSA(pdbfile, path_dssp)
        # xyz = xyz[exposed,:]

        dim = (self.dim_dx[0]*self.dim_dx[1]*self.dim_dx[2]/3, 3)
        grid = gd.Grid(self.files[0]).grid.reshape(dim)
        # result = interp.LinearNDInterpolator(hull, grid)
        # result = spatial.distance.cdist(xyz, grid)
        # return result




    def setMask(self, mask):
        if len(mask) == len(self.mask):
            self.mask = mask
        else:
            print 'Error: unable to set mask, must be an array of length %d' % (len(self.mask))

    def calc(self, method='LD'):

        def symmetrize(a):
            return a + a.T - np.diag(a.diagonal())

        files = self.files
        ids = self.ids
        dim = self.dim_dx[0]*self.dim_dx[1]*self.dim_dx[2]/3
        esd = np.zeros((len(ids), len(ids)))

        indices = it.combinations(range(len(ids)), 2)
        for i, j in indices:
            a = gd.Grid(files[i]).grid.reshape((dim, 3))
            b = gd.Grid(files[j]).grid.reshape((dim, 3))
            if method is 'LD':
                numer = np.linalg.norm(a-b, axis=1)
                denom = dim * np.max(np.hstack((np.linalg.norm(a, axis=1).reshape((dim, 1)), np.linalg.norm(b, axis=1).reshape((dim, 1)))), axis=1)
                esd[i, j] = np.divide(numer, denom).sum()
                print esd[i,j]
        esd = symmetrize(esd)
        self.esd = esd

        # for file, i in zip(files, xrange(len(ids))):
        #     a = gd.Grid(file).grid
        #     a = a.reshape((dim, 3))
        #     for file, j in zip(files, xrange(len(ids))):
        #         b = gd.Grid(file).grid
        #         b = b.reshape((dim, 3))
        #
        #         if method is 'LD':
        #             numer = np.linalg.norm(a-b, axis=1)
        #             denom = dim * np.max(np.hstack((np.linalg.norm(a, axis=1).reshape((dim, 1)), np.linalg.norm(b, axis=1).reshape((dim, 1)))), axis=1)
        #             esd[i, j] = np.divide(numer, denom).sum()
        # self.esd = esd

    def calc_batch(self, method='LD'): # NOT WORKING ... DON'T USE!

        def symmetrize(a):
            return a + a.T - np.diag(a.diagonal())

        files = self.files
        ids = self.ids
        dim = self.dim_dx[0]*self.dim_dx[1]*self.dim_dx[2]/3

        esd = np.zeros((len(ids), len(ids)))
        indices = it.combinations(range(len(ids)), 2)
        list_i = []
        list_j = []
        list_i_files = []
        list_j_files = []
        for i, j in indices:
            list_i.append(i)
            list_j.append(j)
            list_i_files.append(files[i])
            list_j_files.append(files[j])
        kernel = zip(list_i, list_j, list_i_files, list_j_files)
        print'Starting batch ESD calculation ....'
        p = Pool()
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(f_ld, kernel):
            counter += 1
            print '.... Batch ESD %d percent complete ....' % (int(counter * 100 / max_count))
            i = int(result[0])
            j = int(result[1])
            x = result[2]
            esd[i,j] = x
            print '%d, %d, %f' %(i,j,x)
        esd = symmetrize(esd)
        self.esd = esd

######################################################################################################################################################
# Function to run batch esd with LD method
######################################################################################################################################################
def f_ld(kernel):
    i, j, file_i, file_j = kernel
    a = gd.Grid(file_i).grid
    b = gd.Grid(file_j).grid
    dim = a.shape[0] * a.shape[1] * a.shape[2] / 3
    a = a.reshape((dim, 3))
    b = b.reshape((dim, 3))
    numer = np.linalg.norm(a-b, axis=1)
    denom = dim * np.max(np.hstack((np.linalg.norm(a, axis=1).reshape((dim, 1)), np.linalg.norm(b, axis=1).reshape((dim, 1)))), axis=1)
    esd = np.divide(numer, denom).sum()
    return np.array([i, j, esd])

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
def mutatePQR(pqrfile, mutid, resnum, chain=None):  # Only use this function with PARSE for now ...
    parent = pd.parsePQR(pqrfile)
    if chain is None:
        residue = parent.select('resnum %d' % (int(resnum)))
        preceed = parent.select('resnum < %d' % (int(resnum)))
        follow = parent.select('resnum > %d' % (int(resnum)))
    elif chain is not None:
        residue = parent.select('chain %s and resnum %d' % (str(chain), int(resnum)))
        preceed = parent.select('chain %s and resnum < %d' % (str(chain), int(resnum)))
        follow = parent.select('chain %s and resnum > %d' % (str(chain), int(resnum)))
        otherchains = parent.select('not chain %s' % (str(chain)))
    bb = residue.select('not sidechain')
    sc = residue.sidechain
    cg = residue.select('name CG')
    cb = residue.select('name CB')

    # Set charge and radii of side chain to 0, change CG atom to HB1
    residue.setResnames('ALA')
    sc.setCharges(0)
    sc.setRadii(0)
    cb.setRadii(2)
    cg.setNames('HB1')

    # Shorten the HB1-CB bond
    pos_hb1 = (0.7105 * (cg.getCoords() - cb.getCoords())) + cb.getCoords()
    cg.setCoords(pos_hb1)

    # Compile mutated pdb
    ala_atoms = ['N', 'H', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O', 'OXT']
    if chain is None:
        if preceed is None:
            mutant = residue.select('name ' + ' '.join(ala_atoms)) + follow
        if follow is None:
            mutant = preceed + residue.select('name ' + ' '.join(ala_atoms))
        if (preceed is None) and (follow is None):
            mutant = residue.select('name ' + ' '.join(ala_atoms))
        if (preceed is not None) and (follow is not None):
            mutant = preceed + residue.select('name ' + ' '.join(ala_atoms)) + follow
    else:
        if otherchains is None:
            if preceed is None:
                mutant = residue.select('name ' + ' '.join(ala_atoms)) + follow
            if follow is None:
                mutant = preceed + residue.select('name ' + ' '.join(ala_atoms))
            if (preceed is None) and (follow is None):
                mutant = residue.select('name ' + ' '.join(ala_atoms))
            if (preceed is not None) and (follow is not None):
                mutant = preceed + residue.select('name ' + ' '.join(ala_atoms)) + follow
        if otherchains is not None:
            if preceed is None:
                mutant = residue.select('name ' + ' '.join(ala_atoms)) + follow + otherchains
            if follow is None:
                mutant = preceed + residue.select('name ' + ' '.join(ala_atoms)) + otherchains
            if (preceed is None) and (follow is None):
                mutant = residue.select('name ' + ' '.join(ala_atoms)) + otherchains
            if (preceed is not None) and (follow is not None):
                mutant = preceed + residue.select('name ' + ' '.join(ala_atoms)) + follow + otherchains

    # Write mutant pqr
    pd.writePQR(mutid + '.pqr', mutant)


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

    if chain is None:
        for num, chid in zip(resnum, chain):
            print num, chid
            sel = selection(mdl.residue_range(int(num) - 1, int(num) - 1))
            sel.mutate(residue_type=resid)
    else:
        for num, chid in zip(resnum, chain):
            print num, chid
            if ' ' is chid:
                sel = selection(mdl.residue_range(int(num) - 1, int(num) - 1))
            else:
                sel = selection(mdl.residue_range(str(num) + ':' + chid, str(num) + ':' + chid))
            sel.mutate(residue_type=resid)

    aln.append_model(mdl, align_codes='mutant')
    mdl.clear_topology()
    mdl.generate_topology(aln['mutant'])
    mdl.transfer_xyz(aln)

    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.write(file=mutid + '.pdb')

    h = model(env, file=mutid + '.pdb')  # Without this section, chainids and resnums from parent won't be retained!
    m = model(env, file=pdb)
    aln = alignment(env)
    aln.append_model(m, atom_files=pdb, align_codes='parent')
    aln.append_model(h, atom_files=mutid + '.pdb', align_codes='mutant')
    h.res_num_from(m, aln)  # Restore old residue numbering and chain indexing
    h.write(file=mutid + '.pdb')

    ## Working Block
    # env = environ()
    # env.libs.topology.read(file='$(LIB)/top_heav.lib')
    # env.libs.parameters.read(file='$(LIB)/par.lib')
    #
    # aln = alignment(env)
    # mdl = model(env, file=pdb)
    # aln.append_model(mdl, atom_files=pdb, align_codes='parent')
    #
    # if ((chain is None) or (chain.isspace())):
    #     sel = selection(mdl.residue_range(int(resnum) - 1, int(resnum) - 1))
    # else:
    #     sel = selection(mdl.residue_range(str(resnum) + ':' + chain, str(resnum) + ':' + chain))
    #
    # sel.mutate(residue_type=resid)
    #
    # aln.append_model(mdl, align_codes='mutant')
    # mdl.clear_topology()
    # mdl.generate_topology(aln['mutant'])
    # mdl.transfer_xyz(aln)
    #
    # mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    # mdl.write(file=mutid + '.pdb')
    #
    # h = model(env, file=mutid + '.pdb')  # Without this section, chainids and resnums from parent won't be retained!
    # m = model(env, file=pdb)
    # aln = alignment(env)
    # aln.append_model(m, atom_files=pdb, align_codes='parent')
    # aln.append_model(h, atom_files=mutid + '.pdb', align_codes='mutant')
    # h.res_num_from(m, aln)  # Restore old residue numbering and chain indexing
    # h.write(file=mutid + '.pdb')


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
        outfile = os.path.splitext(pdbfile)[0] + '.pqr'
    # os.system('"{0}" {1} {2} {3}'.format(path_pdb2pqr_exe, optargs, pdbfile, outfile))
    (log, err) = runProcess([path_pdb2pqr_exe, '--ff=%s' % (ff), '--chain', pdbfile, outfile])
    return (log, err)


######################################################################################################################################################
# Function to run APBS.exe - should work on any supported OS
######################################################################################################################################################
def execAPBS(path_apbs_exe, pqr_chain, pqr_complex, prefix=None, grid=1.0, ion=0.150, pdie=20.0, sdie=78.54, cfac=1.5, dx=False):
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

    # cfac = 1.5 # hard-coded scaling factor for mesh dimension, for now

    pqr = pd.parsePQR(pqr_complex)
    coords = pqr.getCoords()
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    # Determine mesh dimensions according to Ron's AESOP protocol in the R source file
    fg = np.array((np.ceil(np.max(x) - np.min(x)), np.ceil(np.max(y) - np.min(y)), np.ceil(np.max(z) - np.min(z))))
    fg = np.ceil((fg + 5) * cfac)
    dime_list = (32 * np.linspace(1, 100, 100)) + 1  # list of possible dime values
    dime_ind = np.ceil(
        fg / (32 * grid)) - 1  # index of dime to use from list, subtract one to be consistent with python indexing!

    glen = fg
    dime = np.array((dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])], dime_list[int(dime_ind[2])]))

    # Format APBS input file
    cmd_read = ['read\n',
                '   mol pqr %s\n' % (pqr_chain),
                '   mol pqr %s\n' % (pqr_complex),
                'end\n']
    cmd_solv = ['elec name solv\n',
                '   mg-manual\n',
                '   dime %d %d %d\n' % (dime[0], dime[1], dime[2]),
                '   glen %d %d %d\n' % (glen[0], glen[1], glen[2]),
                '   gcent mol 2\n',
                '   mol 1\n',
                '   lpbe\n',
                '   bcfl sdh\n',
                '   srfm smol\n',
                '   chgm spl2\n',
                '   ion 1 %.2f 2.0\n' % (ion),
                '   ion -1 %.2f 2.0\n' % (ion),
                '   pdie %.2f\n' % (pdie),
                '   sdie %.2f\n' % (sdie),
                '   sdens 10.0\n',
                '   srad 0.0\n',
                '   swin 0.3\n',
                '   temp 298.15\n',
                '   calcenergy total\n']
    if dx is True:
        cmd_solv = cmd_solv + ['   write pot dx %s\n' % (prefix)]
    cmd_solv = cmd_solv + ['end\n']
    cmd_ref = ['elec name ref\n',
               '   mg-manual\n',
               '   dime %d %d %d\n' % (dime[0], dime[1], dime[2]),
               '   glen %d %d %d\n' % (glen[0], glen[1], glen[2]),
               '   gcent mol 2\n',
               '   mol 1\n',
               '   lpbe\n',
               '   bcfl sdh\n',
               '   srfm smol\n',
               '   chgm spl2\n',
               '   pdie %.2f\n' % (pdie),
               '   sdie %.2f\n' % (pdie),
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
    file_apbs_in = prefix + '.in'
    file_apbs_log = prefix + '.log'
    with open(file_apbs_in, 'w') as f:
        for line in apbs_in:
            f.write(line)

    # Execute APBS
    # os.system('"{0}" {1} {2}'.format(path_apbs_exe, '--output-file=%s --output-format=flat'%(file_apbs_log), file_apbs_in))
    # os.system('{0} {1}'.format(path_apbs_exe, file_apbs_in))
    # (log, err) = runProcess([path_apbs_exe, file_apbs_in])
    (log, err) = runProcess([path_apbs_exe, '--output-file=%s' % (file_apbs_log), '--output-format=flat', file_apbs_in])
    pattern = re.compile('(?<=Global net ELEC energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?')
    elec = np.asarray([x.split() for x in re.findall(pattern, log)]).astype(np.float)
    elec = elec.reshape((1, elec.size))

    # return file_apbs_log
    return elec

######################################################################################################################################################
# Function to run multiple APBS processes at once
######################################################################################################################################################
def batchAPBS(kernel):
    path, pqr_chain, pqr_complex, prefix, grid, ion, pdie, sdie, cfac, dx, i, j = kernel
    # print 'Calculating solvation and reference energies for: %s' % (os.path.basename(pqr_chain).split('.')[0])
    energies = execAPBS(path, pqr_chain, pqr_complex, prefix=prefix, grid=grid, ion=ion, pdie=pdie, sdie=sdie,
                        cfac=cfac, dx=dx)
    return np.array([i, j, energies[0][0], energies[0][1]])

######################################################################################################################################################
# Function to run multiple Coulomb processes at once
######################################################################################################################################################
def batchCoulomb(kernel):
    path, pqr_chain, pdie, i, j = kernel
    # print 'Calculating coulombic energies for: %s' % (os.path.basename(pqr_chain).split('.')[0])
    energies = execCoulomb(path, pqr_chain)
    energies = energies / pdie
    return np.array([i, j, energies])

######################################################################################################################################################
# Function to run coulomb.exe - should work on any supported OS
######################################################################################################################################################
def execCoulomb(path_coulomb_exe, pqr):
    (log, err) = runProcess([path_coulomb_exe, pqr])
    pattern = re.compile(
        '(?<=Total energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?')  # May need to update regex
    coul = np.asarray(re.findall(pattern, log)).astype(np.float)
    return coul

######################################################################################################################################################
# Function to run DSSP.exe - should work on any supported OS
######################################################################################################################################################
def execDSSP(pdbfile, dssp):
    (log, err) = runProcess([dssp, pdbfile])
    return log

######################################################################################################################################################
# Function to plot results of Alascan
######################################################################################################################################################
def plotResults(Alascan, filename=None):
    plt.style.use('seaborn-talk')
    figure, axarr = plt.subplots(len(Alascan.mutid) - 1, sharey=True)
    dpi_val = 300
    if len(Alascan.mutid) > 2:
        for i in xrange(1, len(Alascan.mutid)):
            axarr[i - 1].set_title(np.unique(np.array([w.split('_') for w in Alascan.mutid[i]])[:, 0])[0] + ' ddGbind relative to WT')
            axarr[i - 1].set_ylabel('kJ/mol')
            axarr[i - 1].set_xticks(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]])))
            if 100 < len(Alascan.mutid[i]) <= 150:
                axarr[i - 1].set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[i]])[:, 1], rotation='vertical', ha='left', size=6)
            elif len(Alascan.mutid[i]) > 150:
                axarr[i - 1].set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[i]])[:, 1], rotation='vertical', ha='left', size=2)
                dpi_val = 600
            else:
                axarr[i - 1].set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[i]])[:, 1], rotation='vertical', ha='left')
            axarr[i - 1].bar(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]]))[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] > 0],Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] > 0], color='red')
            axarr[i - 1].bar(np.arange(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]]))[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] < 0],Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] < 0], color='blue')
            axarr[i - 1].xaxis.set_ticks_position('bottom')
            axarr[i - 1].yaxis.set_ticks_position('left')
    elif len(Alascan.mutid) == 2:
        axarr.set_title(np.unique(np.array([w.split('_') for w in Alascan.mutid[1]])[:, 0])[0] + ' dGsolv relative to WT')
        axarr.set_ylabel('kJ/mol')
        axarr.set_xticks(np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]])))
        if 100 < len(Alascan.mutid[1]) <= 150:
            axarr.set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1], rotation='vertical', ha='left', size=6)
        elif len(Alascan.mutid[1]) > 150:
            axarr.set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1], rotation='vertical', ha='left', size=2)
            dpi_val = 600
        else:
            axarr.set_xticklabels(np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1], rotation='vertical', ha='left')
        axarr.bar(np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))[Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] > 0],Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] > 0], color='red')
        axarr.bar(np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))[Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] < 0],Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] < 0], color='blue')
        axarr.xaxis.set_ticks_position('bottom')
        axarr.yaxis.set_ticks_position('left')
    plt.tight_layout()
    if filename is not None:
        figure.savefig(filename, dpi=dpi_val)
    return(figure, axarr)

def plotResultsPlotly(Alascan, filename=None):
    """Summary
    
    Parameters
    ----------
    Alascan : TYPE
        Description
    filename : TYPE, optional
        Description
    
    Returns
    -------
    name : TYPE
        Description
    """
    subplot_titles = []
    for i in xrange(1,len(Alascan.mutid)):
         subplot_titles.append(np.unique(np.array([w.split('_') for w in Alascan.mutid[i]])[:,0])[0]+' ddGbind relative to WT')
    fig = tools.make_subplots(rows=len(Alascan.mutid) - 1, cols=1, vertical_spacing=.5, subplot_titles=subplot_titles)
    for i in xrange(1,len(Alascan.mutid)):
        pos_y = np.zeros(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]))
        pos_y[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] > 0] = Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] > 0]
        neg_y = np.zeros(len(Alascan.ddGbind_rel()[Alascan.mask_by_sel[:,i]]))
        neg_y[Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] < 0] = Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]][Alascan.ddGbind_rel()[Alascan.mask_by_sel[:, i]] < 0]
        pos_trace= go.Bar(
            x=np.array([w.split('_') for w in Alascan.mutid[i]])[:,1],
            y=pos_y,
            name = np.unique(np.array([w.split('_') for w in Alascan.mutid[i]])[:,0])[0]+'Loss of binding',
            marker = dict(
                color='rgba(0,136,55,1)'
            )
        )
        neg_trace = go.Bar(
            x=np.array([w.split('_') for w in Alascan.mutid[i]])[:,1],
            y=neg_y,
            name = np.unique(np.array([w.split('_') for w in Alascan.mutid[i]])[:,0])[0]+'Gain in binding',
            marker = dict(
                color='rgba(123,50,148,1)'
            )
        )
        fig.append_trace(pos_trace, i, 1)
        fig.append_trace(neg_trace, i, 1)
        fig['layout']['yaxis'+str(i)].update(title='kJ/mol')

    fig['layout'].update(barmode='stack', hovermode='closest')
    plotly_fig = go.Figure(fig)
    plotly.offline.plot(plotly_fig)
    if filename is not None:
        py.image.save_as(plotly_fig, filename=filename)

######################################################################################################################################################
# Function to plot results of ESD.calc()
######################################################################################################################################################
def plotESD(esd, filename=None, cmap='hot'):
    plt.style.use('seaborn-talk')
    fig, ax = plt.subplots(sharey=True)
    heatmap = ax.pcolor(esd.esd, cmap=cmap, vmin=0, vmax=1)
    ax.set_xlim(0,esd.esd.shape[0])
    ax.set_ylim(0,esd.esd.shape[1])
    ax.set_xticks(np.arange(esd.esd.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(esd.esd.shape[1])+0.5, minor=False)
    ax.set_xticklabels(esd.ids, rotation=90 )
    ax.set_yticklabels(esd.ids)
    fig.colorbar(heatmap)
    plt.tight_layout()
    if filename is not None:
        fig.savefig(filename)

######################################################################################################################################################
# Function to plot ESD dendrogram
######################################################################################################################################################
def plotDend(esd, filename=None):
    plt.style.use('seaborn-talk')
    fig, ax = plt.subplots(sharey=True)
    Z = cluster.linkage(esd.esd, 'ward')
    cluster.dendrogram(
        Z,
        labels=esd.ids,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        ax=ax
    )
    plt.xlabel('Variants')
    plt.ylabel('ESD')
    ax.set_xticklabels(esd.ids, rotation=90 )
    # ax.set_ylim(0,1)
    plt.tight_layout()
    if filename is not None:
        fig.savefig(filename)

def plotESDPlotly(esd, filename=None, cmap='YIGnBu'):
    figure = FF.create_dendrogram(esd.esd, orientation='bottom', labels=esd.ids)
    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'
    # Create Side Dendrogram
    dendro_side = FF.create_dendrogram(esd.esd, orientation='right', labels=esd.ids)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    # Add Side Dendrogram Data to Figure
    figure['data'].extend(dendro_side['data'])
    # Create Heatmap
    dendro_side2 = FF.create_dendrogram(esd.esd, orientation='right')
    for i in range(len(dendro_side2['data'])):
        dendro_side2['data'][i]['xaxis'] = 'x2'
    dendro_leaves = dendro_side2['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))

    heat_data = esd.esd
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]

    heatmap = go.Data([
        go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=heat_data,
            colorscale=cmap
        )
    ])
    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = figure['layout']['xaxis']['tickvals']

    # Add Heatmap Data to Figure
    figure['data'].extend(go.Data(heatmap))

    # Edit Layout
    figure['layout'].update({
                             'showlegend':False, 'hovermode': 'closest',
                             })
    figure['layout'].update({'margin':{'b':140,
                                      't':10}})

    # Edit xaxis
    figure['layout']['xaxis'].update({'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks': ""})
    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks': ""}})

    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'side': 'right',
                                      'ticktext': dendro_side['layout']['yaxis']['ticktext'],
                                      'tickvals': dendro_side['layout']['yaxis']['tickvals'],
                                      'ticks': ""})
    # Edit yaxis2
    figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks': ""}})

    # Plot!
    plotly.offline.plot(figure)
    if filename is not None:
        py.image.save_as(figure, filename=filename)

######################################################################################################################################################
# Function to calculate RSA for PDB
######################################################################################################################################################
def calcRSA(pdbfile, dssp):
    # SA from: C. Chotia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14.
    SA_dict = {'C': 135, 'D': 150, 'S': 115, 'Q': 180, 'K': 200,
               'I': 175, 'P': 145, 'T': 140, 'F': 210, 'N': 160,
               'G': 75, 'H': 195, 'L': 170, 'R': 225, 'W': 255,
               'A': 115, 'V': 155, 'E': 190, 'Y': 230, 'M': 185}
    threshold = 0.2 # RSA below this value is considered completely buried
    log = execDSSP(pdbfile, dssp)
    ag = pd.parsePDB(pdbfile)
    n_atoms = ag.numAtoms()
    ACC = np.zeros(n_atoms, float)
    lines = filter(None, log.split('\n'))
    iterator = iter(lines)
    for line in iterator:
        # print line
        if line.startswith('  #  RESIDUE'):
            break
    for line in iterator:
        if line[13] == '!':
            continue
        res = ag[(line[11], int(line[5:10]), line[10].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        ACC[indices] = int(line[35:38])
        ag.setData('dssp_acc', ACC)

    resid = np.asarray([AA_dict[x] for x in ag.getResnames()])
    sasa = ag._getData('dssp_acc')
    rsa =[]
    for res, sa in zip(resid, sasa):
        rsa.append(sa / SA_dict[res])
        # rsa = np.asarray([SA_dict[res] for res, sasa in zip(resid, sasa])
    rsa = np.asarray(rsa)
    return(resid, rsa, rsa>=threshold)

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
           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
