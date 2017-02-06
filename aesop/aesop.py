import os as os
import sys as sys
import glob as glob
import subprocess as sp
import datetime as dt
import timeit as ti
import re as re
import numpy as np
import prody as pd
import scipy.cluster.hierarchy as cluster
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Pool, cpu_count  # , freeze_support
import itertools as it

# Print Licencse on startup
print("""
AESOP: Analysis of Electrostatic Structure of Proteins

Reed E. S. Harrison, Rohith R. Mohan, Dimitrios Morikis
University of California, Riverside; Department of Bioengineering

Correspondence should be directed to Prof. Dimitrios Morikis at
dmorikis@ucr.edu

Copyright (C) 2016  Reed E. S. Harrison, Rohith R. Mohan, Ronald D.
Gorham Jr., Chris A. Kieslich, Dimitrios Morikis

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
ERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
""")


class Grid:
    """Summary

    The grid class facilitates parsing and writing of OpenDX file formats.
    In the current state, the class is quite rudimentary and only supports
    changing vectors for the grid data.

    Attributes
    ----------

    filename : string
        DX file to import
    pot : ndarray
        Vectors at each grid point. For AESOP, these will typically be
        electrostatic potentials or an electrostatic similarity index.
    header : list
        List of grid parameters from the OpenDX format prior to vectors.
    footer : list
        List of grid parameters from the OpenDX format subsequent to the
        vectors.
    """
    import re as re

    def __init__(self, filename=None):
        self.filename = filename
        if filename is not None:
            self.readDX()

    def readDX(self, filename=None):
        """Summary

        Method to parse a DX file

        Parameters
        ----------

        filename : string
            Name for the OpenDX file to be imported. If unspecified, this
            parameter defaults to the class attribute.
        """
        lines = []
        if filename is None:
            try:
                filename = self.filename
                f = open(filename, 'r')
                lines = f.readlines()
                f.close()
            except:
                print '\nError: No DX file provided\n'
        else:
            try:
                f = open(filename, 'r')
                lines = f.readlines()
                f.close()
            except:
                print '\nError: Can\'t open DX file\n'

        # f       = open(fn, 'r')
        # lines   = f.readlines()
        p_vec = '[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
        # p_count
        # p_delta
        # p_rank
        # p_g

        # comment = []
        # obj     = []
        # origin  = []
        # delta   = []
        # attr    = []
        # comp    = []
        # pot     = []

        header = []
        footer = []
        pot = []

        found_pot = False
        for line in lines:
            if len(re.findall('^#', line)) != 0:
                # comment.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall('^object', line)) != 0:
                # obj.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall('^origin', line)) != 0:
                # origin.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall('^delta', line)) != 0:
                # delta.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall('^attribute', line)) != 0:
                # attr.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall('^component', line)) != 0:
                # comp.append(line)
                if found_pot:
                    footer.append(line)
                elif not found_pot:
                    header.append(line)
                continue
            if len(re.findall(p_vec, line)) != 0:
                found_pot = True
                x = re.findall(p_vec, line)
                pot.append(np.asarray(x).astype(float))
                continue
            # else:
            #     if found_pot:
            #         footer.append(line)
            #     elif not found_pot:
            #         header.append(line)
            # f.close()
            # return np.vstack(pot)

        self.pot = np.vstack(pot)
        self.header = header
        self.footer = footer

    def writeDX(self, filename=None):
        """Summary

        Function to write OpenDX files

        Parameters
        ----------
        filename : string
            Name for OpenDX file that will be written. This should be a full
            path if you wish to place the file somewhere other than the
            current working directory.
        """
        if filename is None:
            filename = os.path.splitext(
                os.path.basename(self.filename))[0] + '.modified.dx'

        header = self.header
        try:
            pot = [
                ' '.join([format(x, '7.6e') for x in v]) + '\n'
                for v in self.pot.tolist()
            ]
        except:
            pot = [format(v, '7.6e') + '\n' for v in self.pot.tolist()]
        footer = self.footer

        f = open(filename, 'w')
        f.writelines(header + pot + footer)
        f.close()


##########################################################################
# Container for performing an Alanine Scan with AESOP
#   pdb     -   PDB file for performing Alascan.
#               Must contain all chain selections with
#               standard aminoacid nomenclature
#   selstr  -   List of selection strings for performing mutations
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant
#   sdie    -   Solvent dielectric constant
##########################################################################


class Alascan:
    """Summary
    Summary of internal varialbles in the Alascan class.
    All parameters are set in the Alascan.__init__(...)
    and the analysis is run with Alascan.run().

    Attributes
    ----------
    apbs : str
        Full path to APBS executable. Must be compatible with OS.
    apbs_results : str
        Full path to output from APBS
    cfac : float
        Scaling factor for grid dimensions. We suggest to leave this
        unchanged.
    coulomb : str
        Full path to coulomb executable from APBS package.
        Must be compatible with OS.
    coulomb_results : str
        Full path to folder containing results from coulomb.
    dime : list
        List of three integers. Parameter required for APBS.
        Please see description at:
        http://www.poissonboltzmann.org/docs/apbs-overview/
    disu : bool, optional
        If True, Modeller will guess the patches for disulfide bridges
        within the provided protein structures. Only relevant if minim
        is set to Trueself.
    dx : bool
        Variable that specifies if potential files should be written to disk.
        Default is False.
    dx_files : list
        If written to disk, list of potential files written to disk.
        The folder containing these files is given by Alascan.apbs_results.
    E_ref : ndarray
        Description
    E_solv : ndarray
        Description
    ff : str
        Force-field to use for PDB2PQR
    file_pdb_template : TYPE
        Description
    gcent : list
        List of three integers. Parameter required for APBS.
        Please see description at:
        http://www.poissonboltzmann.org/docs/apbs-overview/
    Gcoul : ndarray
        Coulombic free energies, corresponding to Alascan.mutid.
    glen : list
        List of three integers. Parameter required for APBS.
        Please see description at:
        http://www.poissonboltzmann.org/docs/apbs-overview/
    Gref : ndarray
        Reference-state free energies, corresponding to Alascan.mutid.
    grid : numeric
        Distance spacing of grid points.
        If the grid dimensions are not divisible by three,
        resolution will be increased (smaller grid spacing) until
        grid dimensions are divisible by three.
    Gsolv : ndarray
        Solvated-state free energies, corresponding to Alascan.mutid.
    ion : numeric
        Ionic strength to be used in the solvated-state APBS calculations.
    jobdir : str
        [Optional] Path to folder containing results.
        If not specified, a directory will be generated.
    jobname : str
        [Optional] Name for current job, will be used to create the jobdir.
    list_chids : list
        Chain ID where mutation was made. Corresponds to Alascan.mutid.
    list_resnames : list
        Residue names where mutation was made. Corresponds to Alascan.mutid.
    list_resnums : list
        Residue numbers where mutation was made. Corresponds to Alascan.mutic.
    logs : list
        List of strings that represent the log files from every executable
        called (namely, PDB2PQR and APBS)
    logs_apbs_dir : str
        Folder in jobdir containing output from
        APBS (logs, input files, dx files)
    mask_by_sel : ndarray
        Matrix containing selection masks. The first column corresponds
        to the selection string for the parent and each column
        thereafter corresponds to an element of the selection
        string (selstr) in the same order.
    max_iter : integer, optional
        If minimization is to be performed, this parameter limits the
        maximum number of calls to the objective function before
        minimization is terminated. Default is 1000 iterations.
    min_atom_shift : float, optional
        If minimization is to be performed, this parameter will determine
        the convergence criteria. If the maximum atomic shift for all atoms
        between minimization steps is less than this value, then
        minimization is terminated. Default value is 0.1 angstroms.
    minim : bool, optional
        If True, minimization will be performed with Modeller's conjugate
        gradient descent algorithm. Default is False for the Alanine scan
        class as no clashes should result from mutations.
    mutid : list
        List of mutant IDs. The first element corresponds to the parent.
        Subsequent elements correspond to each element of the
        selection string list (selstr). Please use Alascan.getMutids()
        to get vectorized version.
    output : string, optional
        If minimization is performed, this parameter deterimines what
        output to STDOUT Modeller will use. 'NO_REPORT' results in a
        minimal output to STDOUT, while 'REPORT' results in a more
        verbose output to STDOUT.
    pdb : str
        Path to PDB file with atomic coordinates. Must follow formatting
        conventions of the Protein Databank.
    pdb2pqr : str
        Full path to PDB2PQR executable.
    pdb_complex_dir : str
        Folder name in the job directory that contains the PDB file(s) of
        the complex structures.
    pdie : numeric
        Protein dielectric constant to be used in APBS.
    pqr_complex_dir : str
        Folder name in the job directory containing PQR files for each
        protein complex.
    pqr_sel_dir : list
        List of folder names in the job directory that contain PQR files for
        selection strings (selstr). Each element of pqr_sel_dir corresponds
        to the same element of selstr.
    region : list
        List of additional selection strings specifying the zone where
        mutations should occur. Generally unused, unless a region of interest
        is involved. Each element of region should correspond to the same
        element of selstr. That is, each region selection string will futher
        narrow down the initial selection string.
    sdie : numeric
        Solvent dielectric constant to be used in APBS.
    selstr : list
        List of selection strings. Typically each selection string will
        correspond to a chain in a protein complex. We advise users to
        specify two or more selection strings. (Ex: ['chain A', 'chain B'])

    Deleted Attributes
    ------------------
    dirs : int
        Description
    prefix : TYPE
        Description
    """

    def __init__(self,
                 pdb,
                 pdb2pqr_exe='pdb2pqr',
                 apbs_exe='apbs',
                 coulomb_exe='coulomb',
                 selstr=['protein'],
                 jobname=None,
                 region=None,
                 grid=1,
                 ion=0.150,
                 pdie=20.0,
                 sdie=78.54,
                 ff='parse',
                 cfac=1.5,
                 dx=False,
                 minim=False):
        """Summary
        Constructor for the Alascan class.

        Parameters
        ----------
        pdb : str
            The pdb file for which a computational alanine scan will be
            performed.
        pdb2pqr_exe : str
            Path to executable for PDB2PQR.
        apbs_exe : str
            Path to executable for APBS
        coulomb_exe : str, optional
            Path to executable for Coulomb from the APBS toolbox.
        selstr : list, optional
            List of selection strings. Typically each selection string will
            correspond to a chain in a protein complex. We advise users to
            specify two or more selection strings. The default selection will
            only investigate solvation free energies as there will be no
            binding event. (Ex: ['chain A', 'chain B'])
        jobname : str, optional
            A string representing the full path to the directory where
            results will be stored.
        region : str, optional
            If provided, region should be a selection string that specifies
            where mutations should be made. This must be the same length as
            selstr. Elements of selstr and region must correspond.
        grid : int, optional
            Desired grid spacing to be used in APBS. Depending on the PDB
            file, actual grid spacing may be a finer resolution.
            Defaults to 1 angstrom.
        ion : float, optional
            Ionic strength (M) to be used in APBS.
        pdie : float, optional
            Protein dielectric constant to be used in APBS.
        sdie : float, optional
            Solvent dielectric constant to be used in APBS>
        ff : str, optional
            Force field to assign charges to individual atoms.
            We recommend using PARSE for this analysis, the default value.
        cfac : float, optional
            Scaling factor for grid dimension.
            This generally should not be changed without a reason.
        dx : bool, optional
            Boolean flag stating if potential files should be written to
            disk (True) or not (False).
        minim : bool, optional
            If TRUE, energy minimization of protein structures will be
            performed with Modeller's conjugate gradient descent
            algorithm. Default is False for the Alanine Scan class as
            mutations should not result in atomic clashes.
        """
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
            jobname = '%4d%02d%02d_%02d%02d%02d' % (
                dt.date.today().year, dt.date.today().month,
                dt.date.today().day, dt.datetime.now().hour,
                dt.datetime.now().minute, dt.datetime.now().second)
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
        self.status = 0  # 0 = alascan not run, 1 = alascan previously ran

        self.genDirs()
        self.genMutid()
        self.genParent()
        self.find_grid()

        self.minim = minim
        self.disu = True
        self.min_atom_shift = 0.1
        self.max_iter = 1000
        # self.report_iter = 10
        self.output = 'NO_REPORT'
        self.logs = []

    # def getPDB(self):
    #     """Summary

    #     Returns
    #     -------
    #     string
    #         Description
    #     """
    #     return self.pdb

    # def getSel(self):
    #     """Summary

    #     Returns
    #     -------
    #     TYPE
    #         Description
    #     """
    #     return self.selstr

    # def getDirs(self):
    #     """Summary

    #     Returns
    #     -------
    #     TYPE
    #         Description
    #     """
    #     return self.dirs

    # def getPrefix(self):
    #     """Summary

    #     Returns
    #     -------
    #     TYPE
    #         Description
    #     """
    #     return self.prefix

    def getMutids(self):
        """Summary

        Returns
        -------
        list
            Returns vectorized format of mutids.
        """
        l = self.mutid
        return [item for sublist in l for item in sublist]

    # def getDX(self):
    #     """Summary

    #     Returns
    #     -------
    #     TYPE
    #         Description
    #     """
    #     return self.dx_files

    def genDirs(self):
        """Summary
        This subroutine will generate all directories needed to contain
        structural files, logs, etc. In the future we may implement a
        method to remove such files when outputs are more standardized.

        Returns
        -------
        None
            No output, operates on the class and generates folders in the job
            directory.
        """
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
            if not os.path.exists(
                    os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1))):
                os.makedirs(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1)))
        self.pqr_complex_dir = pqr_complex_dir
        self.pqr_sel_dir = pqr_sel_dir

        # Create necessary directories for APBS files
        logs_apbs_dir = 'apbs_logs'
        if not os.path.exists(os.path.join(self.jobdir, logs_apbs_dir)):
            os.makedirs(os.path.join(self.jobdir, logs_apbs_dir))
        self.logs_apbs_dir = 'apbs_logs'

    def genMutid(self):
        """Summary
        This subroutine reads the input PDB, selects the structure where
        mutations will occur, and saves all mutant IDs in the class.
        If region is specified in the constructor, then the constraint
        will be applied here.

        Returns
        -------
        None
            Operates on the class to generate a list of mutant IDs for each
            selection string. The following class variables will be generated
            (see class description):
                - mutid
                - list_chids
                - list_resnums
                - list_resnames
                - mask_by_sel
        """
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
            combined_selection = parent_pdb.select(' and '.join(
                ['(' + sel + ')', '(' + reg + ')', 'charged', 'calpha']))
            list_chids[i] = combined_selection.getChids().tolist()
            list_resnums[i] = combined_selection.getResnums().tolist()
            list_resnames[i] = combined_selection.getResnames().tolist()
            code = [
                'seg%d' % (i) + '_' + AA_dict[res_id] + res_no + 'A'
                for ch_id, res_no, res_id in zip(
                    list_chids[i], map(str, list_resnums[i]), list_resnames[i])
            ]
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
        """Summary
        Reads PDB file specified in the constructor; applies and combines
        results from the selection strings; and saves the final template
        structure in the job directory.

        Returns
        -------
        None
            Template pdb written in job directory and location saved in
            Alascan.file_pdb_template.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        infile = os.path.join(jobdir, pdb_complex_dir,
                              parent_file_prefix + '.pdb')
        system = parent_pdb.select('((' + ') or ('.join(selstr) + '))')
        # print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)
        self.file_pdb_template = infile

    def find_grid(self):
        """Summary
        Calculate grid dimensions for APBS (dime, glen, gcent)

        Returns
        -------
        TYPE
            Sets class attributes dime, glen, and gcent
        """
        pdbfile = self.file_pdb_template
        grid = self.grid
        cfac = self.cfac

        glen = np.zeros(3)

        # Find bounding box
        pdb = pd.parsePDB(pdbfile)
        max_xyz = pdb.getCoords().max(axis=0)
        min_xyz = pdb.getCoords().min(axis=0)

        # Determine mesh dimensions according to Ron's AESOP protocol in the R
        # source file
        fg = np.ceil(max_xyz - min_xyz)
        fg = np.ceil((fg + 5) * cfac)
        glen = np.vstack((glen, fg)).max(axis=0)
        dime_list = (32 * np.linspace(1, 100, 100)) + \
            1  # list of possible dime values
        # index of dime to use from list, subtract one to be consistent with
        # python indexing!
        dime_ind = np.ceil(glen / (32 * grid)) - 1
        dime = np.array(
            (dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])],
             dime_list[int(dime_ind[2])]))
        ix = 0
        iy = 0
        iz = 0
        counter = 0
        while ((dime[0] * dime[1] * dime[2] % 3 != 0) or (counter <= 5)):
            ix += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iy += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iz += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            counter += 1

        self.dime = dime
        self.glen = glen
        self.gcent = pd.calcCenter(pdb).astype(int)

    def set_grid(self, dime, glen, gcent):
        """Summary
        In the case that the user wishes to manually specify grid dimension,
        this may be accomplished with the set_grid method.
        Typically, this is used when grid space parameters must be consistent
        for many analyses. Please see description
        at: http://www.poissonboltzmann.org/docs/apbs-overview/ for
        description of parameters.

        Parameters
        ----------
        dime : list
            List of three integers.
        glen : list
            List of three integers.
        gcent : list
            List of three integers.

        Returns
        -------
        TYPE
            Description
        """
        self.dime = dime
        self.glen = glen
        self.gcent = gcent

    def genTruncatedPQR(self):
        """Summary
        Generate a structure for each mutant ID by truncating the side chain
        to form alanine.

        Returns
        -------
        None
            Write library of mutant structures to disk for subsequent analysis.
        """

        def minimize_pqr(self, pqrfile):
            minimize_cg(
                pqrfile,
                dest=pqrfile,
                disu=self.disu,
                min_atom_shift=self.min_atom_shift,
                max_iter=self.max_iter,
                output=self.output)
            (pqr_log, pqr_errs) = execPDB2PQR(
                self.pdb2pqr, pqrfile, outfile=pqrfile, ff=self.ff)
            self.logs.append(pqr_log)

        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_pdb2pqr = self.pdb2pqr
        ff = self.ff

        list_mutids = [item for sublist in self.mutid for item in sublist]
        list_chids = [item for sublist in self.list_chids for item in sublist]
        list_resnums = [
            item for sublist in self.list_resnums for item in sublist
        ]

        # Calculate PQR for parent
        infile = os.path.join(jobdir, pdb_complex_dir, list_mutids[0] + '.pdb')
        outfile = os.path.join(jobdir, pqr_complex_dir,
                               list_mutids[0] + '.pqr')
        print '\n%s:\tgenerating PQR for parent: %s' % (self.jobname,
                                                        list_mutids[0])
        (pqr_log, pqr_errs) = execPDB2PQR(
            path_pdb2pqr, infile, outfile=outfile, ff=ff)
        self.logs.append(pqr_log)
        if self.minim == True:
            minimize_pqr(self, outfile)
        complex_pqr = pd.parsePQR(outfile)
        for sel, seldir in zip(selstr, pqr_sel_dir):
            selfile = os.path.join(jobdir, seldir, list_mutids[0] + '.pqr')
            pqr = complex_pqr.select(sel)
            pd.writePQR(selfile, pqr)
            if self.minim == True:
                minimize_pqr(self, selfile)

        for mutid, chain, resnum in zip(list_mutids[1:], list_chids[1:],
                                        list_resnums[1:]):
            outpath = os.path.join(jobdir, pqr_complex_dir, mutid)
            print '\n%s:\tgenerating PQR for mutant: %s' % (self.jobname,
                                                            mutid)
            mutatePQR(outfile, mutid=outpath, resnum=resnum, chain=chain)
            complex_pqr = pd.parsePQR(outpath + '.pqr')
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid + '.pqr')
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)
                if self.minim == True:
                    minimize_pqr(self, selfile)

    def calcAPBS(self):
        """Summary
        Run APBS on each structure in mutant library, in serial.

        Returns
        -------
        None
            Stores energies from APBS in Gref and Gsolv class attributes.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel))
        Gref = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating solvation and reference energies for ' \
                'mutant: %s' % (self.jobname, mutid)
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_prefix_log = os.path.join(jobdir, logs_apbs_dir, mutid)
                if mask_by_sel[i, j]:
                    energies, log = execAPBS(
                        path_apbs,
                        subunit_pqr,
                        self.dime,
                        self.glen,
                        self.gcent,
                        prefix=path_prefix_log,
                        ion=self.ion,
                        pdie=self.pdie,
                        sdie=self.sdie,
                        dx=self.dx)
                    self.logs.append(log)
                    Gsolv[i, j] = energies[0][0]
                    Gref[i, j] = energies[0][1]
                if not mask_by_sel[i, j]:
                    Gsolv[i, j] = Gsolv[0, j]
                    Gref[i, j] = Gref[0, j]
        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcAPBS_parallel(self, n_workers=None):
        """Summary
        Run APBS on each structure in mutant library, in parallel.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all available
            threads.

        Returns
        -------
        None
            Saves solvated-state and reference-state free energies as class
            attributes.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        if n_workers is None:
            n_workers = cpu_count()

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel)).astype(float)
        Gref = np.zeros((dim_mutid, dim_sel)).astype(float)

        path_list = []
        pqr_chain_list = []
        dime_list = []
        glen_list = []
        gcent_list = []
        prefix_list = []
        ion_list = []
        pdie_list = []
        sdie_list = []
        dx_list = []
        i_list = []
        j_list = []

        # Find all calculations to be done
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                if mask_by_sel[i, j]:
                    path_list.append(path_apbs)
                    pqr_chain_list.append(subunit_pqr)
                    dime_list.append(self.dime)
                    glen_list.append(self.glen)
                    gcent_list.append(self.gcent)
                    prefix_list.append(
                        os.path.join(jobdir, logs_apbs_dir,
                                     '%d_%d_' % (i, j) + mutid)
                    )  # added to make sure apbs.in file is unique!
                    ion_list.append(self.ion)
                    pdie_list.append(self.pdie)
                    sdie_list.append(self.sdie)
                    dx_list.append(self.dx)
                    i_list.append(i)
                    j_list.append(j)

        # Organize kernel and run batch process
        kernel = zip(path_list, pqr_chain_list, dime_list, glen_list,
                     gcent_list, prefix_list, ion_list, pdie_list, sdie_list,
                     dx_list, i_list, j_list)

        p = Pool()
        print '%s:\trunning batchAPBS ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchAPBS, kernel):
            counter += 1
            print '.... %s:\tbatch APBS %d percent complete ....' % (
                self.jobname, int(counter * 100 / max_count))
            data, log = result
            i = int(data[0])
            j = int(data[1])
            solv = data[2]
            ref = data[3]
            self.logs.append(log)
            Gsolv[i, j] = solv
            Gref[i, j] = ref
        if self.dx == True:
            self.dx_files = [x + '.dx' for x in prefix_list]

        # Fill in results that are duplicates
        for i in xrange(dim_mutid):
            for j in xrange(dim_sel):
                if not mask_by_sel[i, j]:
                    Gsolv[i, j] = Gsolv[0, j]
                    Gref[i, j] = Gref[0, j]

        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcCoulomb(self):
        """Summary
        Calculates Coulombic free energies with coulomb.exe from the APBS
        toolbox.

        Returns
        -------
        None
            Saves Coulombic free energies as a class attribute.
        """
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

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gcoul = np.zeros((dim_mutid, dim_sel))

        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating coulombic energies for mutant: %s' % (
                self.jobname, mutid)
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                if mask_by_sel[i, j]:
                    subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                    energies, log = execCoulomb(path_coulomb, subunit_pqr)
                    self.logs.append(log)
                    Gcoul[i, j] = energies / pdie

        # Fill in results that are duplicates
        for i in xrange(dim_mutid):
            for j in xrange(dim_sel):
                if not mask_by_sel[i, j]:
                    Gcoul[i, j] = Gcoul[0, j]

        self.Gcoul = Gcoul

    def calcCoulomb_parallel(self, n_workers=None):
        """Summary
        Calculates Coulombic free energies with coulomb.exe from the APBS
        toolbox in a parallel manner.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all available
            threads.

        Returns
        -------
        None
            Saves Coulombic free energies as a class attribute.
        """
        selstr = self.selstr
        region = self.region
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_coulomb = self.coulomb

        if n_workers is None:
            n_workers = cpu_count()

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
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
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
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
        p = Pool(n_workers)
        print '%s:\trunning batchCoulomb ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchCoulomb, kernel):
            counter += 1
            print '.... %s:\tbatch coulomb %d percent complete ....' % (
                self.jobname, int(counter * 100 / max_count))
            data, log = result
            i = int(data[0])
            j = int(data[1])
            coul = data[2]
            coulomb_results.append([i, j, coul])
            Gcoul[i, j] = coul
            self.logs.append(log)
        coulomb_results = np.asarray(coulomb_results)
        self.coulomb_results = coulomb_results

        # Fill in results that are duplicates
        for i in xrange(dim_mutid):
            for j in xrange(dim_sel):
                if not mask_by_sel[i, j]:
                    Gcoul[i, j] = Gcoul[0, j]
        self.Gcoul = Gcoul

    def ddGa_rel(self):
        """Summary
        Calculates and returns the free energy of association relative to the
        parent free energy of association.

        Returns
        -------
        ndarray
            Array of free energies corresponding to the mutant IDs from the
            Alascan.getMutIDs() method.
        """
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
        """Summary
        Calculates and returns the free energy of a solvation for each mutant
        relative to the parent free energy of solvation.

        Returns
        -------
        ndarray
            Array of solvation free energies corresponding to mutant IDs from
            the Alascan.getMutIDs() method.
        """
        Gsolv = self.Gsolv
        Gref = self.Gref
        dGsolv = Gsolv - Gref
        dGsolv = dGsolv - dGsolv[0, 0]
        return dGsolv[:, 0]

    def calcESI(self, idx=-1):
        """Summary

        Compare potential files and calculate the similarity index.
        Values closer to 1 imply similarity while values closer to zero imply
        dissimilarity.

        Parameters
        ----------
        method : str, optional
            This parameter will allow for other metrics to compare
            grid potentials; however, for now only 'AND' is implemented.
        idx : int
            Index of original PDB files supplied containing reference
            structure. Set to None to perform all pairwise comparisons.

        Returns
        -------
        None
            Writes esi files to the esi_files directory within the job
            directory.
        """

        # If dx = true, use previously calculated DX files, else

        safe = True
        if self.dx == False:
            print('Error: please set the dx argument to true when initializing'
                  ' the alanine scan class, then re-run the analysis')
            safe = False
        if self.status == 0:
            print('Error: please run the alanine scan before attempting to'
                  ' calculate ESIs')
            safe = False

        if safe:
            n_seg = len(self.selstr)
            ind_seg = [x + 1 for x in range(n_seg)]

            file_by_ref = []
            for ind in ind_seg:
                files = [
                    x for x in self.dx_files
                    if int(os.path.basename(x).split('_')[1]) == ind
                ]
                file_by_ref.append(files)

            ids = ['selection %d' % x for x in ind_seg]

            esidir = os.path.join(self.jobdir, 'esi_files')
            self.esidir = esidir
            if not os.path.exists(os.path.join(self.esidir)):
                os.makedirs(os.path.join(self.esidir))

            esifiles = []
            esilist = []
            for i in ind_seg:
                ref_name = ids[i - 1]
                ref = Grid(file_by_ref[i - 1][0])
                dim = ref.pot.size
                n = len(self.mutid[i])

                esi = []
                filename = os.path.join(esidir, ref_name + '.dx')
                for j in xrange(n):
                    # print file_by_ref[i-1][j+1]
                    dat = Grid(file_by_ref[i - 1][j + 1])
                    a = ref.pot.astype(float).reshape((dim, ))
                    b = dat.pot.astype(float).reshape((dim, ))

                    diff = np.abs(a - b)
                    maxpot = np.abs(np.vstack((a, b))).max(axis=0)
                    val = np.divide(diff, maxpot)

                    esi.append(val)

                esi = np.vstack(esi)
                esi = np.sum(esi, axis=0) / n
                ref.pot = esi.reshape((dim / 3, 3))
                ref.writeDX(filename)
                esifiles.append(filename)
                esilist.append(esi)
            esilist = np.vstack(esilist)
            self.esifiles = esifiles
            self.esi = esilist

    def run(self):
        """Summary
        Perform a compuational alanine scan on the initialized Alascan class.

        Returns
        -------
        None
            Outputs text to STDOUT when run is complete, will be made optional
            in the future.
        """
        start = ti.default_timer()
        self.logs = []
        self.genTruncatedPQR()
        self.calcAPBS()
        self.calcCoulomb()
        self.status = 1
        stop = ti.default_timer()
        print '%s:\tAESOP alanine scan completed in %.2f seconds' % (
            self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def run_parallel(self, n_workers=None):
        """Summary
        Perform a computational alanine scan on the initialized Alascan class
        using multiple processes in parallel.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all available
            threads.

        Returns
        -------
        None
            Outputs text to STDOUT when run is complete, will be made optional
            in the future.
        """
        start = ti.default_timer()
        self.logs = []
        self.genTruncatedPQR()
        self.calcAPBS_parallel(n_workers)
        self.calcCoulomb_parallel(n_workers)
        self.status = 1
        stop = ti.default_timer()
        print '%s:\tAESOP alanine scan completed in %.2f seconds' % (
            self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def summary(self, filename=None):
        """Summary
        Summarize results from the computational alanine scan once complete.

        Parameters
        ----------
        filename : str, optional
            In order to write summary to a text file, supply the filename
            (full path).

        Returns
        -------
        None
            Prints summary of residues and energies relative to the parent
            structure if no filename is provided. Otherwise, writes to text
            file.
        """
        selstr = self.selstr
        mutids = self.getMutids()
        energies = None
        if len(selstr) > 1:
            energies = self.ddGa_rel()
        elif len(selstr) == 1:
            energies = self.dGsolv_rel()
        lines = ['%s, %f' % (lbl, val) for lbl, val in zip(mutids, energies)]
        lines.append('\nKey\n---\n')
        for i, seg in enumerate(selstr):
            lines.append('sel%d = %s' % (int(i + 1), str(seg)))
        if filename is None:
            print(lines)
        if filename is not None:
            np.savetxt(filename, lines, fmt='%s')

    def writeLogs(self, filename=None):
        if filename is None:
            filename = os.path.join(self.jobdir, 'AESOP_logs.txt')
        with open(filename, 'w') as f:
            logs = '\n==== Log Instance ====\n'.join(self.logs)
            f.write(logs)

    def viewLogs(self):
        logs = '\n==== Log Instance ====\n'.join(self.logs)
        print logs

    def checkwarnings(self):
        lines = ''.join(self.logs)
        matches = re.findall('[W][Aa][Rr][Nn][Ii][Nn][Gg]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1

    def checkerrors(self):
        lines = ''.join(self.logs)
        matches = re.findall('[E][Rr][Rr][Oo][Rr]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1


##########################################################################
# Container for performing an Directed Mutagenesis Scan with AESOP
#   pdb     -   PDB file for performing Alascan. Must contain all chain
#               selections with standard amino acid nomenclature
#   selstr  -   List of selection strings for chain selection. This does not
#               change what is mutated.
#   target  -   List of selection strings for mutants, each will element of
#               list will be mutated simultaneously
#   mutation -  Mutate each target residue to the corresponding element of
#               this list. (3 letter code)
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant
#   sdie    -   Solvent dielectric constant
##########################################################################


class DirectedMutagenesis:
    """Summary

    Attributes
    ----------
    apbs : str
        Full path to APBS executable. Must be compatible with OS.
    apbs_results : str
        Full path to output from APBS
    cfac : numeric
        Scaling factor for grid dimensions. We suggest to leave this unchanged.
    coulomb : str
        Full path to coulomb executable from APBS package.
        Must be compatible with OS.
    coulomb_results : str
        Full path to folder containing results from coulomb.
    dime : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    disu : bool, optional
        If true, Modeller will guess patches for disulfide bridges. Only
        relevant if minim is set to True.
    dx : bool
        Variable that specifies if potential files should be written to disk.
        Default is False.
    dx_files : list
        If written to disk, list of potential files written to disk.
        The folder containing these files is given by Alascan.apbs_results.
    E_ref : ndarray
        Reference state energy values for each structure from APBS.
    E_solv : ndarray
        Solvated state energy values for each structure from APBS.
    ff : string
        Force-field to use for PDB2PQR
    file_pdb_template : string
        Full path to template PDB file.
    gcent : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    Gcoul : ndarray
        Coulombic free energies, corresponding to Alascan.mutid.
    glen : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    Gref : ndarray
        Reference-state free energies, corresponding to Alascan.mutid.
    grid : numeric
        Distance spacing of grid points. If the grid dimensions are not
        divisible by three, resolution will be increased (smaller grid
        spacing) until grid dimensions are divisible by three.
    Gsolv : ndarray
        Solvated-state free energies, corresponding to Alascan.mutid.
    ion : numeric
        Ionic strength to be used in the solvated-state APBS calculations.
    jobdir : str
        [Optional] Path to folder containing results. If not specified, a
        directory will be generated.
    jobname : str
        [Optional] Name for current job, will be used to create the jobdir.
    list_chids : list
        Chain ID where mutation was made. Corresponds to Alascan.mutid.
    list_resnames : list
        Residue names where mutation was made. Corresponds to Alascan.mutid.
    list_resnums : list
        Residue numbers where mutation was made. Corresponds to Alascan.mutic.
    logs : list
        List of strings that represent the log files from every executable
        called (namely, PDB2PQR and APBS)
    logs_apbs_dir : str
        Folder in jobdir containing output from APBS (logs, input files, dx
        files)
    mask_by_sel : ndarray
        Matrix containing selection masks. The first column corresponds to
        the selection string for the parent and each column thereafter
        corresponds to an element of the selection string (selstr)
        in the same order.
    max_iter : integer, optional
        Maximum number of calls to the objective function. If this value
        is reached, then minimization is terminated. Default value is
        1000 iterations.
    min_atom_shift : float, optional
        If the maximimum atomic shift between minimization steps is less
        thant this value, convergence is reached and minimization is
        terminated. Default value is 0.1 angstroms.
    minim : bool, optional
        If true, structures will be minimzed with Modeller's conjugate
        gradient descent algorithm.
    mutation : list
        Identity of amino acid for mutation of corresponding target.
        Must be the same length as target and each residue must use the
        standard 3-letter amino acid code.
    mutid : list
        List of mutant IDs. The first element corresponds to the parent.
        Subsequent elements correspond to each element of the selection
        string list (selstr). Please use Alascan.getMutids() to get
        vectorized version.
    output : string, optional
        Modeller option specifying whether to print verbose output to
        STDOUT ('REPORT'), or to print minimal output to STDOUT ('NO_
        REPORT')
    pdb : str
        Path to PDB file with atomic coordinates. Must follow formatting
        conventions of the Protein Databank.
    pdb2pqr : str
        Full path to PDB2PQR executable.
    pdb_complex_dir : str
        Folder name in the job directory that contains the PDB file(s)
        of the complex structures.
    pdie : numeric
        Protein dielectric constant to be used in APBS.
    pqr_complex_dir : str
        Folder name in the job directory containing PQR files for each
        protein complex.
    pqr_sel_dir : list
        List of folder names in the job directory that contain PQR files
        for selection strings (selstr). Each element of pqr_sel_dir
        corresponds to the same element of selstr.
    sdie : numeric
        Solvent dielectric constant to be used in APBS.
    selstr : list
        List of selection strings. Typically each selection string will
        correspond to a chain in a protein complex. We advise users to
        specify two or more selection strings. (Ex: ['chain A', 'chain B'])
    target : list
        List of residue numbers to mutate. Must correspond element-wise
        to mutation attribute.
    """

    def __init__(self,
                 pdb,
                 target,
                 mutation,
                 pdb2pqr_exe='pdb2pqr',
                 apbs_exe='apbs',
                 coulomb_exe='coulomb',
                 selstr=['protein'],
                 jobname=None,
                 grid=1,
                 ion=0.150,
                 pdie=20.0,
                 sdie=78.54,
                 ff='parse',
                 cfac=1.5,
                 dx=False,
                 minim=True):
        """Summary

        Parameters
        ----------
        pdb : str
            The pdb file for which a computational alanine scan will be
            performed.
        target : list
            List of residue numbers to mutate. Must correspond element-wise to
            mutation attribute.
        mutation : list
            Identity of amino acid for mutation of corresponding target.
            Must be the same length as target and each residue must use the
            standard 3-letter amino acid code.
        pdb2pqr_exe : str
            Path to executable for PDB2PQR.
        apbs_exe : str
            Path to executable for APBS
        coulomb_exe : str, optional
            Path to executable for Coulomb from the APBS toolbox.
        selstr : list, optional
            List of selection strings. Typically each selection string will
            correspond to a chain in a protein complex. We advise users to
            specify two or more selection strings. The default selection will
            only investigate solvation free energies as there will be no
            binding event. (Ex: ['chain A', 'chain B'])
        jobname : str, optional
            A string representing the full path to the directory where
            results will be stored.
        grid : int, optional
            Desired grid spacing to be used in APBS.
            Depending on the PDB file, actual grid spacing may be a finer
            resolution. Defaults to 1 angstrom.
        ion : float, optional
            Ionic strength (M) to be used in APBS.
        pdie : float, optional
            Protein dielectric constant to be used in APBS.
        sdie : float, optional
            Solvent dielectric constant to be used in APBS>
        ff : str, optional
            Force field to assign charges to individual atoms.
            We recommend using PARSE for this analysis, the default value.
        cfac : float, optional
            Scaling factor for grid dimension.
            This generally should not be changed without a reason.
        dx : bool, optional
            Boolean flag stating if potential files should be written to
            disk (True) or not (False).
        minim : bool, optional
            If true, structures will be minimized with Modeller's conjugate
            gradient descent algorithm. Default is True for the
            DirectedMutagenesis class as mutations may result in
            unfavorable conformations.
        """
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
        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        if jobname is None:
            jobname = '%4d%02d%02d_%02d%02d%02d' % (
                dt.date.today().year, dt.date.today().month,
                dt.date.today().day, dt.datetime.now().hour,
                dt.datetime.now().minute, dt.datetime.now().second)
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

        self.genDirs()
        self.genMutid()
        self.genParent()
        self.find_grid()

        self.minim = minim
        self.disu = True
        self.min_atom_shift = 0.1
        self.max_iter = 1000
        self.output = 'NO_REPORT'
        self.logs = []

    def getMutids(self):
        """Summary

        Returns
        -------
        list
            Returns vectorized format of mutids.
        """
        l = self.mutid
        return [item for sublist in l for item in sublist]

    def genDirs(self):
        """Summary
        This subroutine will generate all directories needed to contain
        structural files, logs, etc. In the future we may implement a method
        to remove such files when outputs are more standardized.

        Returns
        -------
        None
            No output, operates on the class and generates folders in the
            job directory.
        """
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
            if not os.path.exists(
                    os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1))):
                os.makedirs(os.path.join(self.jobdir, 'seg%d_pqr' % (i + 1)))
        self.pqr_complex_dir = pqr_complex_dir
        self.pqr_sel_dir = pqr_sel_dir

        # Create necessary directories for APBS files
        logs_apbs_dir = 'apbs_logs'
        if not os.path.exists(os.path.join(self.jobdir, logs_apbs_dir)):
            os.makedirs(os.path.join(self.jobdir, logs_apbs_dir))
        self.logs_apbs_dir = 'apbs_logs'

    def genMutid(self):
        """Summary
        This subroutine reads the input PDB, selects the structure where
        mutations will occur, and saves all mutant IDs in the class.
        If region is specified in the constructor, then the constraint will
        be applied here.

        Returns
        -------
        None
            Operates on the class to generate a list of mutant IDs for each
            selection string. The following class variables will be
            generated (see class description):
                - mutid
                - list_chids
                - list_resnums
                - list_resnames
                - mask_by_sel
        """
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
        list_chids[0] = list([''])
        list_resnums[0] = [np.zeros(0).tolist()]
        list_resnames[0] = list([''])

        for i, region in zip(xrange(len(selstr)), selstr):
            for j, sel, mut in zip(xrange(len(target)), target, mutation):
                combined_selection = parent_pdb.select('(' + ' and '.join(
                    ['(' + sel + ')', '(' + region + ')', 'calpha']) + ')')

                if combined_selection is not None:
                    chids = combined_selection.getChids().tolist()
                    resnums = combined_selection.getResnums().tolist()
                    resnames = combined_selection.getResnames().tolist()

                    label = 'sel%d' % (i + 1)
                    for resnum, resname in zip(resnums, resnames):
                        label = '_'.join([
                            label,
                            AA_dict[resname] + str(resnum) + AA_dict[mut]
                        ])

                    list_mutids[i + 1].append(label)
                    list_chids[i + 1].append(chids)
                    list_resnums[i + 1].append(resnums)
                    list_resnames[i + 1].append(resnames)
                    mask_by_sel[j + 1, i + 1] = True

        self.mutid = list_mutids
        self.list_chids = list_chids
        self.list_resnums = list_resnums
        self.list_resnames = list_resnames
        self.mask_by_sel = mask_by_sel

    def genParent(self):
        """Summary
        Reads PDB file specified in the constructor; applies and combines
        results from the selection strings; and saves the final template
        structure in the job directory.

        Returns
        -------
        None
            Template pdb written in job directory and location saved in
            Alascan.file_pdb_template.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        infile = os.path.join(jobdir, pdb_complex_dir,
                              parent_file_prefix + '.pdb')
        system = parent_pdb.select('((' + ') or ('.join(selstr) + '))')
        # print '(('+') or ('.join(selstr)+'))'
        pd.writePDB(infile, system)

        self.file_pdb_template = infile

    def find_grid(self):
        """Summary
        Calculate grid dimensions for APBS (dime, glen, gcent)

        Returns
        -------
        TYPE
            Sets class attributes dime, glen, and gcent
        """
        pdbfile = self.file_pdb_template
        grid = self.grid
        cfac = self.cfac

        glen = np.zeros(3)

        # Find bounding box
        pdb = pd.parsePDB(pdbfile)
        max_xyz = pdb.getCoords().max(axis=0)
        min_xyz = pdb.getCoords().min(axis=0)

        # Determine mesh dimensions according to Ron's AESOP protocol in the R
        # source file
        fg = np.ceil(max_xyz - min_xyz)
        fg = np.ceil((fg + 5) * cfac)
        glen = np.vstack((glen, fg)).max(axis=0)
        dime_list = (32 * np.linspace(1, 100, 100)) + \
            1  # list of possible dime values
        # index of dime to use from list, subtract one to be consistent with
        # python indexing!
        dime_ind = np.ceil(glen / (32 * grid)) - 1
        dime = np.array(
            (dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])],
             dime_list[int(dime_ind[2])]))
        ix = 0
        iy = 0
        iz = 0
        counter = 0
        while ((dime[0] * dime[1] * dime[2] % 3 != 0) or (counter <= 5)):
            ix += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iy += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iz += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            counter += 1

        self.dime = dime
        self.glen = glen
        self.gcent = pd.calcCenter(pdb).astype(int)

    def set_grid(self, dime, glen, gcent):
        """Summary
        In the case that the user wishes to manually specify grid dimension,
        this may be accomplished with the set_grid method.
        Typically, this is used when grid space parameters must be consistent
        for many analyses. Please see description
        at: http://www.poissonboltzmann.org/docs/apbs-overview/ for
        description of parameters.

        Parameters
        ----------
        dime : list
            List of three integers.
        glen : list
            List of three integers.
        gcent : list
            List of three integers.

        Returns
        -------
        TYPE
            Sets class attributes for dime, glen, gcent
        """
        self.dime = dime
        self.glen = glen
        self.gcent = gcent

    def genPDB(self, minim=True):
        """Summary
        Generates mutations by calling function to mutate PDB with modeller

        Returns
        -------
        None
            Write PDB library to expected path according to class attributes.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pdb_complex_dir = self.pdb_complex_dir
        mutation = self.mutation

        parent_file_prefix = 'wt'
        parent_pdb = pd.parsePDB(self.pdb)

        list_mutids = [item for sublist in self.mutid for item in sublist]
        list_chids = [item for sublist in self.list_chids for item in sublist]
        list_resnums = [
            item for sublist in self.list_resnums for item in sublist
        ]

        infile = os.path.join(jobdir, pdb_complex_dir,
                              parent_file_prefix + '.pdb')
        system = parent_pdb.select('(' + ') or ('.join(selstr) + ')')
        pd.writePDB(infile, system)

        if self.minim == True:
            minimize_cg(
                infile,
                dest=infile,
                disu=self.disu,
                min_atom_shift=self.min_atom_shift,
                max_iter=self.max_iter,
                output=self.output)
        for mutid, chain, resnum, mut in zip(list_mutids[1:], list_chids[1:],
                                             list_resnums[1:], mutation):
            outpath = os.path.join(jobdir, pdb_complex_dir, mutid)
            print '\n%s:\tgenerating PDB for mutant: %s' % (self.jobname,
                                                            mutid)
            mutatePDB(
                pdb=infile,
                mutid=outpath,
                resnum=resnum,
                chain=chain,
                resid=mut)
            if self.minim == True:
                minimize_cg(
                    outpath + '.pdb',
                    dest=outpath + '.pdb',
                    disu=self.disu,
                    min_atom_shift=self.min_atom_shift,
                    max_iter=self.max_iter,
                    output=self.output)

    def genPQR(self):
        """Summary
        Generates PQR for each PDB in library.

        Returns
        -------
        None
            Calls PDB2PQR and writes PQR to expected path according to
            class attributes.
        """
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
            print '\n%s:\tgenerating PQR for mutant: %s' % (self.jobname,
                                                            mutid)
            infile = os.path.join(jobdir, pdb_complex_dir, mutid + '.pdb')
            outfile = os.path.join(jobdir, pqr_complex_dir, mutid + '.pqr')
            (pqr_log, pqr_errs) = execPDB2PQR(
                path_pdb2pqr, infile, outfile=outfile, ff=ff)
            self.logs.append(pqr_log)
            complex_pqr = pd.parsePQR(outfile)
            for sel, seldir in zip(selstr, pqr_sel_dir):
                selfile = os.path.join(jobdir, seldir, mutid + '.pqr')
                pqr = complex_pqr.select(sel)
                pd.writePQR(selfile, pqr)

    def calcAPBS(self):
        """Summary
        Call apbs to calculate reference-state and solvation-state energies
        for all structures in library.

        Returns
        -------
        None
            Sets class attributes Gsolv and Gref
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel))
        Gref = np.zeros((dim_mutid, dim_sel))

        complex_pqr = os.path.join(jobdir, pqr_complex_dir,
                                   list_mutids[0] + '.pqr')
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            print '\n%s:\tcalculating solvation and reference' \
                ' energies for mutant: %s' % (self.jobname, mutid)
            # complex_pqr = os.path.join(jobdir, pqr_complex_dir, mutid+'.pqr')
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_prefix_log = os.path.join(jobdir, logs_apbs_dir, mutid)
                energies, log = execAPBS(
                    path_apbs,
                    subunit_pqr,
                    self.dime,
                    self.glen,
                    self.gcent,
                    prefix=path_prefix_log,
                    ion=self.ion,
                    pdie=self.pdie,
                    sdie=self.sdie,
                    dx=self.dx)
                self.logs.append(log)
                Gsolv[i, j] = energies[0][0]
                Gref[i, j] = energies[0][1]

        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcAPBS_parallel(self, n_workers=None):
        """Summary
        Run APBS on each structure in mutant library, in parallel.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all
            available threads.

        Returns
        -------
        None
            Saves solvated-state and reference-state free energies
            as class attributes.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        logs_apbs_dir = self.logs_apbs_dir
        path_apbs = self.apbs

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
        mask_by_sel[0, :] = np.ones(dim_sel).astype(bool)
        mask_by_sel[:, 0] = np.ones(dim_mutid).astype(bool)

        Gsolv = np.zeros((dim_mutid, dim_sel)).astype(float)
        Gref = np.zeros((dim_mutid, dim_sel)).astype(float)

        path_list = []
        pqr_chain_list = []
        dime_list = []
        glen_list = []
        gcent_list = []
        prefix_list = []
        ion_list = []
        pdie_list = []
        sdie_list = []
        dx_list = []
        i_list = []
        j_list = []

        # Find all calculations to be done
        complex_pqr = os.path.join(jobdir, pqr_complex_dir,
                                   list_mutids[0] + '.pqr')
        for i, mutid in zip(xrange(dim_mutid), list_mutids):
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_list.append(path_apbs)
                pqr_chain_list.append(subunit_pqr)
                dime_list.append(self.dime)
                glen_list.append(self.glen)
                gcent_list.append(self.gcent)
                prefix_list.append(
                    os.path.join(
                        jobdir, logs_apbs_dir, '%d_%d_' % (i, j) +
                        mutid))  # added to make sure apbs.in file is unique!
                ion_list.append(self.ion)
                pdie_list.append(self.pdie)
                sdie_list.append(self.sdie)
                dx_list.append(self.dx)
                i_list.append(i)
                j_list.append(j)

        # Organize kernel and run batch process
        kernel = zip(path_list, pqr_chain_list, dime_list, glen_list,
                     gcent_list, prefix_list, ion_list, pdie_list, sdie_list,
                     dx_list, i_list, j_list)
        p = Pool(n_workers)
        print '%s:\trunning batchAPBS ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchAPBS, kernel):
            counter += 1
            print '.... %s:\tbatch APBS %d percent complete ....' % (
                self.jobname, int(counter * 100 / max_count))
            data, log = result
            i = int(data[0])
            j = int(data[1])
            solv = data[2]
            ref = data[3]
            Gsolv[i, j] = solv
            Gref[i, j] = ref
            self.logs.append(log)
        if self.dx == True:
            self.dx_files = [x + '.dx' for x in prefix_list]

        self.Gsolv = Gsolv
        self.Gref = Gref

    def calcCoulomb(self):
        """Summary
        Calculates Coulombic free energies with coulomb.exe from the APBS
        toolbox.

        Returns
        -------
        None
            Saves Coulombic free energies as a class attribute.
        """
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
            print '\n%s:\tcalculating coulombic energies for mutant: %s' % (
                self.jobname, mutid)
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                energies, log = execCoulomb(path_coulomb, subunit_pqr)
                self.logs.append(log)
                Gcoul[i, j] = energies / pdie

        self.Gcoul = Gcoul

    def calcCoulomb_parallel(self, n_workers=None):
        """Summary
        Calculates Coulombic free energies with coulomb.exe from the APBS
        toolbox in a parallel manner.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all
            available threads.

        Returns
        -------
        None
            Saves Coulombic free energies as a class attribute.
        """
        selstr = self.selstr
        jobdir = self.jobdir
        pqr_complex_dir = self.pqr_complex_dir
        pqr_sel_dir = self.pqr_sel_dir
        path_coulomb = self.coulomb

        list_mutids = self.getMutids()

        dim_mutid = len(list_mutids)
        dim_sel = len(selstr) + 1

        # Mask parts that are true will be run with APBS
        mask_by_sel = np.copy(self.mask_by_sel)
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
            for j, seldir in zip(
                    xrange(dim_sel), [pqr_complex_dir] + pqr_sel_dir):
                subunit_pqr = os.path.join(jobdir, seldir, mutid + '.pqr')
                path_list.append(path_coulomb)
                pqr_chain_list.append(subunit_pqr)
                pdie_list.append(self.pdie)
                i_list.append(i)
                j_list.append(j)

        # Organize kernel and run batch process
        kernel = zip(path_list, pqr_chain_list, pdie_list, i_list, j_list)
        p = Pool(n_workers)
        print '%s:\trunning batchCoulomb ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchCoulomb, kernel):
            counter += 1
            print '.... %s:\tbatch coulomb %d percent complete ....' % (
                self.jobname, int(counter * 100 / max_count))
            data, log = result
            self.logs.append(log)
            i = int(data[0])
            j = int(data[1])
            coul = data[2]
            Gcoul[i, j] = coul

        self.Gcoul = Gcoul

    def ddGa_rel(self):
        """Summary
        Calculates and returns the free energy of association relative
        to the parent free energy of association.

        Returns
        -------
        ndarray
            Array of free energies corresponding to the mutant IDs from
            the Alascan.getMutIDs() method.
        """
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
        """Summary
        Calculates and returns the free energy of a solvation for each mutant
        relative to the parent free energy of solvation.

        Returns
        -------
        ndarray
            Array of solvation free energies corresponding to mutant IDs
            from the Alascan.getMutIDs() method.
        """
        Gsolv = self.Gsolv
        Gref = self.Gref
        dGsolv = Gsolv - Gref
        dGsolv = dGsolv - dGsolv[0, 0]
        return dGsolv[:, 0]

    def run(self):
        """Summary
        Perform a directed mutagenesis scan on the initialized class.

        Returns
        -------
        None
            Outputs text to STDOUT when run is complete, will be made optional
            in the future.
        """
        start = ti.default_timer()
        self.logs = []
        self.genPDB()
        self.genPQR()
        self.calcAPBS()
        self.calcCoulomb()
        stop = ti.default_timer()
        print '%s:\tAESOP directed mutagenesis scan completed' \
            ' in %.2f seconds' % (self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def run_parallel(self, n_workers=None):
        """Summary
        Perform a computational directed mutagenesis scan on the initialized
        class using multiple processes in parallel.

        Parameters
        ----------
        n_workers : int
            Number of processes to run. If None, method will use all
            available threads.

        Returns
        -------
        None
            Outputs text to STDOUT when run is complete, will be made
            optional in the future.
        """
        start = ti.default_timer()
        self.logs = []
        self.genPDB()
        self.genPQR()
        self.calcAPBS_parallel()
        self.calcCoulomb_parallel()
        stop = ti.default_timer()
        print '%s:\tAESOP directed mutagenesis scan completed' \
            ' in %.2f seconds' % (self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def summary(self, filename=None):
        """Summary
        Summarize results from the computational alanine scan once complete.

        Parameters
        ----------
        filename : str, optional
            In order to write summary to a text file, supply the
            filename (full path).

        Returns
        -------
        None
            Prints summary of residues and energies relative to the parent
            structure if no filename is provided.
            Otherwise, writes to text file.
        """
        selstr = self.selstr
        mutids = self.getMutids()
        energies = None
        if len(selstr) > 1:
            energies = self.ddGa_rel()
        elif len(selstr) == 1:
            energies = self.dGsolv_rel()
        lines = ['%s, %f' % (lbl, val) for lbl, val in zip(mutids, energies)]
        lines.append('\nKey\n---\n')
        for i, seg in enumerate(selstr):
            lines.append('sel%d = %s' % (int(i + 1), str(seg)))
        if filename is None:
            print(lines)
        if filename is not None:
            np.savetxt(filename, lines, fmt='%s')

    def writeLogs(self, filename=None):
        if filename is None:
            filename = os.path.join(self.jobdir, 'AESOP_logs.txt')
        with open(filename, 'w') as f:
            logs = '\n==== Log Instance ====\n'.join(self.logs)
            f.write(logs)

    def viewLogs(self):
        logs = '\n==== Log Instance ====\n'.join(self.logs)
        print logs

    def checkwarnings(self):
        lines = ''.join(self.logs)
        matches = re.findall('[W][Aa][Rr][Nn][Ii][Nn][Gg]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1

    def checkerrors(self):
        lines = ''.join(self.logs)
        matches = re.findall('[E][Rr][Rr][Oo][Rr]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1


##########################################################################
# Container for performing ESD analysis on set of PDB files
#   alascan     -   Alascan class with certain class functions.
##########################################################################


class ElecSimilarity:  # PLEASE SUPERPOSE SYSTEM BEFORE USING THIS METHOD!
    # Coordinates must be consistent!
    """Summary

    Attributes
    ----------
    apbs : str
        Full path to APBS executable. Must be compatible with OS.
    cfac : numeric
        Scaling factor for grid dimensions. We suggest to leave this unchanged.
    dim_dx : list
        Dimensions of grid space
    dime : list
        List of three integers. Parameter required for APBS.Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    disu : bool, optional
        If True, Modeller will guess the patches for disulfide bridges
        within the provided protein structures. Only relevant if minim
        is set to True.
    dx : bool
        Variable that specifies if potential files should be written to disk.
        Default is False.
    dxdir : str
        Folder in job directory where potential files are stored.
    dxfiles : list
        List of all potential files.
    edges : ndarray
        Edges of grid space.
    esd : ndarray
        Matrix of pairwise electrostatic similarities.
    ff : str
        Forcefield to use in assigning charges to PDB files.
    gcent : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    glen : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    grid : int
        Desired grid spacing in Angstroms. Actual grid spacing may be
        slightly lower.
    ids : list
        List of IDs for each structure being compared with the ESD metric.
    ion : numeric
        Ionic strength to be used in the solvated-state APBS calculations.
    jobdir : TYPE
        Description
    jobname : str
        [Optional] Name for current job, will be used to create the jobdir.
    logs : list
        List of strings that represent the log files from every executable
        called (namely, PDB2PQR and APBS)
    max_iter : integer, optional
        Maximum number of calls to the objective function. If this value
        is reached, then minimization is terminated. Default value is
        1000 iterations.
    midpoints : ndarray
        Midpoints of grid space
    min_atom_shift : float, optional
        If the maximimum atomic shift between minimization steps is less
        thant this value, convergence is reached and minimization is
        terminated. Default value is 0.1 angstroms.
    minim : bool, optional
        If true, structures will be minimzed with Modeller's conjugate
        gradient descent algorithm.
    output : string, optional
        Modeller option specifying whether to print verbose output to
        STDOUT ('REPORT'), or to print minimal output to STDOUT ('NO_
        REPORT')
    pdb2pqr : str
        Full path to PDB2PQR executable.
    pdbdir : str
        Folder in job directory where PDB files are located.
    pdbfiles : list
        List of PDB file names.
    pdie : numeric
        Protein dielectric constant to be used in APBS.
    pqrdir : str
        Folder in job directory containing PQR files.
    pqrfiles : list
        List of PQR files in pqrdir
    sdie : numeric
        Solvent dielectric constant to be used in APBS.
    """

    def __init__(self,
                 pdbfiles,
                 pdb2pqr_exe='pdb2pqr',
                 apbs_exe='apbs',
                 selstr=None,
                 jobname=None,
                 grid=1,
                 ion=0.150,
                 pdie=20.0,
                 sdie=78.54,
                 ff='parse',
                 cfac=1.5,
                 minim=False):
        """Summary
        Constructor for ElecSimilarity class. Responsible for preliminary
        parameterization.

        Parameters
        ----------
        pdbfiles : list
            List of file names for all PDBs to be compared.
            Each string should be a full path if not in the current
            working directory.
        pdb2pqr_exe : str
            Full path to PDB2PQR executable.
        apbs_exe : TYPE
            Full path to APBS executable. Must be compatible with OS.
        selstr : list, optional
            List of selection strings. Each element should correspond
            element-wise with the PDB from pdbfiles.
        jobname : None, optional
            [Optional] Name for current job, will be used to create the jobdir.
        grid : int, optional
            Distance spacing of grid points. If the grid dimensions are not
            divisible by three, resolution will be
            increased (smaller grid spacing) until grid dimensions are
            divisible by three.
        ion : float, optional
            Ionic strength (M) to be used in APBS.
        pdie : float, optional
            Protein dielectric constant to be used in APBS.
        sdie : float, optional
            Solvent dielectric constant to be used in APBS.
        ff : str, optional
            Force field to assign charges to individual atoms.
            We recommend using PARSE for this analysis, the default value.
        cfac : float, optional
            Scaling factor for grid dimensions.
            We suggest to leave this unchanged.
        minim : bool, optional
            If TRUE, energy minimization of protein structures will be
            performed with Modeller's conjugate gradient descent
            algorithm. Default is False for the Alanine Scan class as
            mutations should not result in atomic clashes.
        """
        self.pdbfiles = [os.path.basename(pdbfile) for pdbfile in pdbfiles]
        self.ids = [
            os.path.splitext(os.path.basename(pdbfile))[0]
            for pdbfile in pdbfiles
        ]
        self.pdb2pqr = pdb2pqr_exe
        self.apbs = apbs_exe
        self.dx = True

        if jobname is None:
            jobname = '%4d%02d%02d_%02d%02d%02d' % (
                dt.date.today().year, dt.date.today().month,
                dt.date.today().day, dt.datetime.now().hour,
                dt.datetime.now().minute, dt.datetime.now().second)
        self.jobname = jobname
        self.jobdir = jobname
        if not os.path.exists(os.path.join(self.jobdir)):
            os.makedirs(os.path.join(self.jobdir))
        self.pdbdir = os.path.join(self.jobdir, 'pdb_files')
        self.pqrdir = os.path.join(self.jobdir, 'pqr_files')
        self.dxdir = os.path.join(self.jobdir, 'dx_files')
        if not os.path.exists(self.pdbdir):
            os.makedirs(self.pdbdir)
        if not os.path.exists(self.pqrdir):
            os.makedirs(self.pqrdir)
        if not os.path.exists(self.dxdir):
            os.makedirs(self.dxdir)

        for i, pdbfile in zip(xrange(len(pdbfiles)), pdbfiles):
            pdb = pd.parsePDB(pdbfile)
            if selstr is None:
                pd.writePDB(
                    os.path.join(self.pdbdir, os.path.basename(pdbfile)), pdb)
            elif selstr is not None:
                pd.writePDB(
                    os.path.join(self.pdbdir, os.path.basename(pdbfile)),
                    pdb.select(selstr[i]))

        self.grid = grid
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        self.ff = ff
        self.cfac = cfac
        self.logs = []
        self.pqrfiles = []

        self.minim = minim
        self.disu = True
        self.min_atom_shift = 0.1
        self.max_iter = 1000
        self.output = 'NO_REPORT'

    def centerPDB(self):
        """Summary
        Re-reads PDB file in pdbdir and centers coordinates at (0, 0, 0)

        Returns
        -------
        TYPE
            Overwrites previous PDB files in pdbdir
        """
        pdbdir = self.pdbdir
        pdbfiles = self.pdbfiles
        for pdbfile in pdbfiles:
            pdb = pd.parsePDB(os.path.join(pdbdir, pdbfile))
            pd.moveAtoms(pdb, to=np.zeros(3))
            pd.writePDB(os.path.join(pdbdir, pdbfile), pdb)

    def minimPDB(self):
        """Summary
        Re-reads PDB files in pdbdir and performs energy minimization.

        Returns
        -------
        TYPE
            Overwrites previous PDB files in pdbdir
        """
        pdbdir = self.pdbdir
        pdbfiles = self.pdbfiles
        for pdbfile in pdbfiles:
            fname = os.path.join(pdbdir, pdbfile)
            minimize_cg(
                fname,
                dest=fname,
                disu=self.disu,
                min_atom_shift=self.min_atom_shift,
                max_iter=self.max_iter,
                output=self.output)

    def superposePDB(self):
        """Summary
        Superpose each structure in pdbfiles with first element in
        pdbfiles list. This uses Modeller to perform the superpositioning.

        Returns
        -------
        TYPE
            Overwrites PDB files in pdbdir with new coordinate information.
        """
        pdbdir = self.pdbdir
        pdbfiles = self.pdbfiles

        reffile = pdbfiles[0]
        for pdbfile in pdbfiles[1:]:
            ref = os.path.join(pdbdir, reffile)
            pdb = os.path.join(pdbdir, pdbfile)
            superpose(ref=ref, pdb=pdb, atype='CA', output=None)

    def initializeGrid(self):
        """Summary
        Method to find grid parameters and ensure that the product of
        dimensions is divisible by three.

        Returns
        -------
        None
            Sets class attributes dime, glen, gcent.
        """
        pdbdir = self.pdbdir
        pdbfiles = self.pdbfiles
        grid = self.grid
        cfac = self.cfac

        glen = np.zeros(3)

        # Find bounding box
        ref = pd.parsePDB(os.path.join(pdbdir, pdbfiles[0]))
        max_xyz = ref.getCoords().max(axis=0)
        min_xyz = ref.getCoords().min(axis=0)
        for pdbfile in pdbfiles:
            pdb = pd.parsePDB(os.path.join(pdbdir, pdbfile))
            max_xyz = np.vstack(
                (max_xyz, pdb.getCoords().max(axis=0))).max(axis=0)
            min_xyz = np.vstack(
                (min_xyz, pdb.getCoords().min(axis=0))).min(axis=0)

        # Determine mesh dimensions according to Ron's AESOP protocol in the R
        # source file
        fg = np.ceil(max_xyz - min_xyz)
        fg = np.ceil((fg + 5) * cfac)
        glen = np.vstack((glen, fg)).max(axis=0)
        dime_list = (32 * np.linspace(1, 100, 100)) + \
            1  # list of possible dime values
        # index of dime to use from list, subtract one to be consistent with
        # python indexing!
        dime_ind = np.ceil(glen / (32 * grid)) - 1
        dime = np.array(
            (dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])],
             dime_list[int(dime_ind[2])]))
        ix = 0
        iy = 0
        iz = 0
        counter = 0
        while ((dime[0] * dime[1] * dime[2] % 3 != 0) or (counter <= 5)):
            ix += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iy += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iz += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            counter += 1

        self.dime = dime
        self.glen = glen
        self.gcent = pd.calcCenter(ref).astype(int)

    def genPQR(self):
        """Summary
        Convert all PDB files to PQR files with charges allocated according to
        a compatible force-field

        Returns
        -------
        None
            Generates PQR files in the pqrdir
        """
        pdbdir = self.pdbdir
        pqrdir = self.pqrdir
        pdbfiles = self.pdbfiles
        path_pdb2pqr = self.pdb2pqr
        ff = self.ff

        self.pqrfiles = []
        for pdbfile in pdbfiles:
            infile = os.path.join(pdbdir, pdbfile)
            pqrname = os.path.splitext(pdbfile)[0] + '.pqr'
            outfile = os.path.join(pqrdir, pqrname)
            print 'Converting %s to PQR' % (pdbfile)
            (pqr_log, pqr_errs) = execPDB2PQR(
                path_pdb2pqr, infile, outfile=outfile, ff=ff)
            self.logs.append(pqr_log)
            self.pqrfiles.append(pqrname)

    def mutatePQR(self, selstr=['protein'], minim=False):
        """Summary
        Mutate all PQR files, optional method

        Returns
        -------
        None
            Generates PQR files in the pqrdir
        """
        pqrdir = self.pqrdir
        pqrfiles = self.pqrfiles
        sel = ' or '.join(selstr)
        sel = ' and '.join(['(' + sel + ')', 'charged', 'calpha'])

        mutantpqrs = []
        for pqrfile in pqrfiles:
            infile = os.path.join(pqrdir, pqrfile)
            pqr = pd.parsePQR(infile)
            pqrsel = pqr.select(sel)

            resnums = pqrsel.getResnums().tolist()
            resnames = [AA_dict[x] for x in pqrsel.getResnames().tolist()]
            chains = pqrsel.getChids().tolist()

            for resnum, chain, resname in zip(resnums, chains, resnames):
                mutid = os.path.splitext(
                    os.path.basename(infile))[0] + '_%s_%s%sA' % (
                        str(chain), str(resname), str(resnum))
                mutantpqrs.append(mutid + '.pqr')
                print 'Generating mutant: %s' % (mutid)
                mutatePQR(
                    pqrfile=infile,
                    mutid=os.path.join(pqrdir, mutid),
                    resnum=resnum,
                    chain=chain)
                if (self.minim == True) or (minim == True):
                    minimize_cg(
                        os.path.join(pqrdir, mutid + '.pqr'),
                        dest=os.path.join(pqrdir, mutid + '.pqr'),
                        disu=self.disu,
                        min_atom_shift=self.min_atom_shift,
                        max_iter=self.max_iter,
                        output=self.output)
        self.pqrfiles = pqrfiles + mutantpqrs
        self.ids = [
            os.path.splitext(os.path.basename(x))[0] for x in self.pqrfiles
        ]

    def genDX(self):
        """Summary
        Generates potential files using APBS.

        Returns
        -------
        None
            Generates DX files in dxdir
        """
        path_apbs = self.apbs
        pdbfiles = self.pqrfiles
        pqrdir = self.pqrdir
        dxdir = self.dxdir
        grid = self.grid
        ion = self.ion
        pdie = self.pdie
        sdie = self.sdie
        cfac = self.cfac
        glen = self.glen
        gcent = self.gcent
        dime = self.dime

        pqrfiles = [
            os.path.join(pqrdir, os.path.splitext(pdbfile)[0] + '.pqr')
            for pdbfile in pdbfiles
        ]
        apbsfiles = [
            os.path.join(dxdir, os.path.splitext(pdbfile)[0])
            for pdbfile in pdbfiles
        ]
        for pqrfile, apbsfile in zip(pqrfiles, apbsfiles):
            log = calcDX(
                path_apbs,
                pqrfile,
                prefix=apbsfile,
                grid=grid,
                ion=ion,
                pdie=pdie,
                sdie=sdie,
                cfac=cfac,
                glen=glen,
                gcent=gcent,
                dime=dime)
            self.logs.append(log)

    def genDX_parallel(self, n_workers=None):
        """Summary
        Generates multiple potential files in parallel by calling APBS
        multiple times according to how many threads are
        available/specified.

        Parameters
        ----------
        n_workers : int
            Number of processes to run.
            If None, method will use all available threads.

        Returns
        -------
        TYPE
            Generates DX files in dxdir.
        """
        pdbfiles = self.pqrfiles
        pqrdir = self.pqrdir
        dxdir = self.dxdir
        grid = self.grid
        ion = self.ion
        pdie = self.pdie
        sdie = self.sdie
        cfac = self.cfac
        glen = self.glen
        gcent = self.gcent
        dime = self.dime

        if n_workers is None:
            n_workers = cpu_count()

        path_apbs = [self.apbs for pdbfile in pdbfiles]
        pqrfiles = [
            os.path.join(pqrdir, os.path.splitext(pdbfile)[0] + '.pqr')
            for pdbfile in pdbfiles
        ]
        apbsfiles = [
            os.path.join(dxdir, os.path.splitext(pdbfile)[0])
            for pdbfile in pdbfiles
        ]
        grid = [grid for pdbfile in pdbfiles]
        ion = [ion for pdbfile in pdbfiles]
        pdie = [pdie for pdbfile in pdbfiles]
        sdie = [sdie for pdbfile in pdbfiles]
        cfac = [cfac for pdbfile in pdbfiles]
        glen = [glen for pdbfile in pdbfiles]
        gcent = [gcent for pdbfile in pdbfiles]
        dime = [dime for pdbfile in pdbfiles]
        kernel = zip(path_apbs, pqrfiles, apbsfiles, grid, ion, pdie, sdie,
                     cfac, glen, gcent, dime)

        p = Pool(n_workers)
        print '%s:\trunning batchCalcDX ....' % (self.jobname)
        counter = 0
        max_count = len(kernel)
        for result in p.imap_unordered(batchCalcDX, kernel):
            counter += 1
            print '.... %s:\tbatch coulomb %d percent complete ....' % (
                self.jobname, int(counter * 100 / max_count))
            self.logs.append(result)

    def calcESD(self, method='AND'):
        """Summary
        Compare potential files and calculate the similarity distance.
        Smaller distances imply similarity.

        Parameters
        ----------
        method : str, optional
            This parameter will allow for other metrics to compare
            grid potentials; however, for now only 'AND' is implemented.

        Returns
        -------
        None
            Stores esd matrix as class attribute.
        """

        def symmetrize(a):
            """Summary

            Parameters
            ----------
            a : TYPE
                Description

            Returns
            -------
            TYPE
                Description
            """
            return a + a.T - np.diag(a.diagonal())

        pdbfiles = self.pqrfiles
        dxdir = self.dxdir

        self.dxfiles = [
            glob.glob(os.path.join(
                dxdir, os.path.splitext(os.path.basename(pdbfile))[0] + '*.dx'))[0]
            for pdbfile in pdbfiles
        ]
        files = self.dxfiles
        ids = self.ids

        grid = Grid(files[0])
        self.dim_dx = grid.pot.size

        dim = self.dim_dx
        esd = np.zeros((len(ids), len(ids)))

        indices = it.combinations(range(len(ids)), 2)
        for i, j in indices:
            a = Grid(files[i]).pot.reshape((dim, ))
            b = Grid(files[j]).pot.reshape((dim, ))
            if method is 'AND':
                diff = np.abs(a - b)
                maxpot = np.abs(np.vstack((a, b))).max(axis=0)
                val = np.nan_to_num(np.divide(diff, maxpot)).sum() / dim
                esd[i, j] = val
        esd = symmetrize(esd)
        self.esd = esd

    def calcESI(self, method='AND', idx=0):
        """Summary

        Compare potential files and calculate the similarity index.
        Values closer to 1 imply similarity while values closer to zero imply
        dissimilarity.

        Parameters
        ----------
        method : str, optional
            This parameter will allow for other metrics to compare
            grid potentials; however, for now only 'AND' is implemented.
        idx : int
            Index of original PDB files supplied containing reference
            structure. Set to None to perform all pairwise comparisons.

        Returns
        -------
        None
            Writes esi files to the esi_files directory within the job
            directory.
        """
        pdbfiles = self.pqrfiles
        dxdir = self.dxdir

        esidir = os.path.join(self.jobdir, 'esi_files')
        self.esidir = esidir
        if not os.path.exists(os.path.join(self.esidir)):
            os.makedirs(os.path.join(self.esidir))

        self.dxfiles = [
            glob.glob(os.path.join(
                dxdir, os.path.splitext(os.path.basename(pdbfile))[0] + '*.dx'))[0]
            for pdbfile in pdbfiles
        ]
        files = self.dxfiles
        ids = self.ids

        n = len(files) - 1
        if idx is None:
            idx = range(n)
        elif isinstance(idx, (int, long)):
            idx = [idx]

        esifiles = []
        esilist = []
        for i in idx:
            ref_name = ids[i]
            ref = Grid(files[i])
            dim = ref.pot.size

            esi = []
            filename = os.path.join(esidir, ref_name + '.dx')
            for j in xrange(n):

                if (method is 'AND') and (i != j):
                    dat = Grid(files[j])
                    a = ref.pot.astype(float).reshape((dim, ))
                    b = dat.pot.astype(float).reshape((dim, ))

                    diff = np.abs(a - b)
                    maxpot = np.abs(np.vstack((a, b))).max(axis=0)
                    val = np.divide(diff, maxpot)
                    val = np.nan_to_num(val)

                    esi.append(val)

            esi = np.vstack(esi)
            esi = np.ones(esi.shape) - esi
            esi = np.sum(esi, axis=0) / n
            ref.pot = esi.reshape((dim / 3, 3))
            ref.writeDX(filename)
            esifiles.append(filename)
            esilist.append(esi)
        esilist = np.vstack(esilist)
        self.esifiles = esifiles
        self.esi = esilist

    def run(self,
            center=False,
            superpose=False,
            esi=False,
            esd=True,
            selstr=None,
            idx=0,
            minim=False):
        start = ti.default_timer()
        self.logs = []
        if center:
            self.centerPDB()
        if self.minim or minim:
            self.minimPDB()
        if superpose:
            self.superposePDB()
        self.initializeGrid()
        self.genPQR()
        if selstr is not None:
            self.mutatePQR(selstr=selstr)
        if len(self.pdbfiles) == 1 and selstr is None:
            self.mutatePQR()
        self.genDX()
        if esd:
            self.calcESD()
        if esi:
            self.calcESI(idx=idx)
        stop = ti.default_timer()
        print '%s:\tAESOP electrostatic similarity comparison ' \
            'completed in %.2f seconds' % (self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def run_parallel(self,
                     n_workers=None,
                     center=False,
                     superpose=False,
                     esi=False,
                     esd=True,
                     selstr=None,
                     idx=0,
                     minim=False):
        start = ti.default_timer()
        self.logs = []
        if center:
            self.centerPDB()
        if self.minim or minim:
            self.minimPDB()
        if superpose:
            self.superposePDB()
        self.initializeGrid()
        self.genPQR()
        if selstr is not None:
            self.mutatePQR(selstr=selstr)
        if len(self.pdbfiles) == 1 and selstr is None:
            self.mutatePQR()
        if n_workers is None:
            self.genDX_parallel()
        if n_workers is not None:
            self.genDX_parallel(n_workers)
        if esd:
            self.calcESD()
        if esi:
            self.calcESI(idx=idx)
        stop = ti.default_timer()
        print '%s:\tAESOP electrostatic similarity comparison' \
            ' completed in %.2f seconds' % (self.jobname, stop - start)
        warn = self.checkwarnings()
        err = self.checkerrors()
        if warn != 0:
            print 'WARNINGS detected, please view log files!'
        if err != 0:
            print 'ERRORS detected, please view log files!'

    def writeLogs(self, filename=None):
        if filename is None:
            filename = os.path.join(self.jobdir, 'AESOP_logs.txt')
        with open(filename, 'w') as f:
            logs = '\n==== Log Instance ====\n'.join(self.logs)
            f.write(logs)

    def viewLogs(self):
        logs = '\n==== Log Instance ====\n'.join(self.logs)
        print logs

    def checkwarnings(self):
        lines = ''.join(self.logs)
        matches = re.findall('[W][Aa][Rr][Nn][Ii][Nn][Gg]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1

    def checkerrors(self):
        lines = ''.join(self.logs)
        matches = re.findall('[E][Rr][Rr][Oo][Rr]', lines)
        if len(matches) > 0:
            return 1
        elif len(matches) == 0:
            return 0
        else:
            return -1


##########################################################################
# Function to run commands, recording output
##########################################################################


def runProcess(command):
    """Summary
    Simple function intended to capture outputs from processes that write
    to STDOUT.

    Parameters
    ----------
    command : list
        Lists of strings where each element is a part of the entire command.
        Ex: ['script','arg1','arg2',...]

    Returns
    -------
    tuple
        return tuple where first element is output that would have been sent
        to STDOUT and the second element captures errors.
    """

    class runProcess_Exception(Exception):
        pass

    proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    try:
        (out, err) = proc.communicate()
    except:
        raise runProcess_Exception(
            'Unable to execute command - please verify syntax:\n\n\t%s' %
            (command))
        sys.exit(1)
    return (out, err)


##########################################################################
# Function to mutate a single residue in a PDB structure, mutates with
# side-chain truncation
##########################################################################


# Only use this function with PARSE for now ...
def mutatePQR(pqrfile, mutid, resnum, chain=None):
    """Summary
    Mutate PQR file via side-chain truncation scheme (mutate to Alanine)

    Parameters
    ----------
    pqrfile : str
        Full path to PQR file
    mutid : str
        Prefix to use when writing mutated PQR. Should be a full path if
        destination is not in working directory.
    resnum : int
        Residue number to mutate to alanine.
    chain : str, optional
        Chain where residue that will be mutated is located.

    Returns
    -------
    None
        Writes mutated PQR to file specified by the prefix mutid.
    """
    chain = chain.replace(' ', '')
    parent = pd.parsePQR(pqrfile)
    if (chain is None) or (chain is ''):
        residue = parent.select('resnum %d' % (int(resnum)))
        preceed = parent.select('resnum < %d' % (int(resnum)))
        follow = parent.select('resnum > %d' % (int(resnum)))
    elif chain is not None:
        residue = parent.select('chain %s and resnum %d' %
                                (str(chain), int(resnum)))
        preceed = parent.select('chain %s and resnum < %d' %
                                (str(chain), int(resnum)))
        follow = parent.select('chain %s and resnum > %d' %
                               (str(chain), int(resnum)))
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
    ala_atoms = [
        'N', 'H', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O',
        'OXT'
    ]
    if (chain is None) or (chain is ''):
        if preceed is None:
            mutant = residue.select('name ' + ' '.join(ala_atoms)) + follow
        if follow is None:
            mutant = preceed + residue.select('name ' + ' '.join(ala_atoms))
        if (preceed is None) and (follow is None):
            mutant = residue.select('name ' + ' '.join(ala_atoms))
        if (preceed is not None) and (follow is not None):
            mutant = preceed + \
                residue.select('name ' + ' '.join(ala_atoms)) + follow
    else:
        if otherchains is None:
            if preceed is None:
                mutant = residue.select('name ' + ' '.join(ala_atoms)) + follow
            if follow is None:
                mutant = preceed + \
                    residue.select('name ' + ' '.join(ala_atoms))
            if (preceed is None) and (follow is None):
                mutant = residue.select('name ' + ' '.join(ala_atoms))
            if (preceed is not None) and (follow is not None):
                mutant = preceed + \
                    residue.select('name ' + ' '.join(ala_atoms)) + follow
        if otherchains is not None:
            if preceed is None:
                mutant = residue.select('name ' + ' '.join(
                    ala_atoms)) + follow + otherchains
            if follow is None:
                mutant = preceed + \
                    residue.select('name ' + ' '.join(ala_atoms)) + otherchains
            if (preceed is None) and (follow is None):
                mutant = residue.select('name ' + ' '.join(
                    ala_atoms)) + otherchains
            if (preceed is not None) and (follow is not None):
                mutant = preceed + \
                    residue.select('name ' + ' '.join(ala_atoms)
                                   ) + follow + otherchains

    # Write mutant pqr
    pd.writePQR(mutid + '.pqr', mutant)


##########################################################################
# Function to mutate a single residue in a PDB structure, mutates with
# modeller by building internal coordinates of residue
##########################################################################
def complete_structure(pdb, dest=None, disu=False):
    """Summary
    Function to fill in residues with missing atoms. This method simply calls
    complete_pdb from Modeller.

    Parameters
    ----------
    pdb : str
        Full path to pdbfile that will be modified.
    dest : str (optional)
        Full path to destination where completed pdb will be written. If not
        specified, the model object from Modeller will be returned.
    disu : bool (optional)
        If True, complete_pdb will predict and patch all disulfide bridges.
        Default is False.
    """
    class Complete_PDB_Exception(Exception):
        pass

    try:
        from modeller import environ, model, selection
        from modeller.scripts import complete_pdb

        env = environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        mdl = complete_pdb(env, pdb, transfer_res_num=True)
        if disu is True:
            mdl.patch_ss()

        if dest is not None:
            mdl.write(file=dest)
        if dest is None:
            return mdl
    except:
        raise Complete_PDB_Exception('\nUnable to complete missing atoms'
                                     ' in: %s' % (pdb))
        sys.exit(1)


def minimize_cg(struct,
                dest=None,
                disu=True,
                min_atom_shift=0.1,
                max_iter=1000,
                output='NO_REPORT',
                log=None,
                report_iter=10):
    """Summary
    Function to perform conjugate gradient descent minimization in Modeller on
    a user-provided structural file (PDB).

    Parameters
    ----------
    struct : str
        String for path to PDB file
    dest : str
        String for path to location where minimized structure will be written
    disu : bool
        If true, positions of disulfide bridges will be automatically detected
    min_atom_shift : float
        If the max atomic shift between minimization steps is less than this
        value, then convergence is reached and minimization is terminated
    max_iter : int
        Maximum number of calls of objective function before minimization is
        terminated
    output : str
        Valid options are 'NO_REPORT' and 'REPORT'. If set to 'REPORT', then a
        log file during minimation will be printed to screen
    log : str or None
        String for path to location where minimization report will be saved.
        If None, no report will be saved. Report contains only values of
        objective function at after each report interval.
    report_iter : int
        Integer that describes the number of minimization steps to perform
        before reporting the objective function.

    Returns
    -------
    mdl : Model object from Modeller
        If dest is None, the function will return the minimized model.
        If dest is specified, then no model will be returned but the
        minimized model will be written to file.
    """
    from modeller import environ, model, selection
    from modeller.scripts import complete_pdb
    from modeller.optimizers import conjugate_gradients, actions

    class Minimize_CG_Exception(Exception):
        pass

    try:
        if log is not None:
            trcfil = open(log, 'w')

        env = environ()
        env.io.atom_files_directory = ['../atom_files']
        env.edat.dynamic_sphere = True
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        mdl = complete_pdb(env, struct, transfer_res_num=True)

        basename, ext = os.path.splitext(struct)
        if disu is True:
            mdl.patch_ss()

        atmsel = selection(mdl)
        mdl.restraints.make(
            atmsel, restraint_type='STEREO', spline_on_site=False)
        mpdf = atmsel.energy()

        cg = conjugate_gradients(output=output)
        if log is not None:
            cg.optimize(
                atmsel,
                max_iterations=max_iter,
                min_atom_shift=min_atom_shift,
                actions=actions.trace(report, trcfil))
            trcfil.close()
        else:
            cg.optimize(
                atmsel, max_iterations=max_iter, min_atom_shift=min_atom_shift)

        if dest is not None:
            mdl.write(file=dest)
        if dest is None:
            return mdl
    except:
        raise Minimize_CG_Exception('\nCG Minimization failed for: %s' %
                                    (struct))
        sys.exit(1)


def mutatePDB(pdb, mutid, resnum, chain=None, resid='ALA'):
    """Summary
    Function to generate a mutant structure given a local PDB file using
    MODELLER.

    Parameters
    ----------
    pdb : str
        Full path to pdbfile that will be modified.
    mutid : str
        Prefix for mutated structure that will be written. May be a full path
        without file extension if desired output path is not in working
        directory.
    resnum : int, or type that can be forced to int
        Integer number specifying residue number to be mutated.
    chain : str, optional
        Chain ID where specified residue number is to be mutated. This is
        necessary to specify if residue numbers are not unique on each chain.
    resid : str, optional
        Three letter amino acid code specifying the type of mutation. Default
        mutation is to alanine ('ALA').

    Returns
    -------
    None
        Writes mutated structure to file.
    """
    # pdb is the pdb file
    # resnum is the residue number to be mutated
    # chain (optional) can specify what chain the residue to be mutated is
    #                  located on
    # mutid is the prefix for the written mutated structure to be written
    # resid is the residue to mutate to

    try:
        from modeller import environ, model, alignment, selection
    except:
        print('Failed to load modeller: please ensure module is installed'
              ' and license key set')

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
                sel = selection(
                    mdl.residue_range(
                        str(num) + ':' + chid, str(num) + ':' + chid))
            sel.mutate(residue_type=resid)

    aln.append_model(mdl, align_codes='mutant')
    mdl.clear_topology()
    mdl.generate_topology(aln['mutant'])
    mdl.transfer_xyz(aln)

    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.write(file=mutid + '.pdb')

    # Without this section, chainids and resnums from parent won't be retained!
    h = model(env, file=mutid + '.pdb')
    m = model(env, file=pdb)
    aln = alignment(env)
    aln.append_model(m, atom_files=pdb, align_codes='parent')
    aln.append_model(h, atom_files=mutid + '.pdb', align_codes='mutant')
    h.res_num_from(m, aln)  # Restore old residue numbering and chain indexing
    h.write(file=mutid + '.pdb')


def superpose(ref, pdb, atype='CA', output=None):
    """Summary
    Uses Modeller to superpose a PDB file (pdb) to a reference PDB (ref).

    Parameters
    ----------
    ref : str
        Full path to PDB file (or name of file in working directory) that will
        be used as the reference for superpositioning.
    pdb : str
        Full path to PDB file (or name of file in working directory) that will
        be used as the mobile structure for superpositioning.
    atype : str
        Modeller-compatible string selection for atoms to be used in
        superpositioning. We suggest using 'CA'.
    output : str or None
        If output is None, the file specified by pdb will be updated with the
        superposed structure. If specified, output should be a full path where
        the superposed structure will be saved.
    """
    try:
        from modeller import environ, model, alignment, selection
        from modeller.scripts import complete_pdb
    except:
        print('Failed to load modeller: please ensure module is installed and '
              'license key set')

    env = environ()
    env.io.atom_files_directory = '../atom_files'
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    mdl1 = complete_pdb(env, ref, transfer_res_num=True)
    mdl2 = complete_pdb(env, pdb, transfer_res_num=True)

    aln = alignment(env)
    aln.append_model(mdl1, atom_files=ref, align_codes='ref')
    aln.append_model(mdl2, atom_files=pdb, align_codes='pdb')

    aln.malign(gap_penalties_1d=(-600, -400))
    aln.malign3d(
        gap_penalties_3d=(0, 2.0), write_fit=False, write_whole_pdb=False)

    atmsel = selection(mdl1).only_atom_types(atype)
    r = atmsel.superpose(mdl2, aln, superpose_refine=True)

    if output is None:
        output = pdb

    mdl2.write(file=output)


##########################################################################
# Function to run PDB2PQR.exe - should work on any supported OS
##########################################################################


def execPDB2PQR(path_pdb2pqr_exe, pdbfile, outfile=None, ff='parse'):
    """Summary
    Calls the APBS executable according to:
    <path to pdb2pqr appropriate for OS> --ff=parse --chain inputfile outputfile

    Parameters
    ----------
    path_pdb2pqr_exe : str
        Full path to pdb2pqr executable
    pdbfile : str
        PDB file to be converted to a PQR. Should be a full path if not in
        current working directory.
    outfile : str, optional
        File name for PQR file that will be generated. May be a full path if
        desired output is not in current working directory.
    ff : str, optional
        String instructing PDB2PQR what force field to employ. For more
        information visit: http://www.poissonboltzmann.org/docs/pdb2pqr-usage/

    Returns
    -------
    (log, err) : tuple
        When PDB2PQR runs, outputs that would have been sent to STDOUT
        are captured. Log contains the run log and err contains errors.

    Deleted Parameters
    ------------------
    optargs : str, optional
        Description
    """

    class PDB2PQR_Exception(Exception):
        pass

    if outfile is None:
        outfile = os.path.splitext(pdbfile)[0] + '.pqr'
    (log, err) = runProcess([
        path_pdb2pqr_exe, '-v', '--ff=%s' % (ff), '--chain', pdbfile, outfile
    ])
    try:
        pdb = pd.parsePQR(outfile)
    except:
        raise PDB2PQR_Exception(
            '\nPDB2PQR failed for: %s\n\nLogs printed below:\n\n%s' %
            (pdbfile, log))
        sys.exit(1)
    return (log, err)


##########################################################################
# Function to run APBS.exe - should work on any supported OS
##########################################################################


def execAPBS(path_apbs_exe,
             pqr_chain,
             dime,
             glen,
             gcent,
             prefix=None,
             ion=0.150,
             pdie=20.0,
             sdie=78.54,
             dx=False):
    """Summary
    Calls the APBS executable after generating the APBS inputfile.
    Calculates solvation and reference energies.

    Parameters
    ----------
    path_apbs_exe : str
        Full path to APBS executable, EX: 'C:\\APBS\\apbs.exe'.
    pqr_chain : str
        PQR file name containing the segment that will undergo electrostatic
        calculations.
    dime : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    glen : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    gcent : list
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    prefix : str, optional
        Phrase to prepend before any file that is generated before writing.
    ion : float, optional
        Ionic strength for APBS calculation.
    pdie : float, optional
        Protein dielectric constant for APBS calculation.
    sdie : float, optional
        Solvent dielectric constant for APBS calculation.
    dx : bool, optional
        If true, potential files are written.

    Returns
    -------
    file_apbs_log : str
        File name for the log file that APBS generates.
        This file contains results from performed calculations.

    Deleted Parameters
    ------------------
    pqr_complex : STRING
        PQR file name containing the complex that AESOP is analyzing, must
        contain pqr_chain
    grid : float, optional
        Grid spacing for the mesh grid based electrostatic calculations.
        Suggested value of 1 or below
    """

    # path_apbs_exe -   full path to apbs executable ('C:\\APBS\\apbs.exe')
    # pqr_chain     -   path to file with single chain pqr (mutant or parent)
    # pqr_complex   -   path to file with complex pqr, contains mutant or
    #                   parent chain that is used in pqr_chain
    # prefix        -   string to pre-pend to output files (log file, dx file,
    #                   energy file)
    # grid          -   grid spacing using in APBS calculation
    # ion           -   ionic strength for calculation
    # pdie          -   protein dielectric constant
    # sdie          -   solvent dielectric constant
    class APBS_Exception(Exception):
        pass

    if prefix is None:
        prefix = os.path.splitext(pqr_chain)[0]

    # Format APBS input file
    cmd_read = [
        'read\n',
        '   mol pqr %s\n' % (pqr_chain),
        'end\n'
    ]
    cmd_solv = [
        'elec name solv\n', '   mg-manual\n',
        '   dime %d %d %d\n' % (dime[0], dime[1], dime[2]),
        '   glen %d %d %d\n' % (glen[0], glen[1], glen[2]),
        '   gcent %d %d %d\n' % (gcent[0], gcent[1], gcent[2]), '   mol 1\n',
        '   lpbe\n', '   bcfl sdh\n', '   srfm smol\n', '   chgm spl2\n',
        '   ion 1 %.2f 2.0\n' % (ion), '   ion -1 %.2f 2.0\n' % (ion),
        '   pdie %.2f\n' % (pdie), '   sdie %.2f\n' % (sdie),
        '   sdens 10.0\n', '   srad 0.0\n', '   swin 0.3\n',
        '   temp 298.15\n', '   calcenergy total\n'
    ]
    if dx is True:
        cmd_solv = cmd_solv + ['   write pot dx %s\n' % (prefix)]
    cmd_solv = cmd_solv + ['end\n']
    cmd_ref = [
        'elec name ref\n', '   mg-manual\n', '   dime %d %d %d\n' %
        (dime[0], dime[1], dime[2]), '   glen %d %d %d\n' %
        (glen[0], glen[1], glen[2]), '   gcent %d %d %d\n' %
        (gcent[0], gcent[1], gcent[2]), '   mol 1\n', '   lpbe\n',
        '   bcfl sdh\n', '   srfm smol\n', '   chgm spl2\n', '   pdie %.2f\n' %
        (pdie), '   sdie %.2f\n' % (pdie), '   sdens 10.0\n', '   srad 0.0\n',
        '   swin 0.3\n', '   temp 298.15\n', '   calcenergy total\n', 'end\n'
    ]
    cmd_write = [
        'print elecEnergy solv end\n', 'print elecEnergy ref end\n', 'quit\n'
    ]
    apbs_in = cmd_read + cmd_solv + cmd_ref + cmd_write

    # Write APBS input file
    file_apbs_in = prefix + '.in'
    file_apbs_log = prefix + '.log'
    with open(file_apbs_in, 'w') as f:
        for line in apbs_in:
            f.write(line)

    # Execute APBS
    (log, err) = runProcess([
        path_apbs_exe, '--output-file=%s' % (file_apbs_log),
        '--output-format=flat', file_apbs_in
    ])
    pattern = re.compile(
        '(?<=Global net ELEC energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
    )
    elec = np.asarray([x.split()
                       for x in re.findall(pattern, log)]).astype(np.float)
    elec = elec.reshape((1, elec.size))
    if len(elec[0]) != 2:
        # print '\nAPBS failed for: %s' % (file_apbs_in)
        perror = re.compile('[E][Rr][Rr][Oo][Rr]')
        status = 0
        apbs_compiled_error_list = '\n'
        for l in log.split('\n'):
            m = re.findall(perror, l)
            if len(m) > 0:
                status = 1
            if status == 1:
                apbs_compiled_error_list = apbs_compiled_error_list + l + '\n'

        raise APBS_Exception(
            '\nAPBS failed for: %s\n\nLogs printed below:\n\n%s' % (
                file_apbs_in, apbs_compiled_error_list))
        sys.exit(1)

    # return file_apbs_log
    return (elec, log)


##########################################################################
# Function to run APBS.exe to generate a DX file only - should work on any
# supported OS
##########################################################################


def calcDX(path_apbs_exe,
           pqrfile,
           prefix=None,
           grid=1.0,
           ion=0.150,
           pdie=20.0,
           sdie=78.54,
           cfac=1.5,
           glen=None,
           gcent=np.zeros(3),
           dime=None):
    """Summary
    Calls the APBS executable after generating the APBS inputfile. Generates
    a potential file (DX).

    Parameters
    ----------
    path_apbs_exe : str
        Full path to APBS executable, EX: 'C:\\APBS\\apbs.exe'.
    pqrfile : str
        The PQR file that will be used to generate a grid of electrostatic
        potentials. Must be a full path if file is not in current path.
    prefix : str, optional
        Phrase to prepend before any file that is generated before writing.
    grid : float, optional
        Distance spacing of grid points. If the grid dimensions are not
        divisible by three, resolution will be increased (smaller grid
        spacing) until grid dimensions are divisible by three.
    ion : float, optional
        Ionic strength for APBS calculation.
    pdie : float, optional
        Protein dielectric constant for APBS calculation.
    sdie : float, optional
        Solvent dielectric constant for APBS calculation.
    cfac : float, optional
        Scaling factor for grid dimensions. We suggest to leave this unchanged.
    glen : None, optional
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    gcent : TYPE, optional
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/
    dime : None, optional
        List of three integers. Parameter required for APBS. Please see
        description at: http://www.poissonboltzmann.org/docs/apbs-overview/

    Returns
    -------
    (log, err) : tuple
        When APBS runs, outputs that would have been sent to STDOUT are
        captured. Log contains the run log and err contains errors.
    """
    if prefix is None:
        prefix = os.path.splitext(pqrfile)[0]

    # Determine mesh dimensions according to Ron's AESOP protocol in the R
    # source file
    if (dime is None) | (glen is None):
        pqr = pd.parsePQR(pqrfile)
        coords = pqr.getCoords()
        x = coords[:, 0]
        y = coords[:, 1]
        z = coords[:, 2]
        fg = np.array(
            (np.ceil(np.max(x) - np.min(x)), np.ceil(np.max(y) - np.min(y)),
             np.ceil(np.max(z) - np.min(z))))
        fg = np.ceil((fg + 5) * cfac)
        glen = np.zeros(3)
        glen = np.vstack((glen, fg)).max(axis=0)
        dime_list = (32 * np.linspace(1, 100, 100)) + \
            1  # list of possible dime values
        # index of dime to use from list, subtract one to be consistent with
        # python indexing!
        dime_ind = np.ceil(fg / (32 * grid)) - 1
        dime = np.array(
            (dime_list[int(dime_ind[0])], dime_list[int(dime_ind[1])],
             dime_list[int(dime_ind[2])]))
        ix = 0
        iy = 0
        iz = 0
        counter = 0
        while ((dime[0] * dime[1] * dime[2] % 3 != 0) or (counter <= 5)):
            ix += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iy += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            iz += 1
            if (dime[0] * dime[1] * dime[2] % 3 != 0):
                dime = np.array((dime_list[int(dime_ind[0] + ix)],
                                 dime_list[int(dime_ind[1]) + iy],
                                 dime_list[int(dime_ind[2]) + iz]))
            counter += 1
        dime = dime  # .reshape((1, 3))
        glen = glen  # .reshape((1, 3))
        gcent = pd.calcCenter(pqr).astype(int)

    # Format APBS input file
    cmd_read = ['read\n', '   mol pqr %s\n' % (pqrfile), 'end\n']
    cmd_solv = [
        'elec name solv\n', '   mg-manual\n',
        '   dime %d %d %d\n' % (dime[0], dime[1], dime[2]),
        '   glen %d %d %d\n' % (glen[0], glen[1], glen[2]),
        '   gcent %d %d %d\n' % (gcent[0], gcent[1], gcent[2]), '   mol 1\n',
        '   lpbe\n', '   bcfl sdh\n', '   srfm smol\n', '   chgm spl2\n',
        '   ion 1 %.2f 2.0\n' % (ion), '   ion -1 %.2f 2.0\n' % (ion),
        '   pdie %.2f\n' % (pdie), '   sdie %.2f\n' % (sdie),
        '   sdens 10.0\n', '   srad 0.0\n', '   swin 0.3\n',
        '   temp 298.15\n', '   write pot dx %s\n' % (prefix), 'end\n'
    ]

    cmd_write = ['quit\n']
    apbs_in = cmd_read + cmd_solv + cmd_write

    # Write APBS input file
    file_apbs_in = prefix + '.in'
    file_apbs_log = prefix + '.log'
    with open(file_apbs_in, 'w') as f:
        for line in apbs_in:
            f.write(line)

    # Execute APBS
    (log, err) = runProcess([
        path_apbs_exe, '--output-file=%s' % (file_apbs_log),
        '--output-format=flat', file_apbs_in
    ])
    return log


##########################################################################
# Function to run multiple APBS processes at once
##########################################################################


def batchAPBS(kernel):
    """Summary
    Function required to run multiple APBS jobs simultaneously. Not intended
    for general use.

    Parameters
    ----------
    kernel : tuple
        Tuple of parameters required for APBS.

    Returns
    -------
    ndarray
        i, j represent the index in the matrix with which the calculated
        energies correspond. The last two elements are the solvation and
        reference energies, respectively.
    """
    path, pqr_chain, dime, glen, gcent, prefix, ion, pdie, sdie, dx, i, j = kernel
    energies, log = execAPBS(
        path,
        pqr_chain,
        dime,
        glen,
        gcent,
        prefix=prefix,
        ion=ion,
        pdie=pdie,
        sdie=sdie,
        dx=dx)
    return (np.array([i, j, energies[0][0], energies[0][1]]), log)


##########################################################################
# Function to run multiple APBS processes at once for the purpose of
# generating only a DX file
##########################################################################


def batchCalcDX(kernel):
    """Summary
    Function required to run multiple APBS jobs simultaneously. Not
    intended for general use.

    Parameters
    ----------
    kernel : tuple
        Tuple of parameters required for APBS.

    Returns
    -------
    None
        Writes files according to calcDX function.
    """
    path, pqrfile, prefix, grid, ion, pdie, sdie, cfac, glen, gcent, dime = kernel
    log = calcDX(
        path,
        pqrfile,
        prefix=prefix,
        grid=grid,
        ion=ion,
        pdie=pdie,
        sdie=sdie,
        cfac=cfac,
        glen=glen,
        gcent=gcent,
        dime=dime)
    return log


##########################################################################
# Function to run multiple Coulomb processes at once
##########################################################################


def batchCoulomb(kernel):
    """Summary
    Function required to run multiple Coulomb jobs simultaneously. Not
    intended for general use.

    Parameters
    ----------
    kernel : tuple
        Tuple of parameters required for APBS.

    Returns
    -------
    ndarray
        i, j represent the index in the matrix with which the calculated
        energies correspond. The last element is the Coulombic energy.
    """
    path, pqr_chain, pdie, i, j = kernel
    energies, log = execCoulomb(path, pqr_chain)
    energies = energies / pdie
    return (np.array([i, j, energies]), log)


##########################################################################
# Function to run coulomb.exe - should work on any supported OS
##########################################################################


def execCoulomb(path_coulomb_exe, pqr):
    """Summary
    Call Coulomb from APBS tools to calculate Coulombic energies.

    Parameters
    ----------
    path_coulomb_exe : str
        Full path to coulomb executable.
    pqr : TYPE
        Filename for PQR to use for Coulombic energy calculation. Must be
        full path if not in current path.

    Returns
    -------
    float
        Coulombic energy associated with input PQR file.
    """
    (log, err) = runProcess([path_coulomb_exe, pqr])
    pattern = re.compile(
        '(?<=Total energy =)\s+[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'
    )  # May need to update regex
    coul = np.asarray(re.findall(pattern, log)).astype(np.float)
    return (coul, log)


##########################################################################
# Function to run DSSP.exe - should work on any supported OS
##########################################################################
# def execDSSP(pdbfile, dssp):
#     """Summary

#     Parameters
#     ----------
#     pdbfile : TYPE
#         Description
#     dssp : TYPE
#         Description

#     Returns
#     -------
#     TYPE
#         Description
#     """
#     (log, err) = runProcess([dssp, pdbfile])
#     return log

##########################################################################
# Function to plot results of Alascan
##########################################################################


def writePDB(alascan, filename=None):
    """Summary
    Function to write free energies of association/solvation into B-factor
    column of PDB for easy visualization of results.

    Parameters
    ----------
    alascan : scan class
        Alascan or DirectedMutagenesis class after running the complete
        analysis.
    filename : str, optional
        Full path to file where PDB file will be written. Defaults to job
        directory.
    """
    jobdir = alascan.jobdir
    pdbfile = os.path.join(jobdir, alascan.pdb_complex_dir, 'wt.pdb')
    resnums = [item for sublist in alascan.list_resnums
               for item in sublist][1:]
    chids = [item for sublist in alascan.list_chids for item in sublist][1:]
    if isinstance(resnums[0], list):
        resnums = [val for sublist in resnums for val in sublist]
    if isinstance(chids[0], list):
        chids = [val for sublist in chids for val in sublist]
    if len(alascan.selstr) > 1:
        ddG = alascan.ddGa_rel()[1:]
    else:
        ddG = alascan.dGsolv_rel()[1:]
    pdb = pd.parsePDB(pdbfile)
    pdb.setBetas(np.zeros((pdb.numAtoms())))
    pdb.setOccupancies(np.zeros((pdb.numAtoms())))
    for chid, resnum, g in zip(chids, resnums, ddG):
        atom = pdb.select('chain %s and resnum %s' % (chid, resnum))
        atom.setBetas(g)
    if filename is None:
        filename = os.path.join(jobdir, 'wt.ddG.pdb')
    pd.writePDB(filename, pdb)


def plotNetwork(scan,
                filename=None,
                title='',
                dpi=300,
                cutoff=5.,
                E=2.5,
                node_size=1500,
                font_size=12,
                alpha=0.8,
                edge_color='g',
                edge_width=3.,
                layout=None,
                **kwargs):
    """Summary
    Function to visualize electrostatic interactions from a scan class
    (Alascan or Directed Mutagenesis). Requires networkx to be installed.

    Parameters
    ----------
    scan : Alascan or DirectedMutagenesis class
        Scan class where calculation of free energies is complete.
    filename : str or None
        Full path to file where figure will be saved. If None, no figure
        is saved, but the plot is displayed and the graph is returned.
    title : str
        Matplotlib style title for plot.
    dpi : int
        Integer specifying the dots per inch, or image resolution.
    cutoff : float
        Distance cutoff in Angstroms for determining if a electrostatic
        interaction occurs. Default value is 5 Angstroms.
    E : float
        Threshold for determing those nodes that should be included in
        the network based on the value of the free energy perturbation
        that results from mutating the amino acid. If the magnitude of
        the free energy of association relative to the parent structure
        is greater than E, then the node is included. Default is
        2.5 kJ/mol.
    node_size : int
        Parameter to scale size of nodes in network. Larger values
        result in nodes with larger diameter.
    font_size : int
        Font size for text within network. 12 pt font is default.
    alpha : float
        Set transparency of nodes. Default is 0.8. Accepted range
        is [0, 1].
    edge_color : str
        Matplotlib-style specification of line color. Default is 'g'
        (green).
    edge_width : int
        Set the line width for edges. Default is 3 pt font.
    layout : Networkx layout kernel or None
        Network layout from networkx. Extra arguments for this layout
        may be passed as key word arguments to plotNetwork.
    """

    import networkx as nx
    from scipy.spatial.distance import pdist, squareform

    selstr = scan.selstr
    jobdir = scan.jobdir

    pdbfile = os.path.join(scan.jobdir,
                           scan.pqr_complex_dir,
                           'wt.pqr')

    # Match mutids to residues
    if len(selstr) > 1:
        ddGa = scan.ddGa_rel()[1:]
    else:
        ddGa = scan.dGsolv_rel()[1:]
    mutids = scan.getMutids()[1:]

    mutids = np.asarray(mutids)
    mutids = mutids[np.where(np.abs(ddGa) >= E)[0]]
    mutids = mutids.tolist()
    mutids = [w.replace('sel', 'seg') for w in mutids]
    ddGa = np.asarray(ddGa)
    ddGa = ddGa[np.where(np.abs(ddGa) >= E)[0]]
    ddGa = ddGa.tolist()

    n = len(mutids)
    seg2sel = dict(('seg%d' % (key + 1), val)
                   for (key, val) in enumerate(selstr))

    regions = [seg2sel[x.split('_')[0]] for x in mutids]
    numbers = [int(x.split('_')[1][1:-1]) for x in mutids]
    pdb = pd.parsePQR(pdbfile)

    atoms = pdb.select(
        ' and '.join([
            regions[0],
            'resnum %s' % (numbers[0]),
            '(charge >= 0.3 or charge <= -0.3)',
            'sidechain'  # ,
            # 'heavy'
        ])
    )
    chains = [atoms.getChids()[0]]
    resids = [AA_dict[atoms.getResnames()[0]]]

    for region, number in zip(regions[1:], numbers[1:]):
        currsel = pdb.select(
            ' and '.join([
                region,
                'resnum %s' % (number),
                '(charge >= 0.3 or charge <= -0.3)',
                'sidechain'  # ,
                # 'heavy'
            ])
        )
        atoms = atoms + currsel
        chains.append(currsel.getChids()[0])
        resids.append(AA_dict[currsel.getResnames()[0]])

    # Calc distance matrix
    resnums = np.asarray(numbers)
    chains = np.asarray(chains)
    resids = np.asarray(resids)
    lbls = ['%s%d%s' % (rid, res, chid)
            for rid, res, chid in zip(resids, resnums, chains)]
    n = len(resnums)

    atom_dist = squareform(pdist(atoms.getCoords()))
    res_dist = -1.0 * np.ones((n, n))

    for i, i_res, i_chid in zip(range(n), resnums, chains):

        i_idx = np.where(np.logical_and(
            atoms.getResnums() == i_res,
            atoms.getChids() == i_chid
        )
        )[0]

        for j, j_res, j_chid in zip(range(n), resnums, chains):
            if i != j:
                j_idx = np.where(np.logical_and(
                    atoms.getResnums() == j_res,
                    atoms.getChids() == j_chid
                )
                )[0]
                res_dist[i, j] = np.min(atom_dist[np.ix_(i_idx, j_idx)])
            else:
                res_dist[i, j] = 0

    # Find edges
    res_dist *= np.tri(*res_dist.shape, k=-1)
    mask = np.logical_and((res_dist > 0), (res_dist <= cutoff))
    row, col = np.where(mask)

    # Generate graph
    vmax = np.max(np.abs(ddGa))
    norm = mpl.colors.Normalize(vmin=-1 * vmax, vmax=vmax)
    cmap = plt.cm.coolwarm_r

    G = nx.Graph()
    for node, val, lbl in zip(range(n), ddGa, lbls):
        G.add_node(node, ddGa=val, lbl=lbl, alpha=alpha)  # , size=node_size)
    for i, j in zip(row, col):
        G.add_edge(i, j, weight=res_dist[i, j])

    if layout is None:
        pos = nx.shell_layout(G)
    else:
        pos = layout(G, **kwargs)

    energy = nx.get_node_attributes(G, 'ddGa')
    labels = nx.get_node_attributes(G, 'lbl')
    # nodesize = nx.get_node_attributes(G, 'size')
    # d      = nx.degree(G)

    # Draw Network
    plt.figure()
    nx.draw_networkx(G,
                     pos=pos,
                     cmap=cmap,
                     norm=norm,
                     alpha=alpha,
                     node_color=ddGa,
                     vmin=-1.0 * vmax,
                     vmax=vmax,
                     edge_color=edge_color,
                     width=edge_width,
                     labels=labels,
                     font_size=font_size,
                     node_size=node_size)  # [float(v) * float(node_size) for v in d.values()])

    plt.title(title)
    plt.axis('off')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cb = plt.colorbar(sm)
    if len(selstr) > 1:
        cb.set_label('$\Delta\Delta G_A (kJ/mol)$', size=font_size)
    else:
        cb.set_label('$\Delta G_S (kJ/mol)$', size=font_size)

    if filename is None:
        plt.show()
        return G
    else:
        plt.savefig(filename, dpi=dpi)


def plotNetwork_interactive(scan,
                            filename=None,
                            title='',
                            cutoff=5.,
                            E=2.5,
                            font_size=14,
                            node_size=20,
                            edge_color='#888',
                            edge_width=0.5,
                            display_output='external',
                            layout=None,
                            **kwargs):
    """Summary
    Function to visualize electrostatic interactions from a scan class
    (Alascan or Directed Mutagenesis). Figure is more interactive than
    the standard matplotlib figure. Requires networkx to be installed.

    Parameters
    ----------
    scan : Alascan or DirectedMutagenesis class
        Scan class where calculation of free energies is complete.
    filename : str or None
        Full path to file where figure will be saved. If None, no figure
        is saved, but the plot is displayed and the graph is returned.
    title : str
        Matplotlib style title for plot.
    cutoff : float
        Distance cutoff in Angstroms for determining if a electrostatic
        interaction occurs. Default value is 5 Angstroms.
    E : float
        Threshold for determing those nodes that should be included in
        the network based on the value of the free energy perturbation
        that results from mutating the amino acid. If the magnitude of
        the free energy of association relative to the parent structure
        is greater than E, then the node is included. Default is
        2.5 kJ/mol.
    node_size : int
        Parameter to scale size of nodes in network. Larger values
        result in nodes with larger diameter.
    font_size : int
        Font size for text within network. 12 pt font is default.
    edge_color : str
        Matplotlib-style specification of line color. Default is 'g'
        (green).
    edge_width : int
        Set the line width for edges. Default is 3 pt font.
    display_output : str
        Set output to either open local html file in browser or inline
        plot in notebook.
    layout : Networkx layout kernel or None
        Network layout from networkx. Extra arguments for this layout
        may be passed as key word arguments to plotNetwork.
    """
    import networkx as nx
    from scipy.spatial.distance import pdist, squareform
    import plotly.plotly as py
    import plotly
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
    from plotly import tools

    selstr = scan.selstr
    jobdir = scan.jobdir

    pdbfile = os.path.join(scan.jobdir,
                           scan.pqr_complex_dir,
                           'wt.pqr')

    # Match mutids to residues
    if len(selstr) > 1:
        ddGa = scan.ddGa_rel()[1:]
        colorbartitle = 'ddG (kJ/mol)'
    else:
        ddGa = scan.dGsolv_rel()[1:]
        colorbartitle = 'dG (kJ/mol)'
    mutids = scan.getMutids()[1:]

    mutids = np.asarray(mutids)
    mutids = mutids[np.where(np.abs(ddGa) >= E)[0]]
    mutids = mutids.tolist()
    mutids = [w.replace('sel', 'seg') for w in mutids]
    ddGa = np.asarray(ddGa)
    ddGa = ddGa[np.where(np.abs(ddGa) >= E)[0]]
    ddGa = ddGa.tolist()

    n = len(mutids)
    seg2sel = dict(('seg%d' % (key + 1), val)
                   for (key, val) in enumerate(selstr))

    regions = [seg2sel[x.split('_')[0]] for x in mutids]
    numbers = [int(x.split('_')[1][1:-1]) for x in mutids]
    pdb = pd.parsePQR(pdbfile)

    atoms = pdb.select(
        ' and '.join([
            regions[0],
            'resnum %s' % (numbers[0]),
            '(charge >= 0.3 or charge <= -0.3)',
            'sidechain'  # ,
            # 'heavy'
        ])
    )
    chains = [atoms.getChids()[0]]
    resids = [AA_dict[atoms.getResnames()[0]]]

    for region, number in zip(regions[1:], numbers[1:]):
        currsel = pdb.select(
            ' and '.join([
                region,
                'resnum %s' % (number),
                '(charge >= 0.3 or charge <= -0.3)',
                'sidechain'  # ,
                # 'heavy'
            ])
        )
        atoms = atoms + currsel
        chains.append(currsel.getChids()[0])
        resids.append(AA_dict[currsel.getResnames()[0]])

    # Calc distance matrix
    resnums = np.asarray(numbers)
    chains = np.asarray(chains)
    resids = np.asarray(resids)
    lbls = ['%s%d%s' % (rid, res, chid)
            for rid, res, chid in zip(resids, resnums, chains)]
    n = len(resnums)

    atom_dist = squareform(pdist(atoms.getCoords()))
    res_dist = -1.0 * np.ones((n, n))

    for i, i_res, i_chid in zip(range(n), resnums, chains):

        i_idx = np.where(np.logical_and(
            atoms.getResnums() == i_res,
            atoms.getChids() == i_chid
        )
        )[0]

        for j, j_res, j_chid in zip(range(n), resnums, chains):
            if i != j:
                j_idx = np.where(np.logical_and(
                    atoms.getResnums() == j_res,
                    atoms.getChids() == j_chid
                )
                )[0]
                res_dist[i, j] = np.min(atom_dist[np.ix_(i_idx, j_idx)])
            else:
                res_dist[i, j] = 0

    # Find edges
    res_dist *= np.tri(*res_dist.shape, k=-1)
    mask = np.logical_and((res_dist > 0), (res_dist <= cutoff))
    row, col = np.where(mask)

    # Generate graph
    vmax = np.max(np.abs(ddGa))
    G = nx.Graph()
    for node, val, lbl in zip(range(n), ddGa, lbls):
        G.add_node(node, ddGa=val, lbl=lbl)
    for i, j in zip(row, col):
        G.add_edge(i, j, weight=res_dist[i, j])
    if layout is None:
        pos = nx.spring_layout(G)
    else:
        pos = layout(G, **kwargs)
    dmin = 1
    ncenter = 0
    for n in pos:
        x, y = pos[n]
        d = (x - 0.5)**2 + (y - 0.5)**2
        if d < dmin:
            ncenter = n
            dmin = d

    p = nx.single_source_shortest_path_length(G, ncenter)
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=go.Line(width=edge_width, color=edge_color),
        hoverinfo='none',
        mode='lines')

    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace['x'] += [x0, x1, None]
        edge_trace['y'] += [y0, y1, None]

    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        hoverinfo='text',
        marker=go.Marker(
            showscale=True,
            # colorscale options
            # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
            # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' |
            # 'YIGnBu'
            colorscale='RdBu',
            cmin=-1 * vmax,
            cmax=vmax,
            reversescale=True,
            color=[],
            size=node_size,
            colorbar=dict(
                thickness=15,
                title=colorbartitle,
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    for node in G.nodes():
        x, y = pos[node]
        node_trace['x'].append(x)
        node_trace['y'].append(y)

    energy = nx.get_node_attributes(G, 'ddGa')
    for node in G.nodes():
        node_trace['marker']['color'].append(energy[node])
        node_info = "{0:.2f}".format(G.node[node]['ddGa']) + ' kJ/mol'
        node_trace['text'].append(node_info)

    network_annotations = go.Annotations()
    for k in range(len(lbls)):
        network_annotations.append(
            go.Annotation(
                text=lbls[k],
                x=pos[k][0], y=pos[k][1],
                xref='x1', yref='y1',
                font=dict(size=font_size),
                arrowsize=.4,
                showarrow=True)
        )

    network_fig = go.Figure(data=go.Data([edge_trace, node_trace]),
                            layout=go.Layout(
        title='<br>Network Plot',
        titlefont=dict(size=16),
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=40),
        annotations=network_annotations,
        xaxis=go.XAxis(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=go.YAxis(showgrid=False, zeroline=False, showticklabels=False)))
    if display_output == 'external':
        plotly.offline.plot(network_fig)
    elif display_output == 'notebook':
        init_notebook_mode(connected=True)
        plotly.offline.iplot(network_fig)
    if filename is not None:
        py.image.save_as(network_fig, filename=filename)

def plotScan(Alascan, filename=None):
    """Summary
    Function to display results from the computational alanine or directed
    mutagenesis scan.

    Parameters
    ----------
    Alascan : scan class
        Alascan or DirectedMutagenesis class after running the complete
        analysis.
    filename : None, optional
        If the resulting plot should be written to disk, specify a filename.
        Otherwise, the image will only be saved.

    Returns
    -------
    tuple
        Handles to generated figure.
    """
    # plt.style.use('seaborn-talk')
    if filename is not None:
        plt.switch_backend('agg')
    elif os.environ.get('DISPLAY', '') == '':
        print('No display variable found. Supply a filename to generate plot '
              'using non-interactive Agg backend')
    figure, axarr = plt.subplots(len(Alascan.mutid) - 1, sharey=True)
    dpi_val = 300
    if len(Alascan.mutid) > 2:
        for i in xrange(1, len(Alascan.mutid)):
            if len(Alascan.mutid[i]) is not 0:
                axarr[i - 1].set_title(
                    np.unique(
                        np.array([w.split('_') for w in Alascan.mutid[i]
                                  ])[:, 0])[0] + ' ddGa relative to WT')
                axarr[i - 1].set_ylabel('kJ/mol')
                axarr[i - 1].set_xticks(
                    np.arange(
                        len(Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]])))
                if 100 < len(Alascan.mutid[i]) <= 150:
                    axarr[i - 1].set_xticklabels(
                        np.array(
                            [w.split('_') for w in Alascan.mutid[i]])[:, 1],
                        rotation='vertical',
                        ha='left',
                        size=6)
                elif len(Alascan.mutid[i]) > 150:
                    axarr[i - 1].set_xticklabels(
                        np.array(
                            [w.split('_') for w in Alascan.mutid[i]])[:, 1],
                        rotation='vertical',
                        ha='left',
                        size=2)
                    dpi_val = 600
                else:
                    axarr[i - 1].set_xticklabels(
                        np.array(
                            [w.split('_') for w in Alascan.mutid[i]])[:, 1],
                        rotation='vertical',
                        ha='left')
                axarr[i - 1].bar(
                    np.arange(
                        len(Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]]))[
                            Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] > 0],
                    Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]][
                        Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] > 0],
                    color='red')
                axarr[i - 1].bar(
                    np.arange(
                        len(Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]]))[
                            Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] < 0],
                    Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]][
                        Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] < 0],
                    color='blue')
                axarr[i - 1].xaxis.set_ticks_position('bottom')
                axarr[i - 1].yaxis.set_ticks_position('left')
    elif len(Alascan.mutid) == 2:
        axarr.set_title(
            np.unique(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 0])[0] +
            ' dGsolv relative to WT')
        axarr.set_ylabel('kJ/mol')
        axarr.set_xticks(
            np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]])))
        if 100 < len(Alascan.mutid[1]) <= 150:
            axarr.set_xticklabels(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1],
                rotation='vertical',
                ha='left',
                size=6)
        elif len(Alascan.mutid[1]) > 150:
            axarr.set_xticklabels(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1],
                rotation='vertical',
                ha='left',
                size=2)
            dpi_val = 600
        else:
            axarr.set_xticklabels(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1],
                rotation='vertical',
                ha='left')
        axarr.bar(
            np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))[
                Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] > 0],
            Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][Alascan.dGsolv_rel(
            )[Alascan.mask_by_sel[:, 1]] > 0],
            color='red')
        axarr.bar(
            np.arange(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))[
                Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] < 0],
            Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][Alascan.dGsolv_rel(
            )[Alascan.mask_by_sel[:, 1]] < 0],
            color='blue')
        axarr.xaxis.set_ticks_position('bottom')
        axarr.yaxis.set_ticks_position('left')
    plt.tight_layout()
    if filename is not None:
        figure.savefig(filename, dpi=dpi_val)
    return (figure, axarr)


def plotScan_interactive(Alascan, display_output = 'external', filename=None):
    """Summary
    Function to display results from the computational alanine or directed
    mutagenesis scan. Figure is more interactive than the standard matplotlib
    figure.

    Parameters
    ----------
    Alascan : scan class
        Alascan or DirectedMutagenesis class after running the complete
        analysis.
    display_output : str
        Set output to either open local html file in browser or inline
        plot in notebook.
    filename : None, optional
        If the resulting plot should be written to disk, specify a filename.
        Otherwise, the image will only be saved.

    Returns
    -------
    None
        Saves image of figure, if desired.
    """
    import plotly.plotly as py
    import plotly
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
    from plotly import tools

    subplot_titles = []
    for i in range(1, len(Alascan.mutid)):
        subplot_titles.append(
            np.unique(
                np.array([w.split('_') for w in Alascan.mutid[i]])[:, 0])[0] +
            ' ddGa relative to WT')
    fig = tools.make_subplots(
        rows=len(Alascan.mutid) - 1, cols=1, subplot_titles=subplot_titles)
    if len(Alascan.mutid) > 2:
        for i in range(len(Alascan.mutid) - 1, 0, -1):
            pos_y = np.zeros(
                len(Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]]))
            pos_y[Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] >
                  0] = Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]][
                      Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] > 0]
            neg_y = np.zeros(
                len(Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]]))
            neg_y[Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] <
                  0] = Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]][
                      Alascan.ddGa_rel()[Alascan.mask_by_sel[:, i]] < 0]
            pos_trace = go.Bar(
                x=np.array([w.split('_') for w in Alascan.mutid[i]])[:, 1],
                y=pos_y,
                name=np.unique(
                    np.array([w.split('_') for w in Alascan.mutid[i]
                              ])[:, 0])[0] + 'Loss of binding',
                marker=dict(color='rgba(0,136,55,1)'))
            neg_trace = go.Bar(
                x=np.array([w.split('_') for w in Alascan.mutid[i]])[:, 1],
                y=neg_y,
                name=np.unique(
                    np.array([w.split('_') for w in Alascan.mutid[i]
                              ])[:, 0])[0] + 'Gain in binding',
                marker=dict(color='rgba(123,50,148,1)'))
            fig.append_trace(pos_trace, i, 1)
            fig.append_trace(neg_trace, i, 1)
            fig['layout']['yaxis' + str(i)].update(title='kJ/mol')

    if len(Alascan.mutid) == 2:
        pos_y = np.zeros(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))
        pos_y[Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] >
              0] = Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][
                  Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] > 0]
        neg_y = np.zeros(len(Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]]))
        neg_y[Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] <
              0] = Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]][
                  Alascan.dGsolv_rel()[Alascan.mask_by_sel[:, 1]] < 0]
        pos_trace = go.Bar(
            x=np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1],
            y=pos_y,
            name=np.unique(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 0])[0] +
            'Loss of binding',
            marker=dict(color='rgba(0,136,55,1)'))
        neg_trace = go.Bar(
            x=np.array([w.split('_') for w in Alascan.mutid[1]])[:, 1],
            y=neg_y,
            name=np.unique(
                np.array([w.split('_') for w in Alascan.mutid[1]])[:, 0])[0] +
            'Gain in binding',
            marker=dict(color='rgba(123,50,148,1)'))
        fig.append_trace(pos_trace, 1, 1)
        fig.append_trace(neg_trace, 1, 1)
        fig['layout']['yaxis' + str(1)].update(title='kJ/mol')

    fig['layout'].update(barmode='stack', hovermode='closest')
    plotly_fig = go.Figure(fig)
    if display_output == 'external':
        plotly.offline.plot(plotly_fig)
    elif display_output == 'notebook':
        init_notebook_mode(connected=True)
        plotly.offline.iplot(plotly_fig)
    if filename is not None:
        py.image.save_as(plotly_fig, filename=filename)


##########################################################################
# Function to plot results of ESD.calc()
##########################################################################


def plotESD(esd, filename=None, cmap='hot'):
    """Summary
    Function to display an electrostatic similarity heatmap from a previously
    run ElecSimilarity class.

    Parameters
    ----------
    esd : ndarray
        ESD matrix from ElecSimilarity class (ElecSimilarity.esd).
    filename : str, optional
        If the resulting plot should be written to disk, specify a filename.
        Otherwise, the image will only be saved.
    cmap : str, optional
        Colormap from matplotlib to use.

    Returns
    -------
    None
        Writes image to disk, if desired.
    """
    # plt.style.use('seaborn-talk')
    fig, ax = plt.subplots(sharey=True)
    heatmap = ax.pcolor(esd.esd, cmap=cmap, vmin=0, vmax=2)
    ax.set_xlim(0, esd.esd.shape[0])
    ax.set_ylim(0, esd.esd.shape[1])
    ax.set_xticks(np.arange(esd.esd.shape[0]) + 0.5, minor=False)
    ax.set_yticks(np.arange(esd.esd.shape[1]) + 0.5, minor=False)
    ax.set_xticklabels(esd.ids, rotation=90)
    ax.set_yticklabels(esd.ids)
    fig.colorbar(heatmap)
    plt.tight_layout()
    if filename is not None:
        fig.savefig(filename)


##########################################################################
# Function to plot ESD dendrogram
##########################################################################


def plotDend(esd, filename=None):
    """Summary
    Function to display an electrostatic similarity dendrogram from a
    previously run ElecSimilarity class.

    Parameters
    ----------
    esd : ElecSimilarity class
        ElecSimilarity class containing final esd matrix.
    filename : str, optional
        If the resulting plot should be written to disk, specify a filename.
        Otherwise, the image will only be saved.

    Returns
    -------
    None
        Writes image to disk, if desired.
    """
    # plt.style.use('seaborn-talk')
    fig, ax = plt.subplots(sharey=True)
    Z = cluster.linkage(esd.esd)
    cluster.dendrogram(
        Z,
        labels=esd.ids,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        ax=ax)
    plt.xlabel('Variants')
    plt.ylabel('ESD')
    plt.tight_layout()
    if filename is not None:
        fig.savefig(filename)


def plotESD_interactive(esd, filename=None, cmap='YIGnBu', display_output='external'):
    """Summary
    Function to display an electrostatic similarity heatmap from a previously
    run ElecSimilarity class. Figure is more interactive that the standard
    matplotlib figure.

    Parameters
    ----------
    esd : ElecSimilarity class
        ElecSimilarity class containing final esd matrix.
    filename : str, optional
        If the resulting plot should be written to disk, specify a filename.
        Otherwise, the image will only be saved.
    cmap : str, optional
        Colormap from matplotlib to use.
    display_output : str
        Set output to either open local html file in browser or inline
        plot in notebook.

    Returns
    -------
    None
        Writes image to disk, if desired.
    """
    import plotly.plotly as py
    import plotly
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
    from plotly.tools import FigureFactory as FF

    figure = FF.create_dendrogram(
        esd.esd, orientation='bottom', labels=esd.ids)
    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'
    # Create Side Dendrogram
    dendro_side = FF.create_dendrogram(
        esd.esd, orientation='right', labels=esd.ids)
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
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = go.Data([
        go.Heatmap(
            x=dendro_leaves, y=dendro_leaves, z=heat_data, colorscale=cmap)
    ])
    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = figure['layout']['xaxis']['tickvals']

    # Add Heatmap Data to Figure
    figure['data'].extend(go.Data(heatmap))

    # Edit Layout
    figure['layout'].update({
        'showlegend': False,
        'hovermode': 'closest',
    })
    figure['layout'].update({'margin': {'b': 140, 't': 10}})

    # Edit xaxis
    figure['layout']['xaxis'].update({
        'domain': [.15, 1],
        'mirror': False,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'ticks': ""
    })
    # Edit xaxis2
    figure['layout'].update({
        'xaxis2': {
            'domain': [0, .15],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    })

    # Edit yaxis
    figure['layout']['yaxis'].update({
        'domain': [0, .85],
        'mirror': False,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'showticklabels': False,
        'side': 'right',
        'ticktext': dendro_side['layout']['yaxis']['ticktext'],
        'tickvals': dendro_side['layout']['yaxis']['tickvals'],
        'ticks': ""
    })
    # Edit yaxis2
    figure['layout'].update({
        'yaxis2': {
            'domain': [.825, .975],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }
    })

    # Plot!
    if display_output == 'external':
        plotly.offline.plot(figure)
    elif display_output == 'notebook':
        init_notebook_mode(connected=True)
        plotly.offline.iplot(figure)
    if filename is not None:
        py.image.save_as(figure, filename=filename)


# ######################################################################################################################################################
# # Function to calculate RSA for PDB
# ######################################################################################################################################################
# def calcRSA(pdbfile, dssp):
#     """Summary

#     Parameters
#     ----------
#     pdbfile : TYPE
#         Description
#     dssp : TYPE
#         Description

#     Returns
#     -------
#     TYPE
#         Description
#     """
#     # SA from: C. Chotia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14.
#     SA_dict = {'C': 135, 'D': 150, 'S': 115, 'Q': 180, 'K': 200,
#                'I': 175, 'P': 145, 'T': 140, 'F': 210, 'N': 160,
#                'G': 75, 'H': 195, 'L': 170, 'R': 225, 'W': 255,
#                'A': 115, 'V': 155, 'E': 190, 'Y': 230, 'M': 185}
#     threshold = 0.2 # RSA below this value is considered completely buried
#     log = execDSSP(pdbfile, dssp)
#     ag = pd.parsePDB(pdbfile)
#     n_atoms = ag.numAtoms()
#     ACC = np.zeros(n_atoms, float)
#     lines = filter(None, log.split('\n'))
#     iterator = iter(lines)
#     for line in iterator:
#         # print line
#         if line.startswith('  #  RESIDUE'):
#             break
#     for line in iterator:
#         if line[13] == '!':
#             continue
#         res = ag[(line[11], int(line[5:10]), line[10].strip())]
#         if res is None:
#             continue
#         indices = res.getIndices()
#         ACC[indices] = int(line[35:38])
#         ag.setData('dssp_acc', ACC)

#     resid = np.asarray([AA_dict[x] for x in ag.getResnames()])
#     sasa = ag._getData('dssp_acc')
#     rsa =[]
#     for res, sa in zip(resid, sasa):
#         rsa.append(sa / SA_dict[res])
#         # rsa = np.asarray([SA_dict[res] for res, sasa in zip(resid, sasa])
#     rsa = np.asarray(rsa)
#     return(resid, rsa, rsa>=threshold)

##########################################################################
# Dictionary to convert between 3 letter and 1 letter amino acid codes
##########################################################################
"""Summary

Attributes
----------
AA_dict : dict
    Variable used to convert between 1 letter and 3 letter amino acid IDs
"""
AA_dict = {
    'CYS': 'C',
    'ASP': 'D',
    'SER': 'S',
    'GLN': 'Q',
    'LYS': 'K',
    'ILE': 'I',
    'PRO': 'P',
    'THR': 'T',
    'PHE': 'F',
    'ASN': 'N',
    'GLY': 'G',
    'HIS': 'H',
    'LEU': 'L',
    'ARG': 'R',
    'TRP': 'W',
    'ALA': 'A',
    'VAL': 'V',
    'GLU': 'E',
    'TYR': 'Y',
    'MET': 'M'
}
