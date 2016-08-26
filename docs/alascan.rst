.. _alascan::

Alanine Scan
============

Given a PDB structure of atomistic resolution, the alanine scan method iteratively perturbs the native 
structure by mutating single amino acids to alanine one residue at a time. In this manner, the method 
can predict those mutations that are predicted to significantly affect the free energy of association for 
a complex according to the thermodynamic cycle.

.. currentmodule:: aesop
.. autosummary::
    :toctree: api/generated/
    
    Alascan


Example case: Barnase-Barstar
"""""""""""""""""""""""""""""

Open a new python session, import the Alascan class, and import the plotScan function:: 

    from aesop import Alascan, plotScan

Next, you must specify the full paths to your ``apbs``, ``coulomb``, and ``pdb2pqr`` executables, if 
the paths for the directories containing the executables have not already been added to the environment. 
Here is an example for a Windows system::

    path_apbs    = 'C:\\APBS\\apbs.exe'
    path_coulomb = 'C:\\APBS\\coulomb.exe'
    path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr.exe'

Next we will specify the jobname and pdbfile to used in the method. After running the alanine scan, jobname 
will be used to create a folder where files for the method will be generated. You can download the PDB file 
for this example from this link (:download:`download<../data/barnase_barstar.pdb>`). 
Make sure you place the PDB in your working directory::

    jobname = 'alascan'
    pdbfile = 'barnase_barstar.pdb'

The Alascan class will need to know how to define each subunit of the protein complex. To do this, 
the user should specify a list of selection strings. Each element of the list should be a stand-alone 
selection string that fully describes how to select the associated subunit. If the selection string 
list has only 1 element, then be aware that you may only calculate solvation free energies as no 
association of subunits occurs. In this example, barnase is chain A, and barstar is chain B. Thus, 
we specify the selection string in the following manner::

    selstr = ['chain A', 'chain B']

Finally, we may initialize the Alanine scan class::

    alascan = Alascan(pdbfile=pdbfile, pdb2pqr_exe=path_pdb2pqr,
                      apbs_exe=path_apbs, coulomb_exe=path_coulomb,
                      jobname=jobname, selstr=selstr)

After initialization, you can run the analysis in series::

    alascan.run()

... or you can run the analysis in parallel on a certain number of threads::

    alascan.run_parallel(6)

Once complete, you can view the results as a barplot::

    plotScan(alascan)

You may also extract the free energies of association and the associated mutation ids::

    mut_ids  = alascan.getMutids()
    energies = alascan.ddGa_rel()

If you wish, you can use build in function to summarize results. If the file name is not 
specified for the summary, then the summary is simply printed to STDOUT::

    alascan.summary(filename='alascan_summary.txt')


References
""""""""""

Cite