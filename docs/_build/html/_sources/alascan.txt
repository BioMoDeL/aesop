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