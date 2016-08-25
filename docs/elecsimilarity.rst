.. _elecsimilarity::

Electrostatic Similarity
========================

The electrostatic similarity method generates grid potentials for a list of PDB files and compares 
all potential files in a pairwise manner. Here we will provide a test case that compares several 
members of a family of plant proteins. PDB files for the analysis may be downloaded from Github 
by `clicking here <https://raw.githubusercontent.com/rohithmohan/aesop-python/master/tests/elecsimilarity/LTP_pdb.zip>`_.

.. currentmodule:: aesop
.. autosummary::
	:toctree: api/generated/
	
	ElecSimilarity
	
Example case: LTP plant proteins
""""""""""""""""""""""""""""""""

Open a new python session, import the ElecSimilarity class, and import the plotDend and plotESD functions:: 

    from aesop import ElecSimilarity, plotDend, plotESD

Next, you must specify the full paths to your ``apbs`` and ``pdb2pqr`` executables, if 
the paths for the directories containing the executables have not already been added to the environment. 
Here is an example for a Windows system::

    path_apbs    = 'C:\\APBS\\apbs.exe'
    path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr.exe'

Now we will specify what PDB files the method should compare. Here we will use only 3 PDB files 
(`download <https://raw.githubusercontent.com/rohithmohan/aesop-python/master/tests/elecsimilarity/LTP_pdb.zip>`_)::

    pdbfiles = ['1MZL.pdb', 'SCA1.pdb', 'SCA3.pdb']
	
When the method is run, intermediate files will be generated and stored in a folder of the current 
working directory. The user has the option of naming this folder by specifying a job name::

    jobname = 'LTP_test'

Next, the method is initialized by::
	
    family = ElecSimilarity(pdbfiles, pdb2pqr_exe=path_pdb2pqr, apbs_exe=path_apbs, jobname=jobname)

Finally, we are ready to run the analysis. This will take several minutes, so please be patient::

    family.run()
	
You can view results using built-in functions::

    plotDend(family)
	plotESD(family)
	
... or you can access the raw ESD matrix::

    data = family.esd
