.. _elecsimilarity::

Electrostatic Similarity
========================

The electrostatic similarity method generates grid potentials for a list of PDB files and compares 
all potential files in a pairwise manner. Here we will provide a test case that compares several 
members of a family of plant proteins. This example is based on a more comprehensive, published study [Chae2010]_.

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
(:download:`download<../data/LTP_pdb.zip>`). After downloading the PDB files, unzip them and place 
them in the current working directory::

    pdbfiles = ['1MZL.pdb', 'SCA1.pdb', 'SCA3.pdb']

When the method is run, intermediate files will be generated and stored in a folder of the current 
working directory. The user has the option of naming this folder by specifying a job name::

    jobname = 'LTP_test'

Next, the method is initialized by::

    family = ElecSimilarity(pdbfiles, pdb2pqr_exe=path_pdb2pqr, apbs_exe=path_apbs, 
                            jobname=jobname)

Finally, we are ready to run the analysis. This will take several minutes, so please be patient::

    family.run()

You can view results using built-in functions::

    plotDend(family)
    plotESD(family)

... or you can access the raw ESD matrix::

    data = family.esd

Other modules such as numpy or pandas will allow exporting of the ESD matrix to file::

    import numpy as np
    np.savetxt(fname='esd_matrix.txt', data, fmt='%.4f')

References
""""""""""

.. [Chae2010] `Chae, K., B.J. Gonong, S.C. Kim, C.A. Kieslich, D. Morikis, S. Balasubramanian, and E.M. Lord. 2010. A multifaceted study of stigma/style cysteine-rich adhesin (SCA)-like Arabidopsis lipid transfer proteins (LTPs) suggests diversified roles for these LTPs in plant growth and reproduction. J. Exp. Bot. 61: 4277–4290 <https://doi.org/10.1093/jxb/erq228>`_.