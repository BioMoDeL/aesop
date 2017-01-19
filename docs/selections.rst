.. _selections::

Atomic Selections
=================

In order to allow for advanced atomic selections with protein structures, AESOP utilizes selection 
macros from ProDy. Using these macros, the end user can easily string together boolean statements based 
on protein chains, residue numbers, amino acid, atom name, physicochemical properties, or even distance 
criteria that describe the portion of the protein structure file to select. For more examples and 
explanations concerning selection strings please see the `ProDy webpage for atomic selections 
<http://prody.csb.pitt.edu/manual/reference/atomic/select.html>`_.

Basic Examples
""""""""""""""

Selection string for chain A of a PDB file::

    'chain A'
	
Selection string for chain A and residue numbers 1 to 100::

    'chain A and resnum 1 to 100'
	
Selection string for chain A and protein::

    'chain A and protein'

Selection for chain A or chain C::

    'chain A or chain C'
	
Alanine scan example
""""""""""""""""""""

For the alanine scan, we suggest each element of the ``selstr`` list contains a separate chain::

    selstr = ['chain A', 'chain B', 'chain C']
	
Ideally, all chains used should comprise a single protein complex. In case the protein complex is 
quite large and only a handful of mutations are of interest to the end-user, then the region argument 
may be used to restrict ``selstr`` to some subset of residues. For this reason, ``region`` should be a list 
of the same length as ``selstr``. In the current example, if we only wish to mutate residues within 10 
angstroms of chain C, then you could specify ``region`` in the following manner::

    region = ['within 10 of chain C', 'within 10 of chain C', 'within 10 of chain C']

Once again, for in depth discussion of more complicated selection strings, please refer to the `ProDy 
website <http://prody.csb.pitt.edu/manual/reference/atomic/select.html>`_.

DirectedMutagenesis scan example
""""""""""""""""""""""""""""""""

For the directed mutagenesis scan, the user must specify the subunits of the protein complex (``selstr``), 
the targeted residue(s) to mutate (``target``), and mutation to perform (``mutation``). As in the alanine scan, 
we suggest each element of ``selstr`` to contain a separate chain of the protein complex::

    selstr = ['chain A', 'chain B']
	
In order to specify targeted residues to mutate, each element of ``target`` must contain all residues 
that will be mutated. Since residue numbers may overlap between chains of the protein structre, the user 
may need to additionally specify a chain. For example::

    target = ['resnum 50', 'resnum 50 in chain B', '(resnum 50 or resnum 60) and chain B']
	
In the first element ('resnum 50), all residues with residue number 50 will be mutated. In the second element 
('resnum 50 in chain B'), only the residue with number 50 in chain B will be mutated. In the third element 
('(resnum 50 or resnum 60) and chain B'), only residues numbered 50 or 60 in chain B will be mutated. 

Next, the user must specify how to mutate each element of ``target`` by specifying a 3 letter amino acid code for 
each element of the target. These codes should be stored in a list (here, we use the variable name ``mutation``)::

    mutation = ['ALA', 'ARG', 'ASP']
	
Since each element of ``target`` corresponds to each element of ``mutation``, the mutations specified above will perform several different 
mutations. Namely, residues selected by the first element of ``target`` will be mutated to alanine; residues selected by the 
second element of ``target`` will be mutated to arginine; and residues selected by the third element of ``target`` will be mutated 
to aspartatic acid. Currently AESOP does not support mutation of two amino acids to two different amino acids simultaneously, 
though this may be added as a feature in the future. In general, we prefer to mutate one amino acid at a time to prevent 
significantly changing the structure of the native protein throughout the analysis. 
