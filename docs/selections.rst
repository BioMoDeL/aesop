.. _selections::

Atomic Selections
=================

In order to allow for advanced atomic selections with protein structures, AESOP utilizes selection 
macros from ProDy. Using these macros, the end user can easily string together boolean statements based 
on protein chains, residue numbers, amino acid, atom name, physicochemical properties, or even distance 
criteria that describe the portion of the protein structure file to select. For more examples and 
explanations concerning selection strings please see the `ProDy webpage for atomic selections 
<http://prody.csb.pitt.edu/manual/reference/atomic/select.html>`_.

Examples
"""""""""

Selection string for chain A of a PDB file::

    'chain A'
	
Selection string for chain A and residue numbers 1 to 100::

    'chain A and resnum 1 to 100'
	
Selection string for chain A and protein::

    'chain A and protein'

Selection for chain A or chain C::

    'chain A or chain C'
