.. _pdbpreparation::

Preparing PDB Files
===================

AESOP requires protein structures that comply with the PDB format. Given a structure from the 
`Protein Data Bank <www.rcsb.org>`_, the user needs to consider that coordinates for some atoms 
may not be resolved in the deposited structure. Thus, residues may be missing from a protein 
even though the sequence may be known. To fix such issues, the user must perform homology 
modelling to model and refine gaps. A number of computational tools exist to add these missing 
residues including Modeller, PDBFixer, UCSF Chimera, Pymol, and `SWISS-MODEL <https://swissmodel.expasy.org/>`_. 
If you are willing to install OpenMM and all other libraries within the `OMNIA <www.omnia.md>`_ channel, 
then `PDBFixer <https://github.com/pandegroup/pdbfixer>`_ may be installed through anaconda quite easily 
and offers a simple graphical user interface for PDB preparation.

If the PDB contains all residues but is missing a few atoms in one or more residues, AESOP has a 
function that will call complete_pdb from Modeller to fill in missing atoms. You may use the function 
as follows::

	from aesop import complete_structure
	
	pdbfile = 'input.pdb'
	outfile = 'output.pdb'
	
	complete_structure(pdb=pdbfile, dest=outfile, disu=False)
	
The above snipped of code reads ``input.pdb``, fills in missing atoms (not missing residues), and 
writes the completed structure to ``output.pdb``. If ``input.pdb`` has disulfide bridges, simply 
set ``disu=True`` to predict and patch disulfides.

While AESOP will handle protein structures where residue numbering overlaps between protein chains, 
we advise users to make sure only a single model is present in the PDB file to prevent unforseen 
complications. Additionally, each chain should be represented by a unique identifier that complies 
with the PDB format.
