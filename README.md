# aesop-python
Rewriting AESOP in python.

## Framework
- **Inputs**
	- PDB, selections, ionic strength, dielectric constants
- **Generate PQR of parent**
	- PDB2PQR
- **Generate Mutant PQR**
	- Modeller
		- Use model.mutate()
		- Clean up extra files
	- PDB2PQR
- **Batch APBS**
- **ESD**
- **Report results**

## Classes
**Alascan**
- input pdb
- selstr = vector of selctions
- ionstr
- pdie
- sdie
- dirs
- prefix

### Functions
- {filename} = genPQR(pdb_in, pqr_out, arg) *where arg is --chain, --ff, etc.*
- {filename} = genMut(pdb_in, selstr, dir )
- {id_energy} = batchAPBS(pqrs, args) *where arg is apbs args*
