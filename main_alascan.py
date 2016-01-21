from aesop_alascan import *
import pickle as p

path_apbs = 'C:\\APBS\\apbs.exe'
path_coul = 'C:\\APBS\\coulomb.exe'
path_pdb2pqr = 'C:\\pdb2pqr-2.1.0\\pdb2pqr.exe'

pdbfile = 'barnase_barstar.pdb'

alascan = Alascan(pdbfile, path_pdb2pqr, path_apbs, coulomb_exe=path_coul, selstr=['chain A', 'chain B'], 
	jobname='test1', region=None, grid=1, ion=0.150, pdie=20.0, sdie=78.54)

alascan.run()

p.dump(alascan, open('bn_bs_alascan.p', 'wb'))