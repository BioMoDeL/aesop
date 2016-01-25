import aesop_alascan as ala
import pickle as p
import prody as pd
reload(ala)

path_apbs = 'C:\\APBS\\apbs.exe'
path_coul = 'C:\\APBS\\coulomb.exe'
path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr-windows-bin64-2.0.0\\pdb2pqr.exe'

jobname = '3OXU_B_D_truncation'
pdbfile = '3OXU_B_D.pdb'

alascan = ala.Alascan(pdbfile, path_pdb2pqr, path_apbs, coulomb_exe=path_coul, selstr=['chain B', 'chain D'],
	jobname=jobname, region=None, grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse')

alascan.run_truncated()

p.dump(alascan, open(jobname+'_alascan.p', 'wb'))

# ala.plotResults(alascan, filename=jobname+'.png')









# alascan = p.load(open(jobname+'_alascan.p', 'rb'))
# import numpy as np
# np.savetxt('mutids.txt',alascan.getMutids(), fmt='%s')
# np.savetxt('ddGbind_rel_parse.txt', alascan.ddGbind_rel())