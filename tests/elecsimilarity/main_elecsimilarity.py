from aesop import ElecSimilarity, plotDend, plotESD
import pickle as p

import glob as glob
import os as os

path_apbs = 'C:\\APBS\\apbs.exe'
path_coul = 'C:\\APBS\\coulomb.exe'
path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr-windows-bin64-2.0.0\\pdb2pqr.exe'

pdb_dir = 'LTP_pdb'
jobname = 'elecsimilarity'
pdbfiles = glob.glob(os.path.join(pdb_dir, '*.pdb'))

family = ElecSimilarity(pdbfiles, path_pdb2pqr, path_apbs, jobname=jobname)
family.centerPDB()
family.superposePDB()
family.initializeGrid()
family.genPQR()
family.genDX()
family.calcESD()

p.dump(family, open(jobname+'.p', 'wb'))
family = p.load(open(jobname+'.p', 'rb'))

plotDend(family, filename=jobname+'_dendrogram.png')
plotESD(family, filename=jobname+'_esd.png')