"""
ELECTROSTATIC SIMILARITY
	Example case 1: LTP plant proteins
"""

import numpy as np
from aesop import ElecSimilarity, plotDend, plotESD

path_apbs    = 'apbs'			# Make sure these are correct
path_pdb2pqr = 'pdb2pqr'   		# 	for your system

pdbfiles = ['1MZL.pdb', 'SCA1.pdb', 'SCA3.pdb']
jobname = 'LTP_test1'

family = ElecSimilarity(pdbfiles=pdbfiles, 
						pdb2pqr_exe=path_pdb2pqr, 
						apbs_exe=path_apbs,
						jobname=jobname)

family.run(superpose=True, center=False)
# family.run_parallel(superpose=True, center=False) # Uncomment for parallel run

family.viewLogs()
family.writeLogs(filename="family_logs.txt")

plotDend(family, filename='dend.png')
plotESD(family, filename='esd.png')

data = family.esd
np.savetxt('esd_matrix.txt', data, fmt='%.4f')

family.calcESI()

# family.run(esi=True, esd=False, superpose=True) 			# Commented out to 
# family.run_parallel(esi=True, esd=False, superpose=True)	# 	prevent extra runs