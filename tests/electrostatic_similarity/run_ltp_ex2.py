"""
ELECTROSTATIC SIMILARITY
	Example case 2: Alascan of a LTP plant protein
"""

from aesop import ElecSimilarity, plotDend, plotESD

if __name__ == '__main__': # Protect entry point of Windows application for parallel processes

	path_apbs    = 'apbs'		# Make sure these are correct
	path_pdb2pqr = 'pdb2pqr'	# 	for your system

	pdbfiles     = ['1MZL.pdb']
	jobname      = 'LTP_test2'

	family = ElecSimilarity(pdbfiles=pdbfiles, 
							pdb2pqr_exe=path_pdb2pqr, 
							apbs_exe=path_apbs, 
							jobname=jobname,
							minim=False)

	family.run(superpose=False, 
			   esd=True, 
			   esi=True, 
			   selstr=['protein'])