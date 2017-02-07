"""
ALANINE SCAN
	Example case: Barnase-Barstar
"""

from aesop import Alascan, plotScan, writePDB
try:
	from aesop import plotNetwork
except:
	print 'Unable to import plotNetwork, is the NetworkX library installed?'

if __name__ == '__main__': # Protect entry point of Windows application for parallel processes

	path_apbs    = 'apbs'		# Make sure these
	path_coulomb = 'coulomb'	#	are correct 
	path_pdb2pqr = 'pdb2pqr'	#	for your sytem

	jobname = 'alascan'
	pdbfile = 'barnase_barstar.pdb'

	selstr = ['chain A', 'chain B']

	alascan = Alascan(pdb=pdbfile, 
					  pdb2pqr_exe=path_pdb2pqr,
					  apbs_exe=path_apbs,
					  coulomb_exe=path_coulomb,
					  jobname=jobname, 
					  selstr=selstr, 
					  minim=False)
	# alascan = Alascan(pdb=pdbfile, jobname=jobname, selstr=selstr)

	alascan.run()
	# alascan.run_parallel(6) # Uncomment to run in parallel on 6 threads
	# alascan.run_parallel()  # Uncomment to run on half of available threads

	alascan.viewLogs()
	alascan.writeLogs(filename="alascan_logs.txt")

	plotScan(alascan, filename='alascan.png')

	try:
		plotNetwork(alascan, filename='network.png')
	except:
		print 'Skipping plotNetwork example!'

	mut_ids  = alascan.getMutids()
	energies = alascan.ddGa_rel()

	alascan.summary(filename='alascan_summary.txt')

	writePDB(alascan, filename='alascan.ddGa.pdb')
