"""
DIRECTED MUTAGENESIS
	Example case: Barnase-Barstar
"""

from aesop import Alascan, plotScan, writePDB

path_apbs    = 'apbs'		# Make sure these
path_coulomb = 'coulomb'	#	are correct 
path_pdb2pqr = 'pdb2pqr'	#	for your sytem

jobname = 'directedscan'
pdbfile = 'barnase_barstar.pdb'

selstr = ['chain A', 'chain B']

target = ['resnum 27',  'resnum 73',  'resnum 83',  'resnum 87',  # mutations in chain A
          'resnum 145', 'resnum 149', 'resnum 164', 'resnum 186'] # mutations in chain B
mutation = ['ASP', 'LYS', 'GLU', 'GLU', # mutations in chain A
            'ARG', 'ARG', 'ASP', 'LYS'] # mutations in chain B

mutscan = DirectedMutagenesis(pdb=pdbfile, 
							  pdb2pqr_exe=path_pdb2pqr,
							  apbs_exe=path_apbs, 
							  coulomb_exe=path_coulomb,
							  jobname=jobname, 
							  selstr=selstr, 
							  target=target,
							  mutation=mutation,
							  minim=True)

mutscan.run()
# mutscan.run_parallel(6) # Uncomment to run in parallel on 6 threads
# mutscan.run_parallel()  # Uncomment to run on half of available threads

mutscan.viewLogs()
mutscan.writeLogs(filename="mutscan_logs.txt")

plotScan(mutscan, filename='directedmutagenesis.png')

mut_ids  = mutscan.getMutids()
energies = mutscan.ddGa_rel()

mutscan.summary(filename='mutscan_summary.txt')

writePDB(mutscan, filename='mutscan.ddGa.pdb')