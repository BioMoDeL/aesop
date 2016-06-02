from aesop import DirectedMutagenesis, plotScan
import pickle as p

path_apbs = 'C:\\APBS\\apbs.exe'
path_coul = 'C:\\APBS\\coulomb.exe'
path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr-windows-bin64-2.0.0\\pdb2pqr.exe'

jobname = 'mutscan'
pdbfile = 'barnase_barstar.pdb'

selstr = ['chain A', 'chain B']

target = ['resnum 27', 'resnum 27', 'resnum 73', 'resnum 73', 'resnum 142', 'resnum 142', 'resnum 145', 'resnum 145']
mutation = ['ALA', 'ASP', 'ALA', 'LYS', 'ALA', 'LYS', 'ALA', 'LYS']

mutscan = DirectedMutagenesis(pdbfile, target, mutation, path_pdb2pqr, path_apbs, coulomb_exe=path_coul,
                                  selstr=selstr, jobname=jobname,grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse',
                                  cfac=1.5, dx=False)
mutscan.run_parallel()
p.dump(mutscan, open(jobname+'.p', 'wb'))
plotScan(mutscan, filename=jobname+'.png')