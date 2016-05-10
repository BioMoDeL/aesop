from aesop import Alascan, plotScan
import pickle as p

path_apbs = 'C:\\APBS\\apbs.exe'
path_coul = 'C:\\APBS\\coulomb.exe'
path_pdb2pqr = 'C:\\PDB2PQR\\pdb2pqr-windows-bin64-2.0.0\\pdb2pqr.exe'
path_dssp = 'C:\\Python27\\Scripts\\dssp-2.0.4-win32.exe'

jobname = 'alascan'
pdbfile = 'barnase_barstar.pdb'

selstr = ['chain A', 'chain B']

alascan = Alascan(pdbfile, path_pdb2pqr, path_apbs, coulomb_exe=path_coul, selstr=selstr, jobname=jobname,
                      region=None, grid=1, ion=0.150, pdie=20.0, sdie=78.54, ff='parse', cfac=1.5, dx=False)

if __name__=="__main__":
    alascan.run_parallel()

p.dump(alascan, open(jobname+'.p', 'wb'))

plotScan(alascan, filename=jobname+'.png')