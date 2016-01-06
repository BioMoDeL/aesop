
import os as os
import sys as sys
import datetime as dt
import numpy as np
import prody as pd

######################################################################################################################################################
# Container for performing an Alanine Scan with AESOP
#   pdb     -   PDB file for performing Alascan. Must contain all chain selections with standard aminoacid nomenclature
#   selstr  -   List of selection strings for performing mutations
#   ion     -   Ionic strength
#   pdie    -   Protein dielectric constant   
#   sdie    -   Solvent dielectric constant
######################################################################################################################################################
class Alascan:
    def __init__(self, pdb, selstr=['Protein'], ion=0.150, pdie=20.0, sdie=78.54):
        self.pdb = pdb
        self.selstr = selstr
        self.ion = ion
        self.pdie = pdie
        self.sdie = sdie
        # Insert code to instantiate dirs and prefix
        self.dirs = 0
        self.prefix = '%4d%02d%02d'%(dt.date.today().year, dt.date.today().month, dt.date.today().day)

    def getPDB(self):
        return self.pdb
    def getSel(self):
        return self.selstr
    def getDirs(self):
        return self.dirs
    def getPrefix(self):
        return self.prefix
    def getEnergies():
        return 0
    def getMutID():
        return 0

    def genPQR():
        0
    def genMut():
        0
    def batchAPBS():
        0
    def run():
        0
