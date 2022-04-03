# GOMC Example for the Gibbs Ensemble (GEMC) using MoSDeF [1, 2, 5-10, 13-17]

# Note: In this specific example, we will be using the GEMC_NVT ensemble.


# Import the required packages and specify the force field (FF) being used. 

# Note: For GOMC, the residue names are treated as molecules, so the residue names must be unique for each different molecule. [1, 2, 13-17]

# Note: Each residue can be set to a different FF, which is done by setting the residue name to a FF in a dictionary (FF_Dict).  The FF selection can be a FF name (set from foyer FF repositor) or a specified FF xml file. [1, 2, 13-17]

import shutil
import pathlib
import random
from pathlib import Path

import os

WolfDefaultKind = "Cassandra"
WolfDefaultPotential = "DSF"
WolfDefaultAlpha = [0.21]

WolfCutoffBoxList = [0]

WolfCutoffCoulombLowerBoundList = [10]
WolfCutoffCoulombUpperBoundList = [15]
WolfCutoffCoulombIntervalList = [0.5]

WolfAlphaLowerBoundList = [0.0]
WolfAlphabUpperBoundList = [0.5]
WolfAlphaIntervalList = [0.01]

for root, dirs, files in os.walk(".", topdown=False):
   for name in files:
      if(name == "NVT_Cal_water_ethanol_fe.conf"):
         path2File = os.path.join(root, name)
         with open(path2File, "a") as myfile:
            defPotLine = "Wolf\tTrue\t{pot}\n".format(pot=WolfDefaultPotential)
            myfile.write(defPotLine)
            defKindLine = "WolfKind\t{kind}\n".format(kind=WolfDefaultKind)
            myfile.write(defKindLine)
            myfile.write("WolfCalibration\tTrue\n")
            for box, wolfAlphaLower, wolfAlphaUpper, wolfAlphaInterval, wolfCutoffLower, wolfCutoffUpper, wolfCutoffInterval, defaultAlpha \
            in zip(WolfCutoffBoxList, WolfCutoffCoulombLowerBoundList, WolfCutoffCoulombUpperBoundList, WolfCutoffCoulombIntervalList, \
            WolfAlphaLowerBoundList, WolfAlphabUpperBoundList, WolfAlphaIntervalList, WolfDefaultAlpha):
               defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=defaultAlpha)
               myfile.write(defAlphaLine)

               CutoffLine = "WolfCutoffCoulombRange\t{box}\t{lb}\t{ub}\t{inter}\n".format(box=box, lb=wolfCutoffLower, ub=wolfCutoffUpper, inter=wolfCutoffInterval)
               myfile.write(CutoffLine)

               alphaLine = "WolfAlphaRange\t{box}\t{lb}\t{ub}\t{inter}\n".format(box=box, lb=wolfAlphaLower, ub=wolfAlphaUpper, inter=wolfAlphaInterval)
               myfile.write(alphaLine)

