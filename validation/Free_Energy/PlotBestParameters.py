# GOMC Example for the Gibbs Ensemble (GEMC) using MoSDeF [1, 2, 5-10, 13-17]

# Note: In this specific example, we will be using the GEMC_NVT ensemble.


# Import the required packages and specify the force field (FF) being used. 

# Note: For GOMC, the residue names are treated as molecules, so the residue names must be unique for each different molecule. [1, 2, 13-17]

# Note: Each residue can be set to a different FF, which is done by setting the residue name to a FF in a dictionary (FF_Dict).  The FF selection can be a FF name (set from foyer FF repositor) or a specified FF xml file. [1, 2, 13-17]

import shutil
import pathlib
import random
from pathlib import Path
import re
import os
from os.path import exists
from Surface import find_minimum
import pickle as pickle
p = re.compile("Wolf_Calibration_(\w+?)_(\w+?)_BOX_(\d+)_(\w+?).dat")


box = "0"

bestValueFileName = "bestWolfParameters"
if (exists(bestValueFileName+".pickle")):
    with open(bestValueFileName+".pickle", 'rb') as handle:
          model2BestWolfAlphaRCut = pickle.load(handle)

    WolfMethods = ["VLUGT", "GROSS", "VLUGTWINTRACUTOFF"]
    WolfPotentials = ["DSF", "DSP"]
    WolfModels = []
    box = str(0)
    import itertools
    for root, dirs, files in os.walk(".", topdown=False):
       for direc in dirs:
          if ("Calibration" == direc):
             model = os.path.basename(root)
             print("Model ", model)
             if (model not  in WolfModels):
                WolfModels.append(model)
    for model in WolfModels:
        for element in itertools.product(WolfMethods, WolfPotentials):
              for element in itertools.product(WolfMethods, WolfPotentials):
                    print(model, element[0], element[1], box)
                    print(model2BestWolfAlphaRCut[model][(element[0], element[1], box)])
                    alphaRcutRelErrTuple = model2BestWolfAlphaRCut[model][(element[0], element[1], box)]
                    print("Cutoff" , alphaRcutRelErrTuple[0])
                    print("Alpha" , alphaRcutRelErrTuple[1])
    """
                    with open(path2File, "a") as myfile:
                        defPotLine = "Wolf\tTrue\t{pot}\n".format(pot=WolfDefaultPotential)
                        myfile.write(defPotLine)
                        defKindLine = "WolfKind\t{kind}\n".format(kind=WolfDefaultKind)
                        myfile.write(defKindLine)
    """


