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




bestValueFileName = "bestWolfParameters"
if (not exists(bestValueFileName+".pickle")):
    model2BestWolfAlphaRCut = dict()
    for root, dirs, files in os.walk(".", topdown=False):
       for direc in dirs:
          if ("Calibration" == direc):
             model = os.path.basename(root)
             print("Model ", model)
             if (model not  in model2BestWolfAlphaRCut):
                model2BestWolfAlphaRCut[model] = dict()
    for root, dirs, files in os.walk(".", topdown=False):
       for name in files:
          if(p.match(name)):
             groups = p.search(name)
             wolfKind = groups.group(1)
             potential = groups.group(2)
             box = groups.group(3)
             head_tail = root.split(os.sep)
             for direc in head_tail:
                 for key in model2BestWolfAlphaRCut.keys():
                     if(key == direc):
                        print ("model" , key)
                        print ("wolf Kind" , wolfKind)
                        print ("potential Kind" , potential)
                        print ("box" , box)

                        tupleMin = find_minimum(os.path.join(root, name), key, wolfKind, potential, box, True)
                        # Use smaller error, either BF or Grad Desc
                        model2BestWolfAlphaRCut[key][(wolfKind, potential, box)] = [tupleMin[3], tupleMin[4], tupleMin[5], tupleMin[6], tupleMin[7], tupleMin[8]]
             print(model2BestWolfAlphaRCut)


       with open(bestValueFileName+".pickle", 'wb') as handle:
          pickle.dump(model2BestWolfAlphaRCut, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
   with open(bestValueFileName+".pickle", 'rb') as handle:
      model2BestWolfAlphaRCut = pickle.load(handle)

WolfMethods = ["VLUGT", "GROSS", "VLUGTWINTRACUTOFF"]
WolfPotentials = ["DSF", "DSP"]
box = str(0)
import itertools

for root, dirs, files in os.walk(".", topdown=False):
   for name in files:
      if(name == "NVT_Prod_water_ethanol_fe.conf"):
         head_tail = root.split(os.sep)
         for direc in head_tail:
            for key in model2BestWolfAlphaRCut.keys():
               if(key == direc):
                  for element in itertools.product(WolfMethods, WolfPotentials):
                     if(element[0]+"_"+element[1] in root):
                        path2File = os.path.join(root, name)
                        print(path2File)
                        print(key, element[0], element[1], box)
                        print(model2BestWolfAlphaRCut[key][(element[0], element[1], box)])
                        alphaRcutRelErrTuple = model2BestWolfAlphaRCut[key][(element[0], element[1], box)]
                        print("Cutoff" , alphaRcutRelErrTuple[0])
                        print("Alpha" , alphaRcutRelErrTuple[1])
"""
                        with open(path2File, "a") as myfile:
                            defPotLine = "Wolf\tTrue\t{pot}\n".format(pot=WolfDefaultPotential)
                            myfile.write(defPotLine)
                            defKindLine = "WolfKind\t{kind}\n".format(kind=WolfDefaultKind)
                            myfile.write(defKindLine)
"""


