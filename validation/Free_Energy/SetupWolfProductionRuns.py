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

from Surface import find_minimum

WolfMethods = ["Vlugt", "Gross", "VlugtWIntraCutoff"]
potentials = ["DSP", "DSF"]

p = re.compile("Wolf_Calibration_(\w+?)_(\w+?)_BOX_(\d+)_(\w+?).dat")

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
         head_tail = os.path.split(root)
         for direc in head_tail:
             for key in model2BestWolfAlphaRCut.keys():
                 if(key in direc):
                    print ("model" , key)
                    print ("wolf Kind" , wolfKind)
                    print ("potential Kind" , potential)
                    print ("box" , box)

                    tupleMin = find_minimum(os.path.join(root, name), True)
                    # Use smaller error, either BF or Grad Desc
                    if(tupleMin[2] < tupleMin[5]):
                       model2BestWolfAlphaRCut[key][(wolfKind, potential, box)] = [tupleMin[0], tupleMin[1]]
                    else:
                       model2BestWolfAlphaRCut[key][(wolfKind, potential, box)] = [tupleMin[3], tupleMin[4]]
         print(model2BestWolfAlphaRCut)
         quit()

"""
for root, dirs, files in os.walk(".", topdown=False):
for name in files:
  if(name == "NVT_Prod_water_ethanol_fe.conf"):
     print(os.path.join(root, name))
     print("DSF :" ,"DSF" in root)
     print("DSP :" ,"DSP" in root)
"""


