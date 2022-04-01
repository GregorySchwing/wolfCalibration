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

WolfMethods = ["Vlugt", "Gross", "VlugtWIntraCutoff"]

for root, dirs, files in os.walk(".", topdown=False):
   for directory in dirs:
      if("DSP" in directory):
         print("DSP Directory" ,os.path.join(root, directory))
      if("DSF" in directory):
         print("DSF Directory" ,os.path.join(root, directory))
"""
   for name in files:
      if(name == "NVT_Prod_water_ethanol_fe.conf"):
         print(os.path.join(root, name))
         print("DSF :" ,"DSF" in root)
         print("DSP :" ,"DSP" in root)
"""


