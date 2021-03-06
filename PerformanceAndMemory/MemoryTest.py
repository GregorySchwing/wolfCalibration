# GOMC Example for the Gibbs Ensemble (GEMC) using MoSDeF [1, 2, 5-10, 13-17]

# Note: In this specific example, we will be using the GEMC_NVT ensemble.


# Import the required packages and specify the force field (FF) being used. 

# Note: For GOMC, the residue names are treated as molecules, so the residue names must be unique for each different molecule. [1, 2, 13-17]

# Note: Each residue can be set to a different FF, which is done by setting the residue name to a FF in a dictionary (FF_Dict).  The FF selection can be a FF name (set from foyer FF repositor) or a specified FF xml file. [1, 2, 13-17]
import subprocess
import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
import shutil
import pathlib
import random
from pathlib import Path
import itertools
import os
#Trappe SPCE
FF_file_water = 'xml/tip3p.xml'
water = mb.load('O', smiles=True)
water.name = 'H2O'
water.energy_minimize(forcefield=FF_file_water, steps=10**5)

# FF descriptors
FF_dict = {water.name: FF_file_water}
residues_list = [water.name]
fix_bonds_angles_residues = [water.name]

"""

Liquid phase systems contained one solute in a solvent
box of 200 1-octanol, 150 n-hexadecane, or 1000 water
molecules. Initial cubic box sizes were selected to produce
densities that were close to equilibrium, with a side length
of 37.6, 41.6, and 31.3 Å for 1-octanol, n-hexadecane,
and water, respectively.

"""


Molecule_Type_List = [water]
Molecule_Num_List = [1000]


# Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]
#charmm.write_inp()
#charmm.write_psf()
#charmm.write_pdb()


# Charmm writer needs files to actually exist
# So we write them at this location then move the folder
#this into the cwd

systems = [25, 35, 50]
for boxLength in range(1, 7):
    systems.append(boxLength*100)

system2Mem = dict()
for sys in systems:
    if(sys < 200):
        system2Mem[sys] = 1
    elif(sys <= 300):
        system2Mem[sys] = 3
    elif(sys == 400):
        system2Mem[sys] = 7
    elif(sys == 500):
        system2Mem[sys] = 15
    elif(sys == 600):
        system2Mem[sys] = 25
numSystems = len(systems)

calibPathSuffix = "_a"

methods = ["Wolf", "Ewald"]

runsdir = Path("memoryRuns")
sysdir = Path("systems")
runsdir.mkdir(parents=True, exist_ok=True)

box = "0"
WolfDefaultKind = "Vlugt"
WolfDefaultPotential = "DSP"
defaultAlpha = "0.21"
#for r in range(0, numReplicates, 1):

numReplicates = 1

random.seed(123)
randomSeeds = []

for r in range(0, numReplicates, 1):
    randomSeeds.append(random.randint(0,9999999))


replicatePathPrefix = "TI_"
replicatePaths = []
for r in range(0, numReplicates, 1):
    replicatePaths.append(Path(replicatePathPrefix + str(r)))
print(replicatePaths)


for boxLength, method in itertools.product(systems, methods):

    path2System ="../../../systems/" + str(boxLength) + "_a"

    liquid_box_length_Ang = boxLength
    # Build the main simulation liquid box (box 0) [1, 2, 13-17]
    water_ethanol_box_liq = mb.fill_box(compound=Molecule_Type_List,
                                        n_compounds=Molecule_Num_List,
                                        box=[liquid_box_length_Ang / 10,
                                             liquid_box_length_Ang / 10,
                                             liquid_box_length_Ang / 10])

    ## Build the Charmm object, which is required to write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

    ## The reorder_res_in_pdb_psf command reorders the psf and pdb to the order residues variable (i.e., the residues_list in this case) [1, 2, 13-17].  

    ## The fix_res_bonds_angles command fixes the angles and bonds for water in the Charmm FF file.  Note: This is specific to GOMC, as it sets the bond and angle k-values to 999999999999 [1, 2, 5-10, 13-17].
    charmm = mf_charmm.Charmm(water_ethanol_box_liq,
                              str(boxLength)+"_a",
                              ff_filename="NVT_water_FF",
                              forcefield_selection=FF_dict,
                              residues=residues_list,
                              bead_to_atom_name_dict=None,
                              fix_residue=None,
                              gomc_fix_bonds_angles=fix_bonds_angles_residues,
                              reorder_res_in_pdb_psf=True
                              )

    for r in range(0, numReplicates, 1):
        path = runsdir / replicatePaths[r] / Path(str(boxLength) + "_a_" + method)

        path.mkdir(parents=True, exist_ok=True)

        NVT_conf_name = "NVT_water.conf"
        OutputName = "NVT_water"

        NumSteps = 100
        Temp_in_K = 298
        Pressure_in_bar = 1.0

        input_variables_dict_val={"VDWGeometricSigma": True,
                           "DisFreq": 0.50,
                           "RotFreq": 0.20, 
                           "RegrowthFreq": 0.20,
                           "CrankShaftFreq": 0.10,
                           "CBMC_First" : 10,
                           "CBMC_Nth" : 10,
                           "CBMC_Ang" : 100,
                           "CBMC_Dih" : 50,
                           "Rcut": 12,
                           "RcutLow": 0,
                           "Ewald": method == "Ewald",
                           "LRC": True,
                           "RcutCoulomb_box_0": 12,
                           "Tolerance" : 0.00005,
                           "OutputName" : OutputName,
                           "PRNG" : randomSeeds[r],
                           "RestartFreq" : [False,1],
                           "CheckpointFreq" : [False,1],
                           "CoordinatesFreq" : [False,1],
                           "ConsoleFreq" : [True,1],
                           "BlockAverageFreq" : [False,1],
                           "HistogramFreq" : [False,1]
                           }

        orig_ff_file_path = Path("GOMC_toppar_water_ions_namd.str")
        renamed_ff_file_path = Path("NVT_water_FF.inp")
        top = Path("top")

        boxpath = sysdir / Path(str(boxLength) + "_a")

        # make optional
        boxpath.mkdir(parents=True, exist_ok=True)

        shutil.copy((top / orig_ff_file_path), (boxpath / renamed_ff_file_path))  # For Python 3.8+.

        gomc_control.write_gomc_control_file(charmm, NVT_conf_name, 'NVT', RunSteps=NumSteps, Temperature=Temp_in_K, Restart=False, check_input_files_exist=False, ff_psf_pdb_file_directory=path2System, input_variables_dict=input_variables_dict_val)

        if(method != "Ewald"):
             with open("NVT_water.conf", "a") as myfile:
                defPotLine = "Wolf\tTrue\t{pot}\n".format(pot=WolfDefaultPotential)
                myfile.write(defPotLine)
                defKindLine = "WolfKind\t{kind}\n".format(kind=WolfDefaultKind)
                myfile.write(defKindLine)
                defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=defaultAlpha)
                myfile.write(defAlphaLine)

        pathConf = Path(NVT_conf_name)
        pathConf.rename(path / pathConf)

        # Read in the file
        with open('memTemplatePotoffCluster.sh', 'r') as file :
          filedata = file.read()

        # Replace the target string
        filedata = filedata.replace('XXX', str(system2Mem[boxLength]))
        # Write the file out again
        with open('mem.sh', 'w') as file:
          file.write(filedata)

        bash = Path('mem.sh')
        bash.rename(path / bash)
        cwd = os.getcwd()
        os.chdir( path )
        bashCommand = "sbatch mem.sh"
        process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)
        output, error = process.communicate()
        os.chdir(cwd)


