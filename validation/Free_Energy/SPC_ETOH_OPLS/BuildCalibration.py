# GOMC Example for the Gibbs Ensemble (GEMC) using MoSDeF [1, 2, 5-10, 13-17]

# Note: In this specific example, we will be using the GEMC_NVT ensemble.


# Import the required packages and specify the force field (FF) being used. 

# Note: For GOMC, the residue names are treated as molecules, so the residue names must be unique for each different molecule. [1, 2, 13-17]

# Note: Each residue can be set to a different FF, which is done by setting the residue name to a FF in a dictionary (FF_Dict).  The FF selection can be a FF name (set from foyer FF repositor) or a specified FF xml file. [1, 2, 13-17]

import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
import shutil
import pathlib
import random
from pathlib import Path
FF_file_water = '../../common/spc_trappe_opls_compatible.xml'
water = mb.load('O', smiles=True)
water.name = 'H2O'
water.energy_minimize(forcefield=FF_file_water, steps=10**5)

# OPLS ETOH
ethanol = mb.load('CCO', smiles=True)
ethanol.name = 'ETO'
FF_file_ethanol = 'oplsaa'
ethanol.energy_minimize(forcefield=FF_file_ethanol, steps=10**9)

# FF descriptors
FF_dict = {water.name: FF_file_water, ethanol.name: FF_file_ethanol}
residues_list = [ethanol.name, water.name]
fix_bonds_angles_residues = [water.name]
"""

Liquid phase systems contained one solute in a solvent
box of 200 1-octanol, 150 n-hexadecane, or 1000 water
molecules. Initial cubic box sizes were selected to produce
densities that were close to equilibrium, with a side length
of 37.6, 41.6, and 31.3 Å for 1-octanol, n-hexadecane,
and water, respectively.

"""


Molecule_Type_List = [water, ethanol]
Molecule_Num_List = [1000, 1]
liquid_box_length_Ang = 31.3
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
                          'NVT_water_ethanol_fe',
                          ff_filename="NVT_water_ethanol_fe_FF",
                          forcefield_selection=FF_dict,
                          residues=residues_list,
                          bead_to_atom_name_dict=None,
                          fix_residue=None,
                          gomc_fix_bonds_angles=fix_bonds_angles_residues,
                          reorder_res_in_pdb_psf=True
                          )


# Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]
#charmm.write_inp()
#charmm.write_psf()
#charmm.write_pdb()


# Charmm writer needs files to actually exist
# So we write them at this location then move the folder
#this into the cwd



"""
NVT ensemble simulations were performed with a move
ratio of 50% displacements, 20% rotations, 20% coupled-
decoupled configurational-bias (CD-CBMC) regrowth
[79], and 10% crankshaft [80,81]. Parameters for the
configurational-bias regrowth move were 100 angle tri-
als, 50 dihedral trials, and 10 trial locations for grown
pseudo-atoms. NPT ensemble simulations were per-
formed with similar move ratios, except for the addition
of 1% volume changes, while the percentage of displace-
ment moves was reduced to 49%. Non-bonded potentials
were truncated at 14 Å [48–50] and analytical tail cor-
rections were applied to the energy [82]. For simulations
with electrostatic interactions, the real space part of elec-
trostatic potential was truncated at 14 Å and an Ewald
convergence tolerance of 1 × 10 −5 was used [83].
The soft-core parameters used for Lennard-
Jones interactions were, α = 0.5, p = 2, and σ min = 3.0
[62,67].


To calculate the free energy of solvation/hydration, all
intermediate λ states were equilibrated independently
in the canonical ensemble (NVT) for 5 × 10 6 Monte
Carlo steps (MCS) at 298 K, followed by a 3 × 10 7 MCS
isobaric-isothermal (NPT) ensemble simulation at 1 bar
and 298 K. Production data were taken from a subse-
quent 5 × 10 7 MCS NPT simulation, which used the
final configuration of the prior NPT simulation as the
initial configuration.

During the
production run, the change in energy (DeltaU i,j ) between
the current lambda state and all other lambda states,
and the derivative of potential with respect to lambda
(dU Coul /dλ Coul , dU LJ /dλ LJ ), were evaluated and stored
for post-simulation analysis every 5 × 10 3 MCS.

To calculate the free energy of solvation in water and
1-octanol, 23 intermediate lambda states, as shown in
Figure 1, were used:
λ coul,LJ ∈ {
(0.0, 0.0), (0.0, 0.05), (0.0, 0.1), (0.0, 0.15),
(0.0, 0.2), (0.0, 0.25), (0.0, 0.3), (0.0, 0.35),
(0.0, 0.4), (0.0, 0.45), (0.0, 0.5), (0.0, 0.6),
(0.0, 0.7), (0.0, 0.8), (0.0, 0.9), (0.0, 1.0),
(0.2, 1.0), (0.4, 1.0), (0.6, 1.0), (0.7, 1.0),
(0.8, 1.0), (0.9, 1.0), (1.0, 1.0) }

followed by a 3 × 10 7 MCS
isobaric-isothermal (NPT) ensemble simulation at 1 bar
and 298 K.
"""

LambdaVDWList = []
LambdaCoulList = []
# Append 16 0.0's
for x in range(0, 16):
    LambdaVDWList.append(0.0)
# Append 0.2, 0.4, 0.6
for x in range(2, 8, 2):
    LambdaVDWList.append(round(x*0.1,2))
# Append 0.7, 0.8, 0.9, 1.0
for x in range(7, 11, 1):
    LambdaVDWList.append(round(x*0.1,2))

# 0.0-0.5, by 0.5
for x in range(0, 55, 5):
    LambdaCoulList.append(round(x*0.01,2))
# 0.6-0.9
for x in range(6, 10, 1):
    LambdaCoulList.append(round(x*0.1,2))
# Append 7 1.0's
for x in range(0, 8, 1):
    LambdaCoulList.append(1.0)

print(LambdaVDWList)
print(LambdaCoulList)

numReplicates = 5
random.seed(123)
randomSeeds = []

for r in range(0, numReplicates, 1):
    replicateSeeds = []
    for x in range(0, len(LambdaVDWList), 1):
        replicateSeeds.append(random.randint(0,9999999))
    randomSeeds.append(replicateSeeds)

print("Random Seeds:")
print(randomSeeds)

calibPathPrefix = "Calibration"

replicatePathPrefix = "TI_"
replicatePaths = []
for r in range(0, 1, 1):
    replicatePaths.append(Path(replicatePathPrefix + str(r)))
print(replicatePaths)

NVT_Eq_Prefix = Path("NVT_Eq")
NPT_Eq_Prefix = Path("NPT_Eq")

#for r in range(0, numReplicates, 1):
for r in range(0, 1, 1):

    NVT_Eq = replicatePaths[r] / NVT_Eq_Prefix
    NPT_Eq = replicatePaths[r] / NPT_Eq_Prefix
    NVT_Calibration = Path(calibPathPrefix)

    RelPathToNVTEq = Path("../../") / NVT_Eq
    RelPathToNPTEq = Path("../../") / NPT_Eq

    NVT_Calibration.mkdir(parents=True, exist_ok=True)

    prefix = "state_"

    NVT_Cal_conf_name = "NVT_Cal_water_ethanol_fe.conf"

    NVT_Eq_OutputName = "NVT_Eq"
    NPT_Eq_OutputName = "NPT_Eq"
    NVT_Cal_OutputName = "NVT_Cal"
    Restart_XSC_Suffix = "_BOX_0_restart.xsc"
    Restart_COOR_Suffix = "_BOX_0_restart.coor"

    NVT_Restart_XSC_path = Path(NVT_Eq_OutputName + Restart_XSC_Suffix)
    NPT_Restart_XSC_path = Path(NPT_Eq_OutputName + Restart_XSC_Suffix)

    NVT_Restart_COOR_path = Path(NVT_Eq_OutputName + Restart_COOR_Suffix)
    NPT_Restart_COOR_path = Path(NPT_Eq_OutputName + Restart_COOR_Suffix)

    NumCalSteps = 1000
    Temp_in_K = 298
    Pressure_in_bar = 1.0

    ff_psf_pdb_file_directory_name = "../../common"

    # Calibration of wolf paramter across all lambda windows, calibration tool will internally cycle through all wolf potentials and formulae
    #for x in range(0, len(LambdaVDWList)):
    for x in range(0, 1):


        stateName = prefix+str(x)
        statePath = Path(stateName)

        NVT_cal_state_path = NVT_Calibration / statePath
        NVT_cal_state_path.mkdir(parents=True, exist_ok=True)

        NPT_restart_files_state_path = RelPathToNPTEq / statePath 
        NPT_restart_coor = NPT_restart_files_state_path / NPT_Restart_COOR_path
        NPT_restart_xsc = NPT_restart_files_state_path / NPT_Restart_XSC_path

        input_variables_dict_NPT_Prod={"Pressure" : Pressure_in_bar,
                           "VDWGeometricSigma": True,
                           "DisFreq": 0.01,
                           "RotFreq": 0.01,
                           "MultiparticleFreq": 0.68,  
                           "RegrowthFreq": 0.20,
                           "CrankShaftFreq": 0.10,
                           "CBMC_First" : 10,
                           "CBMC_Nth" : 10,
                           "CBMC_Ang" : 100,
                           "CBMC_Dih" : 50,
                           "Rcut": 14,
                           "RcutLow": 0,
                           "LRC": True,
                           "Ewald": False,
                           "RcutCoulomb_box_0": 14,
                           "Tolerance" : 0.00005,
                           "LambdaVDW" : LambdaVDWList,
                           "LambdaCoulomb" : LambdaCoulList,
                           "FreeEnergyCalc" : [False, 5000],
                           "PressureCalc" : [True, 5000],
                           "MoleculeType" : [ethanol.name, 1],
                           "ScaleAlpha" : 0.5,
                           "ScalePower" : 2,
                           "MinSigma" : 3,
                           "InitialState" : x,
                           "OutputName" : NVT_Cal_OutputName,
                           "PRNG" : randomSeeds[r][x],
                           "Exclude" : "1-4"
                           }


        gomc_control.write_gomc_control_file(charmm, NVT_Cal_conf_name, 'NVT', RunSteps=NumCalSteps, Temperature=Temp_in_K, ff_psf_pdb_file_directory=ff_psf_pdb_file_directory_name, Restart=True, binCoordinates_box_0=str(NPT_restart_coor),extendedSystem_box_0=str(NPT_restart_xsc),check_input_files_exist=False,input_variables_dict=input_variables_dict_NPT_Prod
                                            )
        NVTCalConfPath = Path(NVT_Cal_conf_name)
        NVTCalConfPath.rename(NVT_cal_state_path / NVTCalConfPath)

