import mbuild as mb
import foyer
import mosdef_cassandra as mc
import unyt as u
import numpy as np
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
from mosdef_cassandra.utils.get_files import get_example_ff_path, get_example_mol2_path
import os
import shutil
# Load water with SPC/E geometry from mol2 file
#molecule = mb.load(get_example_mol2_path("spce"))

Water_res_name = 'H2O'
FF_file_water = 'files/spce.xml'
water = mb.load('O', smiles=True)
water.name = Water_res_name
water.energy_minimize(forcefield = FF_file_water , steps=10**9)
FF_Dict = {water.name: FF_file_water}

residues_List = [water.name]

Fix_bonds_angles_residues = [ water.name ]

# Load force field
spce = foyer.Forcefield(get_example_ff_path("spce"))

# Use foyer to apply force field
molecule_ff = spce.apply(water)


print('Running: filling liquid box')
water_box_liq = mb.fill_box(compound=[water],
                                    n_compounds=50,
                                    box=[3.0, 3.0, 3.0] )
print('Completed: filling liquid box')


print('Running: GOMC FF file, and the psf and pdb files')
charmm = mf_charmm.Charmm(water_box_liq,
                          'NVT_water',
                          structure_box_1 = None,
                          filename_box_1 = None,
                          ff_filename = 'NVT_water',
                          forcefield_selection = FF_Dict,
                          residues= residues_List ,
                          bead_to_atom_name_dict = None,
                          fix_residue = None,
                          gomc_fix_bonds_angles = Fix_bonds_angles_residues,
                          reorder_res_in_pdb_psf = True
                          )

os.chdir("Ewald_GOMC")

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()

os.chdir("../Ewald_Cassandra")

# Create box and species list
box_list = [water_box_liq]
species_list = [molecule_ff]
molecules_in_boxes = [[50]]
# Define the System
system = mc.System(box_list, species_list, mols_in_boxes=molecules_in_boxes)
# Define the MoveSet
moveset = mc.MoveSet("nvt", species_list)

system2 = mc.System(box_list, species_list, mols_in_boxes=molecules_in_boxes)
# Define the MoveSet
moveset = mc.MoveSet("nvt", species_list)

# Note here we need to use the angle_style="fixed" keyword argument
# SPC/E geometry is rigid; default angle style is "harmonic"
custom_args = {"angle_style": ["fixed"]}

# Run a simulation with at 300 K with 10000 MC moveset
mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=10000,
    temperature=300.0 * u.K,
    **custom_args
)

shutil.copyfile("box1.in.xyz", "../Ewald_GOMC/box1.in.xyz")
shutil.copyfile("box1.in.xyz", "../DSF_GOMC/box1.in.xyz")

os.chdir("../Ewald_GOMC")

os.system("vmd < xyz2bincoords.vmd")  
os.system("vmd < pdb2xsc.vmd")  

gomc_control.write_gomc_control_file(charmm, 'in_NVT.conf',  'NVT', 10000, 300, Restart=True, 
                                     binCoordinates_box_0="NVT_water.coor",
                                     extendedSystem_box_0="NVT_water.xsc",
                                     input_variables_dict={"VDWGeometricSigma": False,
                                                           "Potential": "VDW",
                                                           "LRC":  True,
                                                           "RCut": 12,
                                                           "RCut": 12,
                                                           "RcutCoulomb_box_0": 12,
                                                           "RcutLow": 1.0,
                                                           "Ewald": True,
                                                           "Tolerance": 1e-05,
                                                           "Exclude": "1-3"
                                                           }
                                     )

print('Completed: GOMC EWALD FF file, and the psf and pdb files')

os.chdir("../DSF_Cassandra")

mc.run(
    system=system2,
    moveset=moveset,
    run_type="equilibration",
    run_length=10000,
    temperature=300.0 * u.K,
    charge_style='dsf',
    charge_cutoff=12.0 * u.angstrom,
    dsf_damping=0.22
)

os.chdir("../DSF_GOMC")

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()

os.system("vmd < xyz2bincoords.vmd")  
os.system("vmd < pdb2xsc.vmd")  

gomc_control.write_gomc_control_file(charmm, 'in_NVT.conf',  'NVT', 10000, 300, Restart=True, 
                                     binCoordinates_box_0="NVT_water.coor",
                                     extendedSystem_box_0="NVT_water.xsc",
                                     input_variables_dict={"VDWGeometricSigma": False,
                                                           "Potential": "VDW",
                                                           "LRC":  True,
                                                           "RCut": 12,
                                                           "RCut": 12,
                                                           "RcutCoulomb_box_0": 12,
                                                           "RcutLow": 1.0,
                                                           "Ewald": True,
                                                           "Tolerance": 1e-05,
                                                           "Exclude": "1-3"
                                                           }
                                     )

print('Completed: GOMC DSF FF file, and the psf and pdb files')
path2File = "in_NVT.conf"
with open(path2File, "a") as myfile:
    defPotLine = "Wolf\tTrue\t{pot}\n".format(pot="DSF")
    myfile.write(defPotLine)
    defKindLine = "WolfKind\t{kind}\n".format(kind="VLUGTWINTRACUTOFF")
    myfile.write(defKindLine)
    defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=0, val=0.22)
    myfile.write(defAlphaLine)

