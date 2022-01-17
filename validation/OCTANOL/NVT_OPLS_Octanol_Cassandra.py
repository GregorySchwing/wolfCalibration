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
octanol_res_name = 'OCT'
FF_file_octanol = 'oplsaa'
octanol = mb.load('CCCCCCCCO', smiles=True)
octanol.name = octanol_res_name
octanol.energy_minimize(forcefield = FF_file_octanol)
FF_Dict = {octanol.name: FF_file_octanol}
residues_List = [octanol.name]
print('Running: filling liquid box')


box_liq = mb.fill_box(compound=[octanol],
                                    n_compounds=50,
                                    box=[3.0, 3.0, 3.0] )
print('Completed: filling liquid box')
# Load force field
oplsaa = foyer.forcefields.load_OPLSAA()

# Use foyer to apply force field
molecule_ff = oplsaa.apply(octanol)

# Create box and species list
box_list = [box_liq]
species_list = [molecule_ff]
molecules_in_boxes = [[50]]
# Use Cassandra to insert some initial number of species
#mols_to_add = [[50]]

# Define the System
system = mc.System(box_list, species_list, mols_in_boxes=molecules_in_boxes)
# Define the MoveSet
moveset = mc.MoveSet("nvt", species_list)

print('Running: GOMC FF file, and the psf and pdb files')
charmm = mf_charmm.Charmm(box_liq,
                          'NVT_octanol',
                          structure_box_1 =None  ,
                          filename_box_1 = None,
                          ff_filename ="NVT_octanol_FF" ,
                          forcefield_selection  = FF_Dict ,
                          residues= residues_List ,
                          gomc_fix_bonds_angles = None,
                         )

os.chdir("Ewald_GOMC")

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()

os.chdir("../Ewald_Cassandra")

# Run a simulation with at 300 K with 10000 MC moveset
mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=10000,
    temperature=300.0 * u.K#,
    #charge_style='dsf',
    #charge_cutoff=12.0 * u.angstrom,
    #dsf_damping=0.22
)


shutil.copyfile("box1.in.xyz", "../Ewald_GOMC/box1.in.xyz")
shutil.copyfile("box1.in.xyz", "../DSF_GOMC/box1.in.xyz")

os.chdir("../Ewald_GOMC")

os.system("vmd < xyz2bincoords.vmd")  
os.system("vmd < pdb2xsc.vmd")  


gomc_control.write_gomc_control_file(charmm, 'in_NVT.conf',  'NVT', 10000, 300, Restart=True, 
                                     binCoordinates_box_0="NVT_octanol.coor",
                                     extendedSystem_box_0="NVT_octanol.xsc",
                                     input_variables_dict={"VDWGeometricSigma": False,
                                                           "Potential": "VDW",
                                                           "LRC":  True,
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
    system=system,
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
                                     binCoordinates_box_0="NVT_octanol.coor",
                                     extendedSystem_box_0="NVT_octanol.xsc",
                                     input_variables_dict={"VDWGeometricSigma": False,
                                                           "Potential": "VDW",
                                                           "LRC":  True,
                                                           "RCut": 12,
                                                           "RcutCoulomb_box_0": 12,
                                                           "RcutLow": 1.0,
                                                           "Ewald": True,
                                                           "Tolerance": 1e-05,
                                                           "Exclude": "1-3"
                                                           }
                                     )

print('Completed: GOMC DSF FF file, and the psf and pdb files')

