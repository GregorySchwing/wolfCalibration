import mbuild as mb
from mbuild.lattice import load_cif
from mbuild.utils.io import get_fn
import mbuild as mb
import numpy as np
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control

# load water and build the water box reservoir
FF_file_water = '../common/spce.xml'
water = mb.load('O', smiles=True)
water.name = 'H2O'
water.energy_minimize(forcefield = FF_file_water , steps=10**4)

water_box_reservior = mb.fill_box(compound=[water],
                                  density= 700,
                                  compound_ratio = [1],
                                  box=[6, 6, 6])


# load and build the ETV with 3 cell units in every direction (x, y, and z)
lattice_cif_ETV_triclinic = load_cif("../common/ETV_triclinic.cif")
ETV_triclinic_3_cell = lattice_cif_ETV_triclinic.populate(x=3, y=3, z=3)
ETV_triclinic_3_cell.name = 'ETV'

# build the Charmm object and write the FF, PSF, and PDB files
# Note the FF for ETV is not the correct Physics and is made up for testing.
charmm = mf_charmm.Charmm(ETV_triclinic_3_cell,
                          'ETV_triclinic_3_cell_box_0',
                          structure_box_1=water_box_reservior,
                          filename_box_1="water_box_reservior",
                          ff_filename="ETV_triclinic_water_FF",
                          forcefield_selection={
                              ETV_triclinic_3_cell.name: "../common/Charmm_writer_testing_only_zeolite.xml",
                              water.name: "../common/spce.xml"},
                          residues=[ETV_triclinic_3_cell.name,
                                    water.name],
                          bead_to_atom_name_dict=None,
                          fix_residue=[ETV_triclinic_3_cell.name],
                          )

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()


# Write the GOMC control file
gomc_control.write_gomc_control_file(charmm,
                                     'ETV_triclinic_water.conf',
                                     'GCMC',
                                     100000,
                                     300,
                                     input_variables_dict = {'ChemPot': {ETV_triclinic_3_cell.name: 0,
                                                                         water.name: -4000}
                                                             }
                                     )