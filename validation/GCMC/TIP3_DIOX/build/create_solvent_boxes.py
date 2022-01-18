from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import numpy as np
import shutil

f = open("log.txt", "a")

parser = PDBParser(PERMISSIVE=1)
structure_id = "3rgk"
filename = "../1-1-build/MYO_HEME_MUT.pdb"
structure = parser.get_structure(structure_id, filename)
atoms = structure.get_atoms()
listOfCoords = []
for atom in atoms:
	coords = atom.get_coord()
	listOfCoords.append(coords)
coorNP = np.asarray(listOfCoords)
geoCenter = coorNP.mean(axis=0)
log = "Geometric Center: {}\n".format(geoCenter)
f.write(log)

# calculating Euclidean distance
# using linalg.norm()
maxDistMan = 0
maxDistL2 = 0
for atom in coorNP:
	#manDist = np.sqrt((atom[0] - geoCenter[0])**2 + (atom[1] - geoCenter[1])**2 + (atom[2] - geoCenter[2])**2)
	dist = np.linalg.norm(geoCenter - atom)
	#print ("manDist {} dist {}".format(manDist,dist))

	#print(dist)	
	if (dist > maxDistL2):
		maxDistL2 = dist
	#if (manDist > maxDistMan):
	#	maxDistMan = manDist

log = "maxDistL2 {}\n".format(maxDistL2)
f.write(log)

# 2 times the Max internal distance of protein atoms + 2 angstroms on each side
# The padding is in case the radius of gyration of the protein increases.
# Currently the maximum allowed increase in radius of gyration is 2 angstroms.
# This is likely a highly liberal amount for a globular protein at 310 K in minimal Na/Cl.
maxDistL2_padded = maxDistL2+20
log = "maxDistL2_padded {}\n".format(maxDistL2_padded)
f.write(log)

shutil.copyfile("../1-1-build/MYO_HEME.psf", "MYO_HEME_SHIFTED.psf")


import mbuild as mb
import numpy as np
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control

FF_file_O2 = './FFs/charmmD_molecular_oxygen.xml'
O2 = mb.load('./FFs/DIOX.mol2')
O2.name = 'DIOX'
#O2.energy_minimize(forcefield=FF_file_O2, steps=10**5)
FF_file_water = './FFs/charmm_tip3p.xml'
water = mb.load('O', smiles=True)
water.name = 'TIP3'
#water.energy_minimize(forcefield=FF_file_water, steps=10**5)

FF_dict = {water.name: FF_file_water, O2.name: FF_file_O2}
residues_list = [water.name, O2.name]
fix_bonds_angles_residues = [water.name, O2.name]
bead_to_atom_name_dict = { '_ON':'ON', '_OP':'OP'}

# Build the main simulation liquid box (box 0) and the vapor (box 1) for the simulation [1, 2, 13-17]


water_O2_box_liq = mb.fill_box(compound=[water,O2],
                                    density= 950,
                                    compound_ratio=[0.98, 0.02] ,
                                    box=[2*maxDistL2_padded/10, 2*maxDistL2_padded/10, 2*maxDistL2_padded/10])


geoCenterBox = water_O2_box_liq.center
log = "BOX CENTER :  {}\n".format(geoCenterBox*10)
f.write(log)

trueCenter = [maxDistL2_padded/10, maxDistL2_padded/10, maxDistL2_padded/10]
log = "DESIRED BOX CENTER : {}\n".format(trueCenter*10)
f.write(log)

translationVectorBox = trueCenter-geoCenterBox
log = "BOX TRANSLATION VECTOR : {}\n".format(translationVectorBox*10)
f.write(log)

water_O2_box_liq.translate(translationVectorBox)

geoCenterBoxPostTranslate = water_O2_box_liq.center
log = "BOX CENTER POST TRANSLATE : {}\n".format(geoCenterBoxPostTranslate*10)
f.write(log)

water_O2_box_res = mb.fill_box(compound=[water,O2],
                                    density= 950,
                                    compound_ratio=[0.80 0.20] ,
                                    box=[9, 9, 9])


charmmNAMD = mf_charmm.Charmm(water_O2_box_liq,
                          'GCMC_water_O2_liq_NAMD',
                          structure_box_1=water_O2_box_res,
                          filename_box_1='GCMC_water_O2_res_NAMD',
                          ff_filename="GCMC_water_O2_FF_NAMD",
                          forcefield_selection=FF_dict,
                          residues=residues_list,
                          bead_to_atom_name_dict=bead_to_atom_name_dict,
                          fix_residue=None,
                          gomc_fix_bonds_angles=None,
                          reorder_res_in_pdb_psf=True
                          )


charmm = mf_charmm.Charmm(water_O2_box_liq,
                          'GCMC_water_O2_liq',
                          structure_box_1=water_O2_box_res,
                          filename_box_1='GCMC_water_O2_res',
                          ff_filename="GCMC_water_O2_FF",
                          forcefield_selection=FF_dict,
                          residues=residues_list,
                          bead_to_atom_name_dict=bead_to_atom_name_dict,
                          fix_residue=None,
                          gomc_fix_bonds_angles=fix_bonds_angles_residues,
                          reorder_res_in_pdb_psf=True
                          )


charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()

charmmNAMD.write_inp()

gomc_control.write_gomc_control_file(charmm, 'in_GCMC_NVT.conf', 'GCMC', 100, 310,
                                     input_variables_dict={"VDWGeometricSigma": True,
                                                           "Rcut": 12,
                                                           "DisFreq": 0.00,
                                                           "RotFreq": 0.00, 
                                                           "IntraSwapFreq": 0.00,
                                                           "SwapFreq": 1.00,
                                                           "RegrowthFreq": 0.00,
                                                           "CrankShaftFreq": 0.00,
                                                           "VolFreq": 0.00,
                                                           "MultiParticleFreq": 0.00,
                                                           "ChemPot" : {"TIP3" : -4166, "DIOX" : -8000}
                                                           }
                                    )

f.write('Completed: GOMC FF file, and the psf and pdb files')

log = "PROTEIN GEOMETRIC CENTER: {}\n".format(geoCenter)
f.write(log)

log = "BOX GEOMETRIC CENTER: {}\n".format(geoCenterBoxPostTranslate*10)
f.write(log)

translationArrayProt = np.abs(geoCenterBoxPostTranslate*10 - geoCenter)
log = "PROTEIN TRANSLATION VECTOR : {}\n".format(translationArrayProt)
f.write(log)

atoms = structure.get_atoms()
for atom in atoms:
	newCoords = atom.get_coord()+translationArrayProt
	atom.set_coord(newCoords)
io = PDBIO()
io.set_structure(structure)
io.save("MYO_HEME_MUT_SHIFTED.pdb")
