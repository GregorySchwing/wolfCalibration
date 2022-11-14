"""GOMC's setup for signac, signac-flow, signac-dashboard for this study."""
# project.py

import math
import flow
# from flow.environment import StandardEnvironment
import mbuild as mb
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
import numpy as np

from alchemlyb.parsing.gomc import  extract_dHdl,  extract_u_nk
from alchemlyb.estimators import MBAR, BAR, TI
import alchemlyb.preprocessing.subsampling as ss
import pandas as pd
import numpy as np
import os

import unyt as u
from flow import FlowProject, aggregator
from flow.environment import DefaultSlurmEnvironment

from src.utils.forcefields import get_ff_path
from src.utils.forcefields import get_molecule_path
from templates.NAMD_conf_template import generate_namd_equilb_control_file
from templates.sphere_builder_template import get_sphere_builder_path, get_pdb2bincoords_path, get_pdb2xsc_path, get_water_box_builder_path


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()

class Potoff(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""

    #hostname_pattern = r".*reslab32ai8111"
    template = "potoff.sh"

#class Grid(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
#    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""#

#    hostname_pattern = r".*\.grid\.wayne\.edu"
#    template = "grid.sh"

# ******************************************************
# users typical variables, but not all (start)
# ******************************************************
# set binary path to gomc binary files (the bin folder).
# If the gomc binary files are callable directly from the terminal without a path,
# please just enter and empty string (i.e., "" or '')

# WSU grid binary paths
#gomc_binary_path = "/wsu/home/go/go24/go2432/wolf/GOMC/bin"
#namd_binary_path = "/wsu/home/go/go24/go2432/NAMD_2.14_Linux-x86_64-multicore-CUDA"

gomc_binary_path = "/wsu/home/go/go24/go2432/wolfCalibrationLong/validation/GEMC/signac/bin"
namd_binary_path = "/wsu/home/go/go24/go2432/wolfCalibrationLong/validation/GEMC/signac/bin"
# Potoff cluster bin paths
# Potoff cluster bin paths
gomc_binary_path = "/home6/go2432/wolfCalibration/validation/GEMC/signac/bin"
namd_binary_path = "/home6/go2432/wolfCalibration/validation/GEMC/signac/bin"

# local bin paths
#gomc_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/DensityExperiment/signac/bin"
#namd_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/DensityExperiment/signac/bin/NAMD_Git-2022-07-21_Linux-x86_64-multicore-CUDA"

#gomc_binary_path = "/home/greg/Desktop/wolfCalibration/validation/GEMC/signac/bin"
#namd_binary_path = "/home/greg/Desktop/wolfCalibration/validation/GEMC/signac/bin"

#WSL local bin paths
#gomc_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"
#namd_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"

# brads workstation binary paths
#gomc_binary_path = "/home/brad/Programs/GOMC/GOMC_dev_1_21_22/bin"
#namd_binary_path = "/home/brad/Programs/NAMD/NAMD_2.14_RTX_3080_build_Source_CUDA"

# number of simulation steps
gomc_steps_equilb_design_ensemble = 5 * 10**6 #  set value for paper = 60 * 10**6
gomc_steps_production = 100 * 10**6 # set value for paper = 60 * 10**6
gomc_console_output_data_every_X_steps = 50 * 10**3# set value for paper = 100 * 10**3

gomc_output_data_every_X_steps = 50 * 10**3 # # set value for paper = 50 * 10**3

"""
During the
production run, the change in energy (DeltaU i,j ) between
the current lambda state and all other lambda states,
and the derivative of potential with respect to lambda
(dU Coul /dλ Coul , dU LJ /dλ LJ ), were evaluated and stored
for post-simulation analysis every 5 × 10 3 MCS.
"""
gomc_free_energy_output_data_every_X_steps = 5 * 10**3 # set value for paper = 10 * 10**3

# calc MC steps
MC_steps = int(gomc_steps_equilb_design_ensemble)
EqSteps = 1000
Calibration_MC_steps = 500000
Calibration_MC_Eq_Steps = 10000
Wolf_Sanity_MC_steps = 50 * 10**6 # set value for paper = 60 * 10**6


# Free energy calcs: set free energy data in doc
# this number will generate the lamdas
# set the number of lambda spacings, which includes 0 to 1
#number_of_lambda_spacing_including_zero_int = 11
number_of_lambda_spacing_including_zero_int = 23


# force field (FF) file for all simulations in that job
# Note: do not add extensions
namd_ff_filename_str = "in_namd_FF"
gomc_ff_filename_str = "in_gomc_FF"

# initial mosdef structure and coordinates
# Note: do not add extensions
mosdef_structure_box_0_name_str = "mosdef_box_0"
mosdef_structure_box_1_name_str = "mosdef_box_1"

# melt equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
namd_equilb_NVT_control_file_box_0_name_str = "namd_equilb_NVT_box_0"
namd_equilb_NVT_control_file_box_1_name_str = "namd_equilb_NVT_box_1"

# The equilb using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
gomc_equilb_design_ensemble_control_file_name_str = "gomc_equilb_design_ensemble"

# The production run using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
#gomc_production_control_file_name_str = "gomc_production_run"
gomc_production_control_file_name_str = "wolf_sanity"

# Analysis (each replicates averages):
# Output text (txt) file names for each replicates averages
# directly put in each replicate folder (.txt, .dat, etc)
output_replicate_txt_file_name_liq = "analysis_avg_data_box_liq.txt"
output_replicate_txt_file_name_vap = "analysis_avg_data_box_vap.txt"

# Analysis (averages and std. devs. of  # all the replcates): 
# Output text (txt) file names for the averages and std. devs. of all the replcates, 
# including the extention (.txt, .dat, etc)
output_avg_std_of_replicates_txt_file_name_liq = "analysis_avg_std_of_replicates_box_liq.txt"
output_avg_std_of_replicates_txt_file_name_vap = "analysis_avg_std_of_replicates_box_vap.txt"

# Analysis (Critical and boiling point values):
# Output text (txt) file name for the Critical and Boiling point values of the system using the averages
# and std. devs. all the replcates, including the extention (.txt, .dat, etc)
output_critical_data_replicate_txt_file_name = "analysis_critical_points_all_replicates.txt"
output_critical_data_avg_std_of_replicates_txt_file_name = "analysis_critical_points_avg_std_of_replicates.txt"

output_boiling_data_replicate_txt_file_name = "analysis__boiling_point_all_replicates.txt"
output_boiling_data_avg_std_of_replicates_txt_file_name = "analysis_boiling_point_avg_std_of_replicates.txt"


walltime_mosdef_hr = 24
walltime_namd_hr = 24
#CPU
walltime_gomc_equilbrium_hr = 120
#GPU
#walltime_gomc_equilbrium_hr = 48
#CPU
walltime_gomc_production_hr = 240
#GPU
#walltime_gomc_production_hr = 96
walltime_gomc_analysis_hr = 4
memory_needed = 6



# forcefield names dict
# The only difference in FF for waters b/w opls and trappe is coulombic and LJ scaling
forcefield_residue_to_ff_filename_dict = {
    "MSPCE": {"OPLS" : "mspce_opls.xml", "TRAPPE" : "mspce_trappe.xml"},
    "SPCE": {"OPLS" : "spce_opls.xml", "TRAPPE" : "spce_trappe.xml"},
    "SPC": {"OPLS" : "spc_opls.xml", "TRAPPE" : "spc_trappe.xml"},
    "TIP3": {"OPLS" : "tip3_opls.xml", "TRAPPE" : "tip3_trappe.xml"},
    "TIP4": {"OPLS" : "tip4p_opls.xml", "TRAPPE" : "tip4p_trappe.xml"},
    "ETOH": {"OPLS" : "oplsaa", "TRAPPE" : "trappe-ua.xml"},
    "Ne": "nobel_gas_vrabec_LB_mixing.xml",
    "Rn": "nobel_gas_vrabec_LB_mixing.xml",
}

# smiles of mol2 file input a .mol2 file or smiles as a string
smiles_or_mol2_name_to_value_dict = {
    "MSPCE": {"OPLS" : "mspce.mol2", "TRAPPE" : "mspce.mol2"},
    "MSPCE": {"OPLS" : "O", "TRAPPE" : "O"},
    "SPCE": {"OPLS" : "O", "TRAPPE" : "O"},
    "SPC": {"OPLS" : "O", "TRAPPE" : "O"},
    "TIP3": {"OPLS" : "O", "TRAPPE" : "O"},
    "TIP4": {"OPLS" : 'tip4p.mol2', "TRAPPE" : 'tip4p.mol2'},
    "Ne": {"OPLS" : 'Ne', "TRAPPE" : 'Ne'},
    "Rn": {"OPLS" : 'Rn', "TRAPPE" : 'Rn'},
    "ETOH": {"OPLS" : "CCO", "TRAPPE" : "ethanol.mol2"},
}




# ******************************************************
# users typical variables, but not all (end)
# ******************************************************


# ******************************************************
# signac and GOMC-MOSDEF code (start)
# ******************************************************

# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (start).
# ******************************************************
# ******************************************************

# set the default directory
project_directory_path = str(os.getcwd())
print("project_directory_path = " +str(project_directory_path))


# ******************************************************
# ******************************************************
# functions for selecting/grouping/aggregating in different ways (start)
# ******************************************************
# ******************************************************

def statepoint_without_replica(job):
    keys = sorted(tuple(i for i in job.sp.keys() if i not in {"replica_number_int"}))
    return [(key, job.sp[key]) for key in keys]

def statepoint_without_temperature(job):
    keys = sorted(tuple(i for i in job.sp.keys() if i not in {"production_temperature_K"}))
    return [(key, job.sp[key]) for key in keys]

def append_wolf_calibration_parameters(job):
    import math
    WolfDefaultKind = "VlugtWIntraCutoff"
    WolfDefaultPotential = "DSP"
    WolfDefaultAlpha = [0.21]

    if job.doc.production_ensemble in ["GEMC_NVT", "GEMC_NPT"]:
        WolfCutoffBoxList = [0,1]

        WolfCutoffCoulombLowerBoundList = [10,10]
        WolfCutoffCoulombUpperBoundList = [math.floor(job.doc.liq_box_lengths_ang/2.0),math.floor(job.doc.vap_box_lengths_ang/2.0)]
        WolfCutoffCoulombIntervalList = [(math.floor(job.doc.liq_box_lengths_ang/2.0) - 10)/50,(math.floor(job.doc.vap_box_lengths_ang/2.0) - 10)/50]

        WolfAlphaLowerBoundList = [0.0, 0.0]
        WolfAlphabUpperBoundList = [0.5, 0.5]
        WolfAlphaIntervalList = [0.01, 0.01]
    else:
        WolfCutoffBoxList = [0]

        WolfCutoffCoulombLowerBoundList = [10]
        WolfCutoffCoulombUpperBoundList = [math.floor(job.doc.liq_box_lengths_ang/2.0)]
        WolfCutoffCoulombIntervalList = [(math.floor(job.doc.liq_box_lengths_ang/2.0) - 10)/50]

        WolfAlphaLowerBoundList = [0.0]
        WolfAlphabUpperBoundList = [0.5]
        WolfAlphaIntervalList = [0.01]

    wolfCalFreq = gomc_output_data_every_X_steps

    with open(job.fn("wolf_calibration.conf"), "a") as myfile:
        defPotLine = "Wolf\tFalse\n"
        myfile.write(defPotLine)
        defPotLine = "WolfPotential\t{pot}\n".format(pot=WolfDefaultPotential)
        myfile.write(defPotLine)
        defKindLine = "WolfKind\t{kind}\n".format(kind=WolfDefaultKind)
        myfile.write(defKindLine)
        defPotLine = "WolfCalibrationFreq\tTrue\t{freq}\n".format(freq=wolfCalFreq)
        myfile.write(defPotLine)
        for box, wolfCutoffLower, wolfCutoffUpper, wolfCutoffInterval, wolfAlphaLower, wolfAlphaUpper, wolfAlphaInterval, \
        in zip(WolfCutoffBoxList, WolfCutoffCoulombLowerBoundList, WolfCutoffCoulombUpperBoundList, WolfCutoffCoulombIntervalList, \
        WolfAlphaLowerBoundList, WolfAlphabUpperBoundList, WolfAlphaIntervalList):
            CutoffLine = "WolfCutoffCoulombRange\t{box}\t{lb}\t{ub}\t{inter}\n".format(box=box, lb=wolfCutoffLower, ub=wolfCutoffUpper, inter=wolfCutoffInterval)
            print(CutoffLine)
            myfile.write(CutoffLine)

            alphaLine = "WolfAlphaRange\t{box}\t{lb}\t{ub}\t{inter}\n".format(box=box, lb=wolfAlphaLower, ub=wolfAlphaUpper, inter=wolfAlphaInterval)
            print(alphaLine)
            myfile.write(alphaLine)

def append_checkpoint_line(job, config_file_name, path_to_previous_checkpoint_file):
    with open(job.fn("{}.conf".format(config_file_name)), "a") as myfile:
        checkpointLine = "Checkpoint\tTrue\t{}\n".format(path_to_previous_checkpoint_file)
        myfile.write(checkpointLine)
            
# ******************************************************
# ******************************************************
# functions for selecting/grouping/aggregating in different ways (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# functions for free energy calcs MBAR, TI, and BAR for getting delta free energy and delta error (start)
# ******************************************************
# ******************************************************

def get_delta_TI_or_MBAR(TI_or_MBAR_estimate, k_b_T):
    """ Return the change in free energy and standard deviation for the MBAR and TI estimates.

    """
    delta = TI_or_MBAR_estimate.delta_f_.iloc[0, -1] * k_b_T
    std_delta = TI_or_MBAR_estimate.d_delta_f_.iloc[0, -1] * k_b_T
    return delta, std_delta


def get_delta_BAR(BAR_estimate, k_b_T):
    """ Return the change in free energy and standard deviation for the BAR estimates.

    """
    error_estimate = 0.0

    for i in range(len(BAR_estimate.d_delta_f_) - 1):
        error_estimate += BAR_estimate.d_delta_f_.values[i][i + 1] ** 2

    delta = BAR_estimate.delta_f_.iloc[0, -1] * k_b_T
    std_delta = k_b_T * error_estimate ** 0.5
    return delta, std_delta

# ******************************************************
# ******************************************************
# functions for free energy calcs MBAR, TI, and BAR for getting delta free energy and delta error (end)
# ******************************************************
# ******************************************************

@Project.label
def part_1a_initial_data_input_to_json(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    if job.isfile(f"{'signac_job_document.json'}"):
        data_written_bool = True

    return data_written_bool


@Project.post(part_1a_initial_data_input_to_json)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def initial_parameters(job):
    """Set the initial job parameters into the jobs doc json file."""
    # select


    LambdaVDW_list = []
    LambdaCoul_list = []
    InitialState_list = []
    if job.sp.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
        for lamda_i in range(0, int(number_of_lambda_spacing_including_zero_int)):
            lambda_space_increments = 1 / int(number_of_lambda_spacing_including_zero_int - 1)
            LambdaVDW_list.append(np.round(lamda_i * lambda_space_increments, decimals=8))
            InitialState_list.append(lamda_i)

        """
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
        """
    elif job.sp.solute in ["ETOH"]:
        counter = 0
        # Append 16 0.0's
        for x in range(0, 16):
            LambdaCoul_list.append(0.0)
            InitialState_list.append(counter)
            counter = counter + 1
        # Append 0.2, 0.4, 0.6
        for x in range(2, 8, 2):
            LambdaCoul_list.append(round(x*0.1,2))
            InitialState_list.append(counter)
            counter = counter + 1
        # Append 0.7, 0.8, 0.9, 1.0
        for x in range(7, 11, 1):
            LambdaCoul_list.append(round(x*0.1,2))
            InitialState_list.append(counter)
            counter = counter + 1

        # 0.0-0.5, by 0.5
        for x in range(0, 55, 5):
            LambdaVDW_list.append(round(x*0.01,2))
            #InitialState_list.append(counter)
            #counter = counter + 1
        # 0.6-0.9
        for x in range(6, 10, 1):
            LambdaVDW_list.append(round(x*0.1,2))
            #InitialState_list.append(counter)
            #counter = counter + 1
        # Append 7 1.0's
        for x in range(0, 8, 1):
            LambdaVDW_list.append(1.0)    
            #InitialState_list.append(counter)
            #counter = counter + 1
    elif (job.sp.solute in ["solvent_box"]):
        LambdaVDW_list = [0]
        LambdaCoul_list = [0]
        InitialState_list = [0]
    else:
        print("Didnt recognize solute", job.sp.solute)
    print("*********************")
    print("*********************")
    print("LambdaVDW_list = " + str(LambdaVDW_list))
    print("LambdaCoul_list = " + str(LambdaCoul_list))
    print("InitialState_list = " + str(InitialState_list))
    print("*********************")
    print("*********************")
    if LambdaVDW_list[0] != 0 and LambdaVDW_list[-1] != 1 :
        raise ValueError("ERROR: The selected lambda list values do not start with a 0 and end 1.")

    job.doc.LambdaVDW_list = LambdaVDW_list
    job.doc.LambdaCoul_list = LambdaCoul_list
    job.doc.InitialState_list = InitialState_list

    equilibration_ensemble = "GEMC_NVT"
    production_ensemble = "GEMC_NVT"

    # set the GOMC production ensemble temp, pressure, molecule, box dimenstion and residue names
    job.doc.equilibration_ensemble = equilibration_ensemble
    job.doc.production_ensemble = production_ensemble
    job.doc.production_pressure_bar = (1 * u.atm).to('bar')
    job.doc.production_temperature_K = job.sp.production_temperature_K
    g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
    kg_per_m3 = u.kg / (u.m * u.m * u.m)
    
    
    job.doc.liquid_density = (job.sp.liquid_density * g_per_cm3).to(kg_per_m3)
    job.doc.vapor_density = (job.sp.vapor_density * g_per_cm3).to(kg_per_m3)
    job.doc.N_liquid_solvent = job.sp.N_liquid_solvent
    job.doc.N_vapor_solvent = job.sp.N_vapor_solvent

    job.doc.solvent = job.sp.solvent
    """

    Liquid phase systems contained one solute in a solvent
    box of 200 1-octanol, 150 n-hexadecane, or 1000 water
    molecules. Initial cubic box sizes were selected to produce
    densities that were close to equilibrium, with a side length
    of 37.6, 41.6, and 31.3 Å for 1-octanol, n-hexadecane,
    and water, respectively.

    """
    """
    job.doc.N_liquid_solvent = 1000
    if (job.sp.solute == "solvent_box"):
        job.doc.N_liquid_solute = 0
    else:
        job.doc.N_liquid_solute = 1
    """

    job.doc.liq_box_lengths_ang = job.sp.liq_box_lengths_ang
    job.doc.vap_box_lengths_ang = job.sp.vap_box_lengths_ang

    if job.sp.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
        job.doc.Rcut_ang = 13 * u.angstrom  # this is the Rcut for GOMC it is the Rswitch for NAMD
    else:
        job.doc.Rcut_ang = 13 * u.angstrom  # this is the Rcut for GOMC it is the Rswitch for NAMD

    job.doc.Rcut_for_switch_namd_ang = 17 * u.angstrom  # Switch Rcut for NAMD's Switch function
    job.doc.neighbor_list_dist_namd_ang = 22 * u.angstrom # NAMD's neighbor list

    # list replica seed numbers
    replica_no_to_seed_dict = {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
    }

    job.doc.replica_number_int = replica_no_to_seed_dict.get(
        int(job.sp.replica_number_int)
    )

    # set solvent and solute in doc
    job.doc.solvent = job.sp.solvent
    job.doc.solute = job.sp.solute

    job.doc.namd_node_ncpu = 3
    job.doc.namd_node_ngpu = 1
    #job.doc.namd_node_ngpu = 0

    job.doc.gomc_ncpu = 4  # 1 is optimal but I want data quick.  run time is set for 1 cpu
    job.doc.gomc_ngpu = 1
    #job.doc.gomc_ngpu = 0

    # set rcut, ewalds
    job.doc.winningWolfPotential = ""
    job.doc.winningWolfModel = ""

    # get the namd binary paths
    if job.doc.namd_node_ngpu == 0:
        job.doc.namd_cpu_or_gpu = "CPU"

    elif job.doc.namd_node_ngpu == 1:
        job.doc.namd_cpu_or_gpu = "GPU"

    else:
        raise ValueError(
            "Tee NAMD CPU and GPU can not be determined as force field (FF) is not available in the selection, "
            "or GPU selection is is not 0 or 1."
        )

    # get the gomc binary paths
    if job.doc.gomc_ngpu == 0:
        job.doc.gomc_cpu_or_gpu = "CPU"

    elif job.doc.gomc_ngpu == 1:
        job.doc.gomc_cpu_or_gpu = "GPU"

    else:
        raise ValueError(
            "The GOMC CPU and GPU can not be determined as force field (FF) is not available in the selection, "
            "or GPU selection is is not 0 or 1."
        )

    job.doc.namd_equilb_NVT_gomc_binary_file = f"namd2"
    # set the initial iteration number of the simulation
    if equilibration_ensemble == "NPT":
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NPT"
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_NPT"
    elif equilibration_ensemble == "NVT":
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NVT"
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_NVT"
    elif job.doc.production_ensemble in ["GEMC_NVT", "GEMC_NPT"]:
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_GEMC"
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_GEMC"
        job.doc.production_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_GEMC"
    elif job.doc.production_ensemble in ["GCMC"]:
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_GCMC"
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_GCMC"
        job.doc.production_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_GCMC"
    else:
        raise ValueError(
            "ERROR: A wrong ensemble has been specified for the gomc binary file"
        )


# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (end).
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and force field (FF) files were written (start)
# ******************************************************
# ******************************************************

# check if GOMC-MOSDEF wrote the gomc files
# @Project.pre(select_production_ensemble)
@Project.label
@flow.with_job
def mosdef_input_written(job):
    """Check that the mosdef files (psf, pdb, and force field (FF) files) are written ."""
    file_written_bool = False

    if (
        job.isfile(f"{namd_ff_filename_str}.inp")
        and job.isfile(f"{gomc_ff_filename_str}.inp")
        and job.isfile(
            f"{mosdef_structure_box_0_name_str}.psf"
        )
        and job.isfile(
            f"{mosdef_structure_box_0_name_str}.pdb"
        )
    ):
        file_written_bool = True

    return file_written_bool


# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and FF files were written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC control file was written (start)
# ******************************************************
# ******************************************************
# function for checking if the GOMC control file is written
def gomc_control_file_written(job, control_filename_str):
    """General check that the gomc control files are written."""
    file_written_bool = False
    control_file = f"{control_filename_str}.conf"

    if job.isfile(control_file):
        with open(job.fn(f"{control_file}"), "r") as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "OutputName" in line:
                    split_move_line = line.split()
                    if split_move_line[0] == "OutputName":
                        file_written_bool = True

    return file_written_bool

# function for checking if the NAMD control file is written
def namd_control_file_written(job, control_filename_str):
    """General check that the NAMD control files are written."""
    file_written_bool = False
    control_file = f"{control_filename_str}.conf"
    if job.isfile(control_file):
        with open(job.fn(f"{control_file}"), "r") as fp:
            out_namd = fp.readlines()
            for i, line in enumerate(out_namd):
                if "cellBasisVector1" in line:
                    split_move_line = line.split()
                    if split_move_line[0] == "cellBasisVector1":
                        file_written_bool = True

    return file_written_bool


@Project.label
@flow.with_job
def part_2a_wolf_calibration_control_file_written(job):
    """General check that the namd_equilb_NVT_control_file
    (high temperature to set temp NAMD control file) is written."""
    if (job.sp.wolf_model != "Calibrator"):
        return True
    output_name_control_file_name = "wolf_calibration"
    try:
        return gomc_control_file_written(
            job,
            output_name_control_file_name,
        )
    except:
        return False


@Project.label
@flow.with_job
def part_2a_wolf_sanity_control_file_written(job):
    """General check that the namd_equilb_NVT_control_file
    (high temperature to set temp NAMD control file) is written."""
    if (job.sp.wolf_model == "Calibrator" or job.sp.electrostatic_method == "Ewald"):
        return True
    output_name_control_file_name = "wolf_sanity"
    try:
        return gomc_control_file_written(
            job,
            output_name_control_file_name,
        )
    except:
        return False

# checking if the NAMD control file is written for the melt equilb NVT run
@Project.label
@flow.with_job
def part_2a_namd_equilb_NVT_box_0_control_file_written(job):
    """General check that the namd_equilb_NPT_control_file
    (high temperature to set temp NAMD control file) is written."""
    return namd_control_file_written(job, namd_equilb_NVT_control_file_box_0_name_str)

# checking if the NAMD control file is written for the melt equilb NVT run
@Project.label
@flow.with_job
def part_2a_namd_equilb_NVT_box_1_control_file_written(job):
    """General check that the namd_equilb_NPT_control_file
    (high temperature to set temp NAMD control file) is written."""
    return namd_control_file_written(job, namd_equilb_NVT_control_file_box_1_name_str)

# checking if the GOMC control file is written for the equilb run with the selected ensemble
@Project.label
@flow.with_job
def part_2b_gomc_equilb_design_ensemble_control_file_written(job):
    """General check that the gomc_equilb_design_ensemble (run temperature) gomc control file is written."""
    try:
        if (job.sp.solute == "solvent_box"):
            return True
    except:
        return False
    
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                gomc_control_file_written(
                    job,
                    job.doc.gomc_equilb_design_ensemble_dict[
                        str(initial_state_i)
                    ]["output_name_control_file_name"],
                )
            except:
                return False
        return True
    except:
        return False

# checking if the GOMC control file is written for the production run
@Project.label
@flow.with_job
def part_2c_gomc_production_control_file_written(job):
    """General check that the gomc_production_control_file (run temperature) is written."""
    try:
        if (job.sp.solute == "solvent_box"):
            return True
    except:
        return False
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                return gomc_control_file_written(
                    job,
                    job.doc.gomc_production_run_ensemble_dict[
                        str(initial_state_i)
                    ]["output_name_control_file_name"],
                )
            except:
                return False
        return True
    except:
        return False

# ******************************************************
# ******************************************************
# check if GOMC control file was written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC simulations started (start)
# ******************************************************
# ******************************************************
# function for checking if GOMC simulations are started
def gomc_simulation_started(job, control_filename_str):
    """General check to see if the gomc simulation is started."""
    output_started_bool = False
    if job.isfile("out_{}.dat".format(control_filename_str)) and job.isfile(
        "{}_merged.psf".format(control_filename_str)
    ):
        output_started_bool = True

    return output_started_bool

# function for checking if NAMD simulations are started
def namd_simulation_started(job, control_filename_str):
    """General check to see if the namd simulation is started."""
    output_started_bool = False
    if job.isfile("out_{}.dat".format(control_filename_str)) and job.isfile(
        "{}.restart.xsc".format(control_filename_str)
    ):
        output_started_bool = True

    return output_started_bool


# check if melt equilb_NVT namd run is started
@Project.label
@flow.with_job
def part_3a_output_namd_equilb_NVT_box_0_started(job):
    """Check to see if the namd_equilb_NPT_control_file is started
    (high temperature to set temperature in NAMD control file)."""
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_simulation_started(ewald_job, namd_equilb_NVT_control_file_box_0_name_str)
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_simulation_started(ewald_job, namd_equilb_NVT_control_file_box_0_name_str)

    return namd_simulation_started(job, namd_equilb_NVT_control_file_box_0_name_str)


# check if melt equilb_NVT namd run is started
@Project.label
@flow.with_job
def part_3a_output_namd_equilb_NVT_box_1_started(job):
    """Check to see if the namd_equilb_NPT_control_file is started
    (high temperature to set temperature in NAMD control file)."""
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_simulation_started(ewald_job, namd_equilb_NVT_control_file_box_1_name_str)
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_simulation_started(ewald_job, namd_equilb_NVT_control_file_box_1_name_str)

    return namd_simulation_started(job, namd_equilb_NVT_control_file_box_1_name_str)


# check if equilb_with design ensemble GOMC run is started
@Project.label
@flow.with_job
def part_3b_output_gomc_equilb_design_ensemble_started(job):
    """Check to see if the gomc_equilb_design_ensemble simulation is started (set temperature)."""
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                if job.isfile(
                    "out_{}.dat".format(
                        job.doc.gomc_equilb_design_ensemble_dict[
                            str(initial_state_i)
                        ]["output_name_control_file_name"]
                    )
                ):
                    gomc_simulation_started(
                        job,
                        job.doc.gomc_equilb_design_ensemble_dict[
                            str(initial_state_i)
                        ]["output_name_control_file_name"],
                    )

                else:
                    return False
            except:
                return False

        return True
    except:
        return False

# check if equilb_with design ensemble GOMC run is started
@Project.label
@flow.with_job
def part_3b_output_gomc_calibration_started(job):
    """Check to see if the gomc_calibration simulation is started (set temperature)."""
    try:
        ewald_sp = job.statepoint()
        ewald_sp['replica_number_int']=0
        ewald_sp['electrostatic_method']="Wolf"
        ewald_sp['wolf_potential']="Calibrator"
        ewald_sp['wolf_model']="Calibrator"
        jobs = list(pr.find_jobs(ewald_sp))
        for ewald_job in jobs:
            if ewald_job.isfile(
                "Wolf_Calibration_VLUGTWINTRACUTOFF_DSF_BOX_0_wolf_calibration.dat"
            ):
                return True
            else:
                return False


    except:
        return False

# check if equilb_with design ensemble GOMC run is started
@Project.label
@flow.with_job
def part_3b_output_gomc_sseq_started(job):
    """Check to see if the gomc_calibration simulation is started (set temperature)."""
    Single_state_gomc_eq_control_file_name = "single_state_eq"
#This will cause Ewald sims to wait for Wolf calibration to complete.
        #This will cause Ewald sims to wait for Wolf calibration to complete.
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                ewald_job.isfile(f"out_{Single_state_gomc_eq_control_file_name}.dat")
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return ewald_job.isfile(f"out_{Single_state_gomc_eq_control_file_name}.dat")
    else:
        return job.isfile(f"out_{Single_state_gomc_eq_control_file_name}.dat")


# check if equilb_with design ensemble GOMC run is started
@Project.label
@flow.with_job
def part_3b_output_gomc_wolf_sanity_started(job):
    """Check to see if the gomc_calibration simulation is started (set temperature)."""
    wolf_sanity_eq_control_file_name = "wolf_sanity"
    try:
        #This will cause Ewald sims to wait for Wolf calibration to complete.
        if(job.sp.electrostatic_method == "Ewald"):
            wolf_sp = job.statepoint()
            wolf_sp['electrostatic_method']="Wolf"
            jobs = list(pr.find_jobs(wolf_sp))
            for ewald_job in jobs:
                if ewald_job.isfile(f"out_{wolf_sanity_eq_control_file_name}.dat"):
                    return True
                else:
                    return False


        if job.isfile(f"out_{wolf_sanity_eq_control_file_name}.dat"):
            return True
        else:
            return False
    except:
        return False

        return True

# check if production GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@flow.with_job
def part_part_3c_output_gomc_production_run_started(job):
    """Check to see if the gomc production run simulation is started (set temperature)."""
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                if job.isfile(
                    "out_{}.dat".format(
                        job.doc.gomc_production_run_ensemble_dict[
                            str(initial_state_i)
                        ]["output_name_control_file_name"]
                    )
                ):
                    gomc_simulation_started(
                        job,
                        job.doc.gomc_production_run_ensemble_dict[
                            str(initial_state_i)
                        ]["output_name_control_file_name"],
                    )
                else:
                    return False
            except:
                return False
        return True
    except:
        return False

# ******************************************************
# ******************************************************
# check if GOMC simulations started (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC and NAMD simulation are completed properly (start)
# ******************************************************
# ******************************************************
# function for checking if GOMC simulations are completed properly
def gomc_sim_completed_properly(job, control_filename_str):
    """General check to see if the gomc simulation was completed properly."""
    job_run_properly_bool = False
    output_log_file = "out_{}.dat".format(control_filename_str)
    if job.isfile(output_log_file):
        with open(job.fn(f"{output_log_file}"), "r") as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "Move" in line:
                    split_move_line = line.split()
                    if (
                        split_move_line[0] == "Move"
                        and split_move_line[1] == "Type"
                        and split_move_line[2] == "Mol."
                        and split_move_line[3] == "Kind"
                    ):
                        job_run_properly_bool = True
    else:
        job_run_properly_bool = False

    return job_run_properly_bool

# function for checking if NAMD simulations are completed properly
def namd_sim_completed_properly(job, control_filename_str):
    """General check to see if the namd simulation was completed properly."""
    job_run_properly_bool = False

    try:
        if (job.sp.electrostatic_method == "Wolf"):
            output_log_file = job.doc.path_to_namd_console
        else:
            output_log_file = "out_{}.dat".format(control_filename_str)
        if job.isfile(output_log_file):
            with open(job.fn(f"{output_log_file}"), "r") as fp:
                out_namd = fp.readlines()
                for i, line in enumerate(out_namd):
                    if "WallClock:" in line:
                        split_move_line = line.split()
                        if (split_move_line[0] == "WallClock:"
                                and split_move_line[2] == "CPUTime:"
                                and split_move_line[4] == "Memory:"
                        ):
                            job_run_properly_bool = True
        else:
            job_run_properly_bool = False
    except:
        job_run_properly_bool = False
    return job_run_properly_bool

# check if melt equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4a_job_namd_equilb_NVT_box_0_completed_properly(job):
    """Check to see if the  namd_equilb_NPT_control_file was completed properly
    (high temperature to set temperature NAMD control file)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_0_name_str)
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_0_name_str)
    elif (job.sp.replica_number_int == 0):
        return namd_sim_completed_properly(job, namd_equilb_NVT_control_file_box_0_name_str)
    else:
        ewald_sp = job.statepoint()
        ewald_sp['replica_number_int']=0
        jobs = list(pr.find_jobs(ewald_sp))
        for ewald_job in jobs:
            return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_0_name_str)

# check if melt equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4a_job_namd_equilb_NVT_box_1_completed_properly(job):
    """Check to see if the  namd_equilb_NPT_control_file was completed properly
    (high temperature to set temperature NAMD control file)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_1_name_str)
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_1_name_str)
    elif (job.sp.replica_number_int == 0):
        return namd_sim_completed_properly(job, namd_equilb_NVT_control_file_box_1_name_str)
    else:
        ewald_sp = job.statepoint()
        ewald_sp['replica_number_int']=0
        jobs = list(pr.find_jobs(ewald_sp))
        for ewald_job in jobs:
            return namd_sim_completed_properly(ewald_job, namd_equilb_NVT_control_file_box_1_name_str)

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_calibration_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        ewald_sp = job.statepoint()
        ewald_sp['electrostatic_method']="Wolf"
        ewald_sp['solute']="solvent_box"
        ewald_sp['wolf_model']="Calibrator"        
        ewald_sp['wolf_potential']="Calibrator"
        ewald_sp['replica_number_int']=0
        jobs = list(pr.find_jobs(ewald_sp))
        for ewald_job in jobs:
            control_file_name_str = "wolf_calibration"
            if gomc_sim_completed_properly(
                ewald_job,
                control_file_name_str,
            ) is False:
                return False
            elif ewald_job.isfile(
                "Wolf_Calibration_VLUGTWINTRACUTOFF_DSF_BOX_0_wolf_calibration.dat"
            ):
                return True
            else:
                return False
    except:
        return False


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_sseq_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    Single_state_gomc_eq_control_file_name = "single_state_eq"
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    if(job.sp.electrostatic_method == "Wolf"):
        if (job.sp.solute in ["solvent_box"]):
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return gomc_sim_completed_properly(ewald_job, Single_state_gomc_eq_control_file_name)
        else:
            ewald_sp = job.statepoint()
            ewald_sp['electrostatic_method']="Ewald"
            ewald_sp['wolf_model']="Ewald"
            ewald_sp['wolf_potential']="Ewald"
            jobs = list(pr.find_jobs(ewald_sp))
            for ewald_job in jobs:
                return gomc_sim_completed_properly(ewald_job, Single_state_gomc_eq_control_file_name)
    else:
        return gomc_sim_completed_properly(job, Single_state_gomc_eq_control_file_name)

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_wolf_sanity_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    wolf_sanity_control_file_name = "wolf_sanity"
    
    if(job.sp.wolf_model == "Calibrator"):
        return True
    """
    if(job.sp.electrostatic_method == "Ewald"):
        wolf_sp = job.statepoint()
        wolf_sp['electrostatic_method']="Wolf"
        jobs = list(pr.find_jobs(wolf_sp))
        for ewald_job in jobs:
            if gomc_sim_completed_properly(
                ewald_job,
                wolf_sanity_control_file_name,
            ) is False:
                return False
            else:
                return True
    """

    try:
        if gomc_sim_completed_properly(
            job,
            wolf_sanity_control_file_name,
        ) is False:
            #print("gomc_equilb_design_ensemble incomplete state " +  str(initial_state_i))
            return False
        else:
            return True
    except:
        return False


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_wolf_sanity_individual_simulation_averages_completed(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    if(job.sp.wolf_model == "Calibrator"):
        return True
    
    return job.isfile('wolf_sanity_energies_{}.csv'.format(job.id)) and \
        job.isfile('wolf_sanity_full_energies_{}.csv'.format(job.id))



@Project.operation.with_directives(
     {
         "np": 1,
         "ngpu": 0,
         "memory": memory_needed,
         "walltime": walltime_gomc_analysis_hr,
     }
)
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.post(part_4b_wolf_sanity_individual_simulation_averages_completed)
@flow.with_job
def part_4b_wolf_sanity_individual_simulation_averages(job):
    
    import re
    EnRegex = re.compile("ENER_0")
    DensRegex = re.compile("STAT_0")

    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_solute_title = 'solute'  # column title title for temp
    output_column_energy_title = 'total_energy'  # column title title for delta_MBAR
    output_column_energy_std_title = 'total_energy_std'  # column title title for ds_MBAR
    output_column_density_title = 'density'  # column title title for delta_MBAR
    output_column_density_std_title = 'density_std'  # column title title for ds_MBAR


    # get the averages from each individual simulation and write the csv's.
    k_b = 1.9872036E-3  # kcal/mol/K
    temperature = job.sp.production_temperature_K
    k_b_T = temperature * k_b
    dict_of_energies = {}
    dict_of_densities = {}
    dict_of_equilibrated_energies = {}
    dict_of_equilibrated_densities = {}
    dict_of_uncorr_energies = {}
    dict_of_uncorr_densities = {}
    dict_of_full_energies = {}
    dict_of_full_densities = {}
    blk_file = f'out_wolf_sanity.dat'
    steps = []
    energies = []
    densities = []
    with open(blk_file, 'r', encoding='utf8') as f:
        for line in f:
            if EnRegex.match(line):
                try:
                    steps.append(float(line.split()[1]))
                    energies.append(float(line.split()[2]))
                except:
                    print(line)
                    print("An exception occurred") 
            if DensRegex.match(line):
                #print('\n'.join(line.split()[1] for line in f))
                try:
                    if (job.doc.equilibration_ensemble in ["NVT"]):
                        densities.append(float(line.split()[7]))
                    elif (job.doc.equilibration_ensemble in ["NPT"]):
                        densities.append(float(line.split()[8]))      
                    elif (job.doc.equilibration_ensemble in ["GEMC_NVT"]):
                        densities.append(float(line.split()[4]))              
                except:
                    print("An exception occurred") 
    steps_np = np.array(steps)
    energies_np = np.array(energies)
    densities_np = np.array(densities)

    nskip = 10000

    dict_of_full_energies["steps"] = steps_np
    dict_of_full_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = energies_np
    
    df2 = pd.DataFrame.from_dict(dict_of_full_energies)
    df2.to_csv('wolf_sanity_full_energies_{}.csv'.format(job.id), header=True, index=False, sep=' ')
    
    dict_of_full_densities["steps"] = steps_np
    dict_of_full_densities[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = densities_np
    
    df4 = pd.DataFrame.from_dict(dict_of_full_densities)
    df4.to_csv('wolf_sanity_full_densities_{}.csv'.format(job.id), header=True, index=False, sep=' ')

    from pymbar import timeseries
    t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
    A_t_equil = energies_np[t0:]
    A_t_equil_densities = densities_np[t0:]
    A_t_equil_steps = steps_np[t0:]

    dict_of_equilibrated_energies["steps"] = A_t_equil_steps
    dict_of_equilibrated_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = A_t_equil
    dict_of_equilibrated_densities[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = A_t_equil_densities
    dict_of_equilibrated_densities["steps"] = A_t_equil_steps

    dfUC1 = pd.DataFrame.from_dict(dict_of_equilibrated_energies)
    dfUC1.to_csv('wolf_sanity_equilibrated_energies_{}.csv'.format(job.id), header=True, index=False, sep=' ')
    dfUC2 = pd.DataFrame.from_dict(dict_of_equilibrated_densities)
    dfUC2.to_csv('wolf_sanity_equilibrated_densities_{}.csv'.format(job.id), header=True, index=False, sep=' ')

    indices = timeseries.subsampleCorrelatedData(A_t_equil, g=g)
    steps_np = A_t_equil_steps[indices]
    energies_np = A_t_equil[indices]
    densities_np = A_t_equil_densities[indices]

    print("Num uncorrelated equilibrated energy samples",np.shape(energies_np)[0])
    dict_of_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_mean'] = [energies_np.mean()]
    dict_of_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_std'] = [energies_np.std()]
    
    dict_of_uncorr_energies["steps"] = steps_np
    dict_of_uncorr_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = energies_np

    df1 = pd.DataFrame.from_dict(dict_of_energies)
    df1.to_csv('wolf_sanity_energies_{}.csv'.format(job.id))
    
    df2 = pd.DataFrame.from_dict(dict_of_uncorr_energies)
    df2.to_csv('wolf_sanity_uncorr_energies_{}.csv'.format(job.id), header=True, index=False, sep=' ')
    #df2.to_csv('wolf_sanity_full_energies_{}.csv'.format(job.id), header=False, index=False, sep=' ')
    
    #print(densities_np.mean())
    dict_of_densities[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_mean'] = [densities_np.mean()]
    dict_of_densities[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_std'] = [densities_np.std()]

    dict_of_uncorr_densities["steps"] = steps_np
    dict_of_uncorr_densities[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = densities_np

    df3 = pd.DataFrame.from_dict(dict_of_densities)
    df3.to_csv('wolf_sanity_densities_{}.csv'.format(job.id))

    df4 = pd.DataFrame.from_dict(dict_of_uncorr_densities)
    #df4.to_csv('wolf_sanity_full_densities_{}.csv'.format(job.id), header=False, index=False, sep=' ')
    df4.to_csv('wolf_sanity_uncorr_densities_{}.csv'.format(job.id), header=True, index=False, sep=' ')

@Project.label
@flow.with_job
def part_4b_wolf_sanity_analysis_completed(job):
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"        
    ewald_sp['wolf_potential']="Calibrator"   
    ewald_sp['solute']="solvent_box"   
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    try:
        for ewald_job in jobs:
            if (ewald_job.isfile("wolf_statistics.csv")):
                job.doc.winningWolfPotential = ewald_job.doc.winningWolfPotential
                job.doc.winningWolfModel = ewald_job.doc.winningWolfModel
                return True
            else:
                return False
    except:
        return False

@Project.label
@flow.with_job
def part_4b_wolf_sanity_histograms_created(job):
    df1 = pd.DataFrame()
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"        
    ewald_sp['wolf_potential']="Calibrator"   
    ewald_sp['solute']="solvent_box"   
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    try:
        for ewald_job in jobs:
            if (ewald_job.isfile("wolf_sanity_all_energies.csv")):
                df1 = pd.read_csv (ewald_job.fn('wolf_sanity_all_energies.csv'), sep=',', header=0, na_values='NaN', index_col=0)
            else:
                return False
    except:
        return False

    colList = df1.columns.tolist()
    colList.remove("Ewald_Ewald")
    colList.remove("steps")
    try:
        for ewald_job in jobs:
            for col, col_i in zip(colList, range(0, len(colList))):
                try:
                    if (ewald_job.isfile("PotentialEnergyDistribution_Ewald_vs_{}.png".format(col))):
                        continue
                    else:
                        return False
                except:
                    return False
        return True
    except:
        return False
@Project.label
@flow.with_job
def part_4b_is_winning_wolf_model_or_ewald(job):
    try:
        if (job.sp.electrostatic_method == "Ewald"):
            return True
        elif (job.sp.wolf_model == job.doc.winningWolfModel \
            and job.sp.wolf_potential == job.doc.winningWolfPotential):
            return True
        else:
            return False
    except:
        return False
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Calibrator")
@Project.pre(lambda j: j.sp.wolf_model == "Calibrator")
@Project.pre(lambda j: j.sp.solute == "solvent_box")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(lambda *jobs: all(part_4b_wolf_sanity_individual_simulation_averages_completed(j)
                               for j in jobs[0]._project))
@Project.post(part_4b_wolf_sanity_analysis_completed)
@flow.with_job
def part_4b_wolf_sanity_analysis(job):
    df1 = pd.DataFrame()
    df3 = pd.DataFrame()
    df5 = pd.DataFrame()


    jobs = list(pr.find_jobs({"replica_number_int": 0}))
    print(jobs)
    for other_job in jobs:
            print("reading wolf_sanity_equilibrated_energies_{}.csv".format(other_job.id))
            try:
                df6 = pd.read_csv (other_job.fn('wolf_sanity_equilibrated_energies_{}.csv'.format(other_job.id)), sep=' ')
                #print(df2)
                if (df5.empty):
                    df5 = df6
                else:
                    #df1 = df1.merge(df2, on="steps")
                    df5 = pd.merge(df5, df6, on='steps', how='outer')
            except:
                print("failed to read dataframe")
                


    jobs = list(pr.find_jobs({"replica_number_int": 0}))
    print(jobs)
    for other_job in jobs:
            print("reading wolf_sanity_uncorr_energies_{}.csv".format(other_job.id))
            try:
                df2 = pd.read_csv (other_job.fn('wolf_sanity_uncorr_energies_{}.csv'.format(other_job.id)), sep=' ')
                #print(df2)
                if (df1.empty):
                    df1 = df2
                else:
                    #df1 = df1.merge(df2, on="steps")
                    df1 = pd.merge(df1, df2, on='steps', how='outer')
            except:
                print("failed to read dataframe")
                
    for other_job in jobs:
            print("reading wolf_sanity_full_energies_{}.csv".format(other_job.id))
            try:
                df4 = pd.read_csv (other_job.fn('wolf_sanity_full_energies_{}.csv'.format(other_job.id)), sep=' ')
                #print(df2)
                if (df3.empty):
                    df3 = df4
                else:
                    #df1 = df1.merge(df2, on="steps")
                    df3 = pd.merge(df3, df4, on='steps', how='outer')
            except:
                print("failed to read dataframe")
        
    print(df1)
    df1.to_csv('wolf_sanity_uncorr_energies.csv')
    df3.to_csv('wolf_sanity_all_energies.csv')
    df5.to_csv('wolf_sanity_equilibrated_energies.csv')


    statistics_equilibrated = pd.DataFrame()
    import scipy
    from scipy.stats import ttest_ind
    from scipy.spatial.distance import jensenshannon
    listOfWolfMethods = list(df5.columns.values.tolist())
    print(listOfWolfMethods)
    listOfWolfMethods.remove("steps")
    print(listOfWolfMethods)
    ref_mean = df5["Ewald_Ewald"].mean()
    for method in listOfWolfMethods:
        print("Comparing statistical identicallness of Ewald and", method)
        welchs_output = scipy.stats.ttest_ind(df5["Ewald_Ewald"], df5[method], equal_var=False, nan_policy='omit')
        statistics_equilibrated[method] = [df5[method].mean(), df5[method].std(),(df5[method].mean()-ref_mean)/ref_mean, welchs_output[0], welchs_output[1]]

    # Change the row indexes
    statistics_equilibrated.index = ['mean', 'std', 'relative_error', 't-statistic', 'p-value']   
    statistics_equilibrated = statistics_equilibrated.T.sort_values('p-value', ascending=False).T
    statistics_equilibrated.to_csv('wolf_statistics_equilibrated.csv', sep = ' ', )

    statistics = pd.DataFrame()
    listOfWolfMethods = list(df1.columns.values.tolist())
    print(listOfWolfMethods)
    listOfWolfMethods.remove("steps")
    print(listOfWolfMethods)
    ref_mean = df1["Ewald_Ewald"].mean()
    for method in listOfWolfMethods:
        print("Comparing statistical identicallness of Ewald and", method)
        welchs_output = scipy.stats.ttest_ind(df1["Ewald_Ewald"], df1[method], equal_var=False, nan_policy='omit')
        statistics[method] = [df1[method].mean(), df1[method].std(),(df1[method].mean()-ref_mean)/ref_mean, welchs_output[0], welchs_output[1]]

    # Change the row indexes
    statistics.index = ['mean', 'std', 'relative_error', 't-statistic', 'p-value']   
    statistics = statistics.T.sort_values('p-value', ascending=False).T
    statistics.to_csv('wolf_statistics.csv', sep = ' ', )

    job.doc.winningWolfModel = (statistics.columns[1]).split("_")[0]
    job.doc.winningWolfPotential = (statistics.columns[1]).split("_")[1]
    print(statistics)



# ******************************************************
# ******************************************************
# data analysis - get the average data from each individual simulation (start)
# ******************************************************
# ******************************************************

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
#@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(part_4b_job_gomc_calibration_completed_properly)
@flow.with_job
def part_4b_job_gomc_wolf_parameters_found(job):
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"        
    ewald_sp['wolf_potential']="Calibrator"
    ewald_sp['solute']="solvent_box"   
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    for ewald_job in jobs:
        if (not ewald_job.isfile("bestWolfParameters.pickle")):
            return False
        else:
            return True
        
# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
#@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(part_4b_job_gomc_calibration_completed_properly)
@flow.with_job
def part_4b_job_gomc_all_surface_plot_created(job):
    import re
    regex = re.compile("(.*)_all_surfaces.html")
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"        
    ewald_sp['wolf_potential']="Calibrator"
    ewald_sp['solute']="solvent_box"   
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    for ewald_job in jobs:
        for root, dirs, files in os.walk(ewald_job.fn("")):
            for file in files:
                if regex.match(file):
                    return True
    return False
        
# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
# For some reason this is failing on all but replica 0..

@Project.label
@flow.with_job
@Project.pre(part_2a_wolf_sanity_control_file_written)
@Project.pre(part_4b_job_gomc_wolf_parameters_found)
def part_4b_job_gomc_wolf_parameters_appended(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    import re
    regex = re.compile("(\w+?)_initial_state_(\w+?).conf")
    success = True
    atLeastOneMatchExists = False
    if (job.sp.electrostatic_method == "Ewald"):
        return True
    """
    if (job.sp.electrostatic_method == "Ewald"):
        ewald_sp = job.statepoint()
        ewald_sp['electrostatic_method']="Wolf"
        jobs = list(pr.find_jobs(ewald_sp))
        for ewald_job in jobs:
            for root, dirs, files in os.walk(ewald_job.fn("")):
                for file in files:
                    if regex.match(file):
                        atLeastOneMatchExists = True
                        with open(ewald_job.fn(file), "r") as openedfile:
                            last_line = openedfile.readlines()[-1]
                        if ("RcutCoulomb" in last_line):
                            continue
                        else:
                            success = success and False
    """
    if(not job.isfile("wolf_sanity.conf")):
        return False
    regex = re.compile("wolf_sanity.conf")
    for root, dirs, files in os.walk(job.fn("")):
        for file in files:
            if regex.match(file):
                atLeastOneMatchExists = True
                with open(file, "r") as openedfile:
                    last_line = openedfile.readlines()[-1]
                if ("RcutCoulomb" in last_line):
                    continue
                else:
                    success = success and False

    return success and atLeastOneMatchExists

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.pre(lambda j: j.sp.wolf_model != "Calibrator" and j.sp.electrostatic_method == "Wolf")
@Project.pre(part_2b_gomc_equilb_design_ensemble_control_file_written)
@Project.pre(part_2c_gomc_production_control_file_written)
@Project.pre(part_2a_wolf_sanity_control_file_written)
@Project.pre(part_4b_job_gomc_wolf_parameters_found)
@Project.post(part_4b_job_gomc_wolf_parameters_appended)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def part_4b_job_gomc_append_wolf_parameters(job):
    import pickle as pickle
    testEachWolf = True
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"
    ewald_sp['wolf_potential']="Calibrator"
    ewald_sp['solute']="solvent_box"
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    winningWolf = {}
    for ewald_job in jobs:
        if (testEachWolf):
            with open(ewald_job.fn("bestWolfParameters.pickle"), 'rb') as handle:
                winningWolf = pickle.load(handle)
        else:
            try:
                import pickle as pickle
                import re
                with open(ewald_job.fn("bestWolfParameters.pickle"), 'rb') as handle:
                    model2BestWolfAlphaRCut = pickle.load(handle)
                
                bestModel = ""
                smallestRelErr = 1.0
                largestRelErr = 0
                bestRCut = 0
                bestAlpha = 0

                #print("Replica :", job.sp.replica)
                for model in model2BestWolfAlphaRCut:
                    print("Model :", model)
                    print("RelErr :",  model2BestWolfAlphaRCut[model]['GD_relerr'])
                    print("RCut :", model2BestWolfAlphaRCut[model]['GD_rcut'])
                    print("Alpha :", model2BestWolfAlphaRCut[model]['GD_alpha'])
                    if (model2BestWolfAlphaRCut[model]['GD_relerr']  < smallestRelErr):
                        bestModel = model
                        smallestRelErr = model2BestWolfAlphaRCut[model]['GD_relerr']   
                        bestRCut =  model2BestWolfAlphaRCut[model]['GD_rcut']  
                        bestAlpha =  model2BestWolfAlphaRCut[model]['GD_alpha']  
                    if (model2BestWolfAlphaRCut[model]['GD_relerr']  > largestRelErr):
                        worstModel = model
                        largestRelErr = model2BestWolfAlphaRCut[model]['GD_relerr']   
                        worstRCut =  model2BestWolfAlphaRCut[model]['GD_rcut']  
                        worstAlpha =  model2BestWolfAlphaRCut[model]['GD_alpha']  
                print("worstModel :", worstModel)
                print("largestRelErr :", largestRelErr)
                print("worstRCut :", worstRCut)
                print("worstAlpha :", worstAlpha)

                print("bestModel :", bestModel)
                print("smallestRelErr :", smallestRelErr)
                print("bestRCut :", bestRCut)
                print("bestAlpha :", bestAlpha)

                winningWolf = {"WolfKind": bestModel[0], "Potential": bestModel[1], "RCutCoul": bestRCut,
                "Alpha":bestAlpha}
                with open("winningWolfParameters.pickle", 'wb') as handle:
                    pickle.dump(winningWolf, handle, protocol=pickle.HIGHEST_PROTOCOL)

            except:
                return False

    import re
    regex = re.compile("(\w+?)_initial_state_(\w+?).conf")
    if (job.doc.equilibration_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]):  
        box_list = ["0", "1"]
    else:
        box_list = ["0"]
    for root, dirs, files in os.walk(job.fn("")):
        for file in files:
            if regex.match(file):
                with open(file, "a") as myfile:
                    if (testEachWolf):
                        defWolfLine = "Wolf\tTrue\n"
                        myfile.write(defWolfLine)
                        defPotLine = "WolfPotential\t{pot}\n".format(pot=job.sp.wolf_potential)
                        myfile.write(defPotLine)
                        defKindLine = "WolfKind\t{kind}\n".format(kind=job.sp.wolf_model)
                        myfile.write(defKindLine)
                        for box in box_list:
                            defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=winningWolf[(job.sp.wolf_model, job.sp.wolf_potential, box)]["GD_alpha"])
                            myfile.write(defAlphaLine)
                            defRCutLine = "RcutCoulomb\t{box}\t{val}\n".format(box=box, val=winningWolf[(job.sp.wolf_model, job.sp.wolf_potential, box)]["GD_rcut"])
                            myfile.write(defRCutLine)
                    else:
                        defWolfLine = "Wolf\tTrue\n"
                        myfile.write(defWolfLine)
                        defPotLine = "WolfPotential\t{pot}\n".format(pot=winningWolf["Potential"])
                        myfile.write(defPotLine)
                        defKindLine = "WolfKind\t{kind}\n".format(kind=winningWolf["WolfKind"])
                        myfile.write(defKindLine)
                        defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=winningWolf["Alpha"])
                        myfile.write(defAlphaLine)
                        defRCutLine = "RcutCoulomb\t{box}\t{val}\n".format(box=box, val=winningWolf["RCutCoul"])
                        myfile.write(defRCutLine)


    regex = re.compile("wolf_sanity.conf")
    box = "0"
    for root, dirs, files in os.walk(job.fn("")):
        for file in files:
            if regex.match(file):
                with open(file, "a") as myfile:
                    if (testEachWolf):
                        defWolfLine = "Wolf\tTrue\n"
                        myfile.write(defWolfLine)
                        defPotLine = "WolfPotential\t{pot}\n".format(pot=job.sp.wolf_potential)
                        myfile.write(defPotLine)
                        defKindLine = "WolfKind\t{kind}\n".format(kind=job.sp.wolf_model)
                        myfile.write(defKindLine)
                        for box in box_list:
                            defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=winningWolf[(job.sp.wolf_model, job.sp.wolf_potential, box)]["GD_alpha"])
                            myfile.write(defAlphaLine)
                            defRCutLine = "RcutCoulomb\t{box}\t{val}\n".format(box=box, val=winningWolf[(job.sp.wolf_model, job.sp.wolf_potential, box)]["GD_rcut"])
                            myfile.write(defRCutLine)
                    else:
                        defWolfLine = "Wolf\tTrue\n"
                        myfile.write(defWolfLine)
                        defPotLine = "WolfPotential\t{pot}\n".format(pot=winningWolf["Potential"])
                        myfile.write(defPotLine)
                        defKindLine = "WolfKind\t{kind}\n".format(kind=winningWolf["WolfKind"])
                        myfile.write(defKindLine)
                        defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=box, val=winningWolf["Alpha"])
                        myfile.write(defAlphaLine)
                        defRCutLine = "RcutCoulomb\t{box}\t{val}\n".format(box=box, val=winningWolf["RCutCoul"])
                        myfile.write(defRCutLine)

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_equilb_design_ensemble_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                filename_4b_iter = job.doc.gomc_equilb_design_ensemble_dict[
                    str(initial_state_i)
                ]["output_name_control_file_name"]

                if gomc_sim_completed_properly(
                    job,
                    filename_4b_iter,
                ) is False:
                    #print("gomc_equilb_design_ensemble incomplete state " +  str(initial_state_i) + " " + job.fn(""))
                    return False
            except:
                return False
        return True
    except:
        return False

# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4c_job_production_run_completed_properly(job):
    """Check to see if the gomc production run simulation was completed properly (set temperature)."""
    try:
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                filename_4c_iter = job.doc.gomc_production_run_ensemble_dict[
                    str(initial_state_i)
                ]["output_name_control_file_name"]
                if gomc_sim_completed_properly(
                    job,
                    filename_4c_iter,
                ) is False:
                    #print("Isn't finished ",filename_4c_iter)
                    return False

                # check specifically for the FE files
                if job.isfile(f'Free_Energy_BOX_0_{filename_4c_iter}.dat') is False:
                    #print("Isn't finished ",f'Free_Energy_BOX_0_{filename_4c_iter}.dat')
                    return False

            except:
                return False
        return True
    except:
        return False


# check if analysis is done for the individual replicates wrote the gomc files
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.label
@flow.with_job
def part_5a_analysis_individual_simulation_averages_completed(job):
    """Check that the individual simulation averages files are written ."""
    file_written_bool = False
    if (
        job.isfile(
            f"{output_replicate_txt_file_name_liq}"
        )
        and job.isfile(
            f"{output_replicate_txt_file_name_vap}"
        )
    ):
        file_written_bool = True

    return file_written_bool


# check if analysis for averages of all the replicates is completed
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.label
def part_5b_analysis_replica_averages_completed(*jobs):
    """Check that the individual simulation averages files are written ."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{output_avg_std_of_replicates_txt_file_name_liq}"
            )
            and job.isfile(
                f"../../analysis/{output_avg_std_of_replicates_txt_file_name_vap}"
            )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return all_file_written_bool_pass

# check if analysis for critical points is completed
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5b_analysis_replica_averages_completed)
@Project.label
def part_5c_analysis_critical_and_boiling_points_replicate_data_completed(*jobs):
    """Check that the critical and boiling point replicate file is written ."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{output_critical_data_replicate_txt_file_name}"
            ) \
                and
            job.isfile(
                f"../../analysis/{output_boiling_data_replicate_txt_file_name}"
            )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return file_written_bool

# check if analysis for critical points is completed
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5b_analysis_replica_averages_completed)
@Project.label
def part_5d_analysis_critical_and_boiling_points_avg_std_data_completed(*jobs):
    """Check that the avg and std dev critical and boiling point data file is written ."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{output_critical_data_avg_std_of_replicates_txt_file_name}"
            ) \
                and
                job.isfile(
                    f"../../analysis/{output_boiling_data_avg_std_of_replicates_txt_file_name}"
                )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return file_written_bool


# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# build system, with option to write the force field (force field (FF)), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (start)
# ******************************************************
# build system
def build_charmm(job, write_files=True):
    """Build the Charmm object and potentially write the pdb, psd, and force field (FF) files."""
    print("#**********************")
    print("Started: GOMC Charmm Object")
    print("#**********************")
    mbuild_box_seed_no = job.doc.replica_number_int

    # set free energy data in doc
    # Free energy calcs
    # lamda generator
    # get the paths to the smiles or mol2 files

    smiles_or_mol2 = get_molecule_path(smiles_or_mol2_name_to_value_dict[job.sp.solvent][job.sp.forcefield])

    solvent = mb.load(smiles_or_mol2[1],
                      smiles=smiles_or_mol2[0]
                      )
    solvent.name = job.doc.solvent

    solvent_ff = get_ff_path(forcefield_residue_to_ff_filename_dict[job.sp.solvent][job.sp.forcefield])

    # only put the FF molecules in the simulation in the dictionaly input into the Chamm object.
    minimal_forcefield_dict = {solvent.name: solvent_ff
                            }

    #solute.energy_minimize(forcefield=forcefield_dict[job.sp.solute], steps=10 ** 5)

    # for trappe, currently unused'
    if (job.sp.forcefield == "TRAPPE"):
        bead_to_atom_name_dict = { '_CH3':'C', '_CH2':'C',  'O':'O', 'H':'H'}
    else:
        bead_to_atom_name_dict = None

    residues_list = [solvent.name]
    print("residues_list  = " +str(residues_list ))

    #if job.doc.solvent in ["TIP4", "TIP3"]:
    gomc_fix_bonds_angles_residues_list = [solvent.name]
    #else:
    #    gomc_fix_bonds_angles_residues_list  = None
    print('Running: filling liquid box')
    box_0 = mb.fill_box(compound=[solvent],
                        n_compounds=[job.doc.N_liquid_solvent],
                        box=[u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                            ],
                        seed=mbuild_box_seed_no
                        )
    print('Completed: filling liquid box')

    print('Running: filling vapor box')
    box_1 = mb.fill_box(compound=[solvent],
                        n_compounds=[job.doc.N_vapor_solvent],
                        box=[u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            ],
                        seed=mbuild_box_seed_no
                        )
    print('Completed: filling vapor box')
    """
    angstrom3 = (u.angstrom * u.angstrom * u.angstrom)
    cm3 = (u.cm * u.cm * u.cm)
    job.doc.volume = ((job.doc.vap_box_lengths_ang * u.angstrom) * (job.doc.vap_box_lengths_ang * u.angstrom) * (job.doc.vap_box_lengths_ang * u.angstrom)).to(cm3)

    from scipy import constants
    molar_mass_of_solvent = 18.01528 * u.mol
    job.doc.N_vapor_solvent = int((constants.Avogadro * job.sp.vapor_density * job.doc.volume )/ molar_mass_of_solvent)
    

    print('Running: filling vapor box : ', job.doc.N_vapor_solvent)
    box_1 = mb.fill_box(compound=[solvent],
                        n_compounds=job.doc.N_vapor_solvent,
                        box=[u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("nm"),
                            ],
                        seed=mbuild_box_seed_no
                        )
    print('Completed: filling vapor box')
    """

    print('Running: GOMC FF file, and the psf and pdb files')
    gomc_charmm = mf_charmm.Charmm(
        box_0,
        mosdef_structure_box_0_name_str,
        structure_box_1=box_1,
        filename_box_1=mosdef_structure_box_1_name_str,
        ff_filename=  gomc_ff_filename_str,
        forcefield_selection=minimal_forcefield_dict,
        residues=residues_list,
        bead_to_atom_name_dict=bead_to_atom_name_dict,
        gomc_fix_bonds_angles=gomc_fix_bonds_angles_residues_list,
    )

    print('Running: namd_charmm')
    namd_charmm = mf_charmm.Charmm(
        box_0,
        mosdef_structure_box_0_name_str,
        structure_box_1=None,
        filename_box_1=None,
        ff_filename= namd_ff_filename_str,
        forcefield_selection=minimal_forcefield_dict,
        residues=residues_list,
        bead_to_atom_name_dict=bead_to_atom_name_dict,
        gomc_fix_bonds_angles=None,
    )

    namd_charmm_box_1 = mf_charmm.Charmm(
        box_1,
        mosdef_structure_box_1_name_str,
        structure_box_1=None,
        filename_box_1=None,
        ff_filename= namd_ff_filename_str,
        forcefield_selection=minimal_forcefield_dict,
        residues=residues_list,
        bead_to_atom_name_dict=bead_to_atom_name_dict,
        gomc_fix_bonds_angles=None,
    )

    gomc_charmm.write_inp()


    namd_charmm.write_inp()

    #if (job.sp.electrostatic_method == "Ewald" and job.sp.replica_number_int == 0):
    gomc_charmm.write_psf()
    gomc_charmm.write_pdb()
    namd_charmm.write_psf()
    namd_charmm.write_pdb()

    #namd_charmm_box_1.write_inp()
    #namd_charmm_box_1.write_psf()
    #namd_charmm_box_1.write_pdb()

    """
    template = get_pdb2bincoords_path()
    # Read in the file
    with open(template, 'r') as file :
        filedata = file.read()
    
    # convert water shell to namd bin coords file
    # Replace the target string
    filedata = filedata.replace("PDB_FILE", job.fn(mosdef_structure_box_1_name_str))
    filedata = filedata.replace("PSF_FILE", job.fn(mosdef_structure_box_1_name_str))

    # Write the file out again
    with open(job.fn("create_box_1_namdbin.tcl"), 'w') as file:
        file.write(filedata)

    from vmd import evaltcl
    print("Making solvated sphere namd bin file", job)
    ions = evaltcl("source " + job.fn("create_box_1_namdbin.tcl"))
    ionsList = ions.split()

    template = get_pdb2xsc_path()
    # Read in the file
    with open(template, 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("PDB_FILE", job.fn(mosdef_structure_box_1_name_str))
    filedata = filedata.replace("PSF_FILE", job.fn(mosdef_structure_box_1_name_str))
    filedata = filedata.replace("XSC_FILE", job.fn(mosdef_structure_box_1_name_str))

    # Write the file out again
    with open(job.fn("create_box_1_xsc.tcl"), 'w') as file:
        file.write(filedata)

    print("Making solvated sphere", job)
    ions = evaltcl("source " + job.fn("create_box_1_xsc.tcl"))
    ionsList = ions.split()
    """
    print("#**********************")
    print("Completed: GOMC Charmm Object")
    print("#**********************")

    return [namd_charmm, gomc_charmm]


# ******************************************************
# ******************************************************
# build system, with option to write the force field (FF), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (end)
# ******************************************************


# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (start)
# ******************************************************
# ******************************************************
@Project.pre(part_1a_initial_data_input_to_json)
@Project.post(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.post(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.post(part_2b_gomc_equilb_design_ensemble_control_file_written)
@Project.post(part_2c_gomc_production_control_file_written)
@Project.post(mosdef_input_written)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def build_psf_pdb_ff_gomc_conf(job):
    """Build the Charmm object and write the pdb, psd, and force field (FF)
    files for all the simulations in the workspace."""
    [namd_charmm_object_with_files, gomc_charmm_object_with_files] = build_charmm(job, write_files=True)

    namd_restart_pdb_psf_file_box_0_name_str = namd_equilb_NVT_control_file_box_0_name_str
    namd_restart_pdb_psf_file_box_1_name_str = namd_equilb_NVT_control_file_box_1_name_str

    prefix = ""
    if(job.sp.electrostatic_method == "Ewald" and job.sp.replica_number_int == 0):
        prefix = job.fn("")
    else:
        ref_sp = job.statepoint()
        ref_sp['electrostatic_method']="Ewald"
        ref_sp['replica_number_int']=0
        jobs = list(pr.find_jobs(ref_sp))
        for ref_job in jobs:
            prefix = ref_job.fn("")

    Coordinates_box_0 = "{}.pdb".format(
        prefix+mosdef_structure_box_0_name_str
    )
    Structure_box_0 = "{}.psf".format(
        prefix+mosdef_structure_box_0_name_str
    )
    Coordinates_box_1 = "{}.pdb".format(
        prefix+mosdef_structure_box_1_name_str
    )
    Structure_box_1 = "{}.psf".format(
        prefix+mosdef_structure_box_1_name_str
    )
    binCoordinates_box_0 = "{}.restart.coor".format(
        prefix+namd_restart_pdb_psf_file_box_0_name_str
    )
    extendedSystem_box_0 = "{}.restart.xsc".format(
        prefix+namd_restart_pdb_psf_file_box_0_name_str
    )
    binCoordinates_box_1 = "{}.restart.coor".format(
        prefix+namd_restart_pdb_psf_file_box_1_name_str
    )
    extendedSystem_box_1 = "{}.restart.xsc".format(
        prefix+namd_restart_pdb_psf_file_box_1_name_str
    )

    job.doc.path_to_namd_box_0_console =  prefix+f"out_{namd_equilb_NVT_control_file_box_0_name_str}.dat"
    job.doc.path_to_namd_box_1_console =  prefix+f"out_{namd_equilb_NVT_control_file_box_1_name_str}.dat"

    job.doc.path_to_ref_pdb =  Coordinates_box_0
    job.doc.path_to_ref_psf =  Structure_box_0
    job.doc.path_to_ref_binCoordinates =  binCoordinates_box_0
    job.doc.path_to_ref_extendedSystem =  extendedSystem_box_0

    if (job.doc.equilibration_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]):  
        job.doc.path_to_ref_pdb_box_1 =  Coordinates_box_1
        job.doc.path_to_ref_psf_box_1 =  Structure_box_1
        job.doc.path_to_ref_binCoordinates_box_1 =  binCoordinates_box_1
        job.doc.path_to_ref_extendedSystem_box_1 =  extendedSystem_box_1 
    else:
        job.doc.path_to_ref_pdb_box_1 =  None
        job.doc.path_to_ref_psf_box_1 =  None
        job.doc.path_to_ref_binCoordinates_box_1 =  None
        job.doc.path_to_ref_extendedSystem_box_1 =  None         

    Single_state_gomc_eq_control_file_name = "single_state_eq"


    Single_state_gomc_eq_Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Structure_box_0 = "{}_BOX_0_restart.psf".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Coordinates_box_1 = "{}_BOX_1_restart.pdb".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Structure_box_1 = "{}_BOX_1_restart.psf".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_binCoordinates_box_1 = "{}_BOX_1_restart.coor".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_extendedSystem_box_1 = "{}_BOX_1_restart.xsc".format(
        Single_state_gomc_eq_control_file_name
    )

    if (job.sp.electrostatic_method == "Wolf"):
        ref_sp = job.statepoint()
        ref_sp['electrostatic_method']="Ewald"
        if (job.sp.wolf_model == "Calibrator"):
            ref_sp['wolf_model']="Ewald"
            ref_sp['wolf_potential']="Ewald"
        else:
            ref_sp['wolf_model']="Ewald"
            ref_sp['wolf_potential']="Ewald"
        jobs = list(pr.find_jobs(ref_sp))
        for ref_job in jobs:
            #if (ref_job.isfile(f"{Coordinates_box_0}")):
            job.doc.path_to_namd_console =  ref_job.fn(f"out_{namd_equilb_NVT_control_file_box_0_name_str}.dat")
            job.doc.path_to_namd_console =  ref_job.fn(f"out_{namd_equilb_NVT_control_file_box_1_name_str}.dat")

            job.doc.path_to_sseq_pdb =  ref_job.fn(Single_state_gomc_eq_Coordinates_box_0)
            job.doc.path_to_sseq_psf =  ref_job.fn(Single_state_gomc_eq_Structure_box_0)
            job.doc.path_to_sseq_pdb_box_1 =  ref_job.fn(Single_state_gomc_eq_Coordinates_box_1)
            job.doc.path_to_sseq_psf_box_1 =  ref_job.fn(Single_state_gomc_eq_Structure_box_1)
            job.doc.path_to_sseq_binCoordinates =  ref_job.fn(Single_state_gomc_eq_binCoordinates_box_0)
            job.doc.path_to_sseq_extendedSystem =  ref_job.fn(Single_state_gomc_eq_extendedSystem_box_0)
            job.doc.path_to_sseq_binCoordinates_box_1 =  ref_job.fn(Single_state_gomc_eq_binCoordinates_box_1)
            job.doc.path_to_sseq_extendedSystem_box_1 =  ref_job.fn(Single_state_gomc_eq_extendedSystem_box_1)
            job.doc.path_to_sseq_console =  ref_job.fn(f"out_{Single_state_gomc_eq_control_file_name}.dat")
            job.doc.path_to_sseq_checkpoint =  ref_job.fn(f"{Single_state_gomc_eq_control_file_name}_restart.chk")
       
    else:    
        job.doc.path_to_sseq_pdb =  job.fn(Single_state_gomc_eq_Coordinates_box_0)
        job.doc.path_to_sseq_psf =  job.fn(Single_state_gomc_eq_Structure_box_0)
        job.doc.path_to_sseq_pdb_box_1 =  job.fn(Single_state_gomc_eq_Coordinates_box_1)
        job.doc.path_to_sseq_psf_box_1 =  job.fn(Single_state_gomc_eq_Structure_box_1)
        job.doc.path_to_sseq_binCoordinates =  job.fn(Single_state_gomc_eq_binCoordinates_box_0)
        job.doc.path_to_sseq_extendedSystem =  job.fn(Single_state_gomc_eq_extendedSystem_box_0)
        job.doc.path_to_sseq_binCoordinates_box_1 =  job.fn(Single_state_gomc_eq_binCoordinates_box_1)
        job.doc.path_to_sseq_extendedSystem_box_1 =  job.fn(Single_state_gomc_eq_extendedSystem_box_1)
        job.doc.path_to_sseq_console =  job.fn(f"out_{Single_state_gomc_eq_control_file_name}.dat")
        job.doc.path_to_sseq_checkpoint =  job.fn(f"{Single_state_gomc_eq_control_file_name}_restart.chk")

    FreeEnergyCalc = [True, int(gomc_free_energy_output_data_every_X_steps)]
    # This has to be off during calibration
    NoFreeEnergyCalc = [False, int(gomc_free_energy_output_data_every_X_steps)]

    MoleculeType = [job.sp.solute, 1]

    use_ElectroStatics = True

    if (job.sp.forcefield in ["OPLS"]):
        VDWGeometricSigma = True
    else:
        VDWGeometricSigma = False

    
    
    Exclude = "1-4"

    # common variables
    cutoff_style = "VDW"
    if cutoff_style != "VDW":
        raise ValueError("ERROR: this project is only set up for the SWITCH cutoff style for NAMD"
                         "and VDW for GOMC.  Therefore, the cutoff style selected must be VDW. "
                         "Rswitch for namd only so the r_switch_dist_start and "
                         "r_switch_dist_end must be supplied for NAMD. GOMC will then move to VDW "
                         "with the switch dist (r_switch_dist_start) as the cutoff with LRC.")

    production_temperature_K = (job.sp.production_temperature_K * u.K).to_value("K")

    production_pressure_bar = (job.doc.production_pressure_bar * u.bar).to_value("bar")

    box_lengths_liq_ang = [u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       ]

    box_lengths_vap_ang = [u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.vap_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       ]

    seed_no = job.doc.replica_number_int

    namd_template_path_str = os.path.join(project_directory_path, "templates/NAMD_NVT_conf_template.conf")

    if job.doc.solvent in ["TIP3", "SPC", "SPCE", "MSPCE"]:
        namd_uses_water = True
        namd_water_model = 'tip3'
    elif job.doc.solvent in ["TIP4"]:
        namd_uses_water = True
        namd_water_model = 'tip4'
    else:
        namd_uses_water = False
        namd_water_model= None

    # generate the namd file
    # NOTE: the production and melt temps are converted to intergers so they can be ramped down
    # from hot to cool to equilibrate the system.
    generate_namd_equilb_control_file(template_path_filename=namd_template_path_str,
                                      namd_path_conf_filename=namd_equilb_NVT_control_file_box_0_name_str,
                                      namd_path_file_output_names=namd_equilb_NVT_control_file_box_0_name_str,
                                      namd_uses_water=namd_uses_water,
                                      namd_water_model=namd_water_model,
                                      namd_electrostatics_bool=use_ElectroStatics,
                                      namd_vdw_geometric_sigma_bool=VDWGeometricSigma,
                                      namd_psf_path_filename=f"{mosdef_structure_box_0_name_str}.psf",
                                      namd_pdb_path_filename=f"{mosdef_structure_box_0_name_str}.pdb",
                                      namd_ff_path_filename=f"{namd_ff_filename_str}.inp",
                                      namd_production_temp_K= int(production_temperature_K),
                                      namd_production_pressure_bar=production_pressure_bar,
                                      electrostatic_1_4=namd_charmm_object_with_files.coul_1_4,
                                      non_bonded_cutoff=job.doc.Rcut_for_switch_namd_ang,
                                      non_bonded_switch_distance=job.doc.Rcut_ang,
                                      pairlist_distance=job.doc.neighbor_list_dist_namd_ang,
                                      box_lengths=box_lengths_liq_ang,
                                      )

    # generate the namd file
    # NOTE: the production and melt temps are converted to intergers so they can be ramped down
    # from hot to cool to equilibrate the system.
    generate_namd_equilb_control_file(template_path_filename=namd_template_path_str,
                                      namd_path_conf_filename=namd_equilb_NVT_control_file_box_1_name_str,
                                      namd_path_file_output_names=namd_equilb_NVT_control_file_box_1_name_str,
                                      namd_uses_water=namd_uses_water,
                                      namd_water_model=namd_water_model,
                                      namd_electrostatics_bool=use_ElectroStatics,
                                      namd_vdw_geometric_sigma_bool=VDWGeometricSigma,
                                      namd_psf_path_filename=f"{mosdef_structure_box_1_name_str}.psf",
                                      namd_pdb_path_filename=f"{mosdef_structure_box_1_name_str}.pdb",
                                      namd_ff_path_filename=f"{namd_ff_filename_str}.inp",
                                      namd_production_temp_K= int(production_temperature_K),
                                      namd_production_pressure_bar=production_pressure_bar,
                                      electrostatic_1_4=namd_charmm_object_with_files.coul_1_4,
                                      non_bonded_cutoff=job.doc.Rcut_for_switch_namd_ang,
                                      non_bonded_switch_distance=job.doc.Rcut_ang,
                                      pairlist_distance=job.doc.neighbor_list_dist_namd_ang,
                                      box_lengths=box_lengths_vap_ang,
                                      )

    print("#**********************")
    print("Completed: namd_equilb_GEMC GOMC control file writing")
    print("#**********************")
    # ******************************************************
    # namd_equilb_NPT - psf, pdb, force field (FF) file writing and GOMC control file writing  (end)
    # ******************************************************
    MC_steps = int(gomc_steps_equilb_design_ensemble)

    # output all data and calc frequecy
    console_output_true_list_input = [
        True,
        int(gomc_console_output_data_every_X_steps),
    ]    
    # output all data and calc frequecy
    restart_true_list_input = [
        True,
        int(gomc_steps_equilb_design_ensemble//2),
    ]
    output_true_list_input = [
        True,
        int(gomc_output_data_every_X_steps),
    ]
    output_false_list_input = [
        False,
        int(gomc_output_data_every_X_steps),
    ]

    if job.doc.solute in ["ETOH", "ETOH-OPLS", "solvent_box"]:
        useCoul = True
        CBMC_First = (10,)
        CBMC_Nth = (10,)
        CBMC_Ang = (100,)
        CBMC_Dih = (50,)
        if job.doc.equilibration_ensemble in ["NVT"]:
            VolFreq = (0.00,)
            MultiParticleFreq = (None,)
            IntraSwapFreq = (0.0,)
            CrankShaftFreq = (0.1,)
            SwapFreq = (None,)
            DisFreq = (0.50,)
            RotFreq = (0.2,)
            RegrowthFreq = (0.20,)

        elif job.doc.equilibration_ensemble in ["NPT"]:
            VolFreq = (0.01,)
            MultiParticleFreq = (None,)
            IntraSwapFreq = (0.0,)
            CrankShaftFreq = (0.1,)
            SwapFreq = (None,)
            DisFreq = (0.49,)
            RotFreq = (0.2,)
            RegrowthFreq = (0.20,)
            
        elif job.doc.equilibration_ensemble in ["GEMC_NVT"]:
            VolFreq = (0.01,)
            MultiParticleFreq = (0.00,)
            IntraSwapFreq = (0.10,)
            CrankShaftFreq = (0.1,)
            SwapFreq = (0.20,)
            DisFreq = (0.29,)
            RotFreq = (0.20,)
            RegrowthFreq = (0.20,)
            CrankShaftFreq = (0.0,)
        else:
            raise ValueError(
                "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                "in the GOMC control file writer."
            )

    print("#**********************")
    print("Started: equilb NPT NAMD -> NPT GOMC control file writing")
    print("#**********************")
    # output all data and calc frequecy
    restart_true_list_input = [
        True,
        int(MC_steps//2),
    ]
    gomc_control.write_gomc_control_file(
        gomc_charmm_object_with_files,
        Single_state_gomc_eq_control_file_name,
        job.doc.equilibration_ensemble,
        MC_steps,
        production_temperature_K,
        ff_psf_pdb_file_directory=None,
        check_input_files_exist=False,
        Parameters="{}.inp".format(gomc_ff_filename_str),
        Restart=True,
        RestartCheckpoint=True,
        ExpertMode=False,
        Coordinates_box_0=job.doc.path_to_ref_pdb,
        Structure_box_0=job.doc.path_to_ref_psf,
        binCoordinates_box_0=job.doc.path_to_ref_binCoordinates,
        extendedSystem_box_0=job.doc.path_to_ref_extendedSystem,
        binVelocities_box_0=None,
        Coordinates_box_1=job.doc.path_to_ref_pdb_box_1,
        Structure_box_1=job.doc.path_to_ref_psf_box_1,
        binCoordinates_box_1=job.doc.path_to_ref_binCoordinates_box_1,
        extendedSystem_box_1=job.doc.path_to_ref_extendedSystem_box_1,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": production_pressure_bar,
            "Ewald": True,
            "ElectroStatic": use_ElectroStatics,
            "VDWGeometricSigma": VDWGeometricSigma,
            "Rcut": job.doc.Rcut_ang,
            "Exclude": Exclude,
            "VolFreq": VolFreq[-1],
            "MultiParticleFreq": MultiParticleFreq[-1],
            "IntraSwapFreq": IntraSwapFreq[-1],
            "CrankShaftFreq": CrankShaftFreq[-1],
            "SwapFreq": SwapFreq[-1],
            "DisFreq": DisFreq[-1],
            "RotFreq": RotFreq[-1],
            "RegrowthFreq": RegrowthFreq[-1],
            "OutputName": Single_state_gomc_eq_control_file_name,
            "EqSteps": EqSteps,
            "PressureCalc": output_false_list_input,
            "RestartFreq": restart_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_false_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_false_list_input,
            "Potential": cutoff_style,
            "LRC": True,
            "RcutLow": 1.0,
            "RcutCoulomb_box_1" : min((math.floor(job.doc.vap_box_lengths_ang/2.0)-1), 20),
            "CBMC_First": CBMC_First[-1],
            "CBMC_Nth": CBMC_Nth[-1],
            "CBMC_Ang": CBMC_Ang[-1],
            "CBMC_Dih": CBMC_Dih[-1],
        },
    )

    print("#**********************")
    print("Started: equilb NPT NAMD -> NPT GOMC control file writing")
    print("#**********************")

    # output all data and calc frequecy
    restart_true_list_input = [
        True,
        int(Wolf_Sanity_MC_steps//2),
    ]
    print("#**********************")
    print("Started:  Wolf Sanity GOMC control file writing")
    print("#**********************")
    wolf_sanity_control_file_name = "wolf_sanity"
    gomc_control.write_gomc_control_file(
        gomc_charmm_object_with_files,
        wolf_sanity_control_file_name,
        job.doc.equilibration_ensemble,
        Wolf_Sanity_MC_steps,
        production_temperature_K,
        ff_psf_pdb_file_directory=None,
        check_input_files_exist=False,
        Parameters="{}.inp".format(gomc_ff_filename_str),
        Restart=True,
        RestartCheckpoint=True,
        ExpertMode=False,
        Coordinates_box_0=job.doc.path_to_sseq_pdb,
        Structure_box_0=job.doc.path_to_sseq_psf,
        binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
        extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
        binVelocities_box_0=None,
        Coordinates_box_1=job.doc.path_to_sseq_pdb_box_1,
        Structure_box_1=job.doc.path_to_sseq_psf_box_1,
        binCoordinates_box_1=job.doc.path_to_sseq_binCoordinates_box_1,
        extendedSystem_box_1=job.doc.path_to_sseq_extendedSystem_box_1,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": production_pressure_bar,
            "Ewald": job.sp.electrostatic_method == "Ewald",
            "ElectroStatic": use_ElectroStatics,
            "VDWGeometricSigma": VDWGeometricSigma,
            "Rcut": job.doc.Rcut_ang,
            "Exclude": Exclude,
            "VolFreq": VolFreq[-1],
            "MultiParticleFreq": MultiParticleFreq[-1],
            "IntraSwapFreq": IntraSwapFreq[-1],
            "CrankShaftFreq": CrankShaftFreq[-1],
            "SwapFreq": SwapFreq[-1],
            "DisFreq": DisFreq[-1],
            "RotFreq": RotFreq[-1],
            "RegrowthFreq": RegrowthFreq[-1],
            "OutputName": wolf_sanity_control_file_name,
            "EqSteps": EqSteps,
            "PressureCalc": output_true_list_input,
            "RestartFreq": restart_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": console_output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_false_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_false_list_input,
            "Potential": cutoff_style,
            "LRC": True,
            "RcutLow": 1.0,
            "RcutCoulomb_box_1" : min((math.floor(job.doc.vap_box_lengths_ang/2.0)-1), 20) if job.sp.electrostatic_method == "Ewald" else None,
            "CBMC_First": CBMC_First[-1],
            "CBMC_Nth": CBMC_Nth[-1],
            "CBMC_Ang": CBMC_Ang[-1],
            "CBMC_Dih": CBMC_Dih[-1],
        },
    )
    append_checkpoint_line(job, wolf_sanity_control_file_name, job.doc.path_to_sseq_checkpoint)

    print("#**********************")
    print("Finished: Wolf Sanity GOMC control file writing")
    print("#**********************")
    #
    if (job.sp.electrostatic_method == "Wolf"):
        output_name_control_file_calibration_name = "wolf_calibration"

        if job.doc.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
            useCoul = False
            CBMC_First = (12,)
            CBMC_Nth = (10,)
            CBMC_Ang = (50,)
            CBMC_Dih = (50,)
            if job.doc.equilibration_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.4,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)

            elif job.doc.equilibration_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.39,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)
            elif job.doc.equilibration_ensemble in ["GEMC_NVT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (0.00,)
                IntraSwapFreq = (0.10,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (0.20,)
                DisFreq = (0.29,)
                RotFreq = (0.20,)
                RegrowthFreq = (0.20,)
                CrankShaftFreq = (0.0,)
            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )

        if job.doc.solute in ["ETOH", "ETOH-OPLS", "solvent_box"]:
            useCoul = True
            CBMC_First = (10,)
            CBMC_Nth = (10,)
            CBMC_Ang = (100,)
            CBMC_Dih = (50,)
            if job.doc.equilibration_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.50,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)

            elif job.doc.equilibration_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.49,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)
            elif job.doc.equilibration_ensemble in ["GEMC_NVT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (0.00,)
                IntraSwapFreq = (0.10,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (0.20,)
                DisFreq = (0.29,)
                RotFreq = (0.20,)
                RegrowthFreq = (0.20,)
                CrankShaftFreq = (0.0,)
            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )                  

        gomc_control.write_gomc_control_file(
            gomc_charmm_object_with_files,
            output_name_control_file_calibration_name,
            job.doc.equilibration_ensemble,
            Calibration_MC_steps,
            production_temperature_K,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(gomc_ff_filename_str),
            Restart=True,
            RestartCheckpoint=True,
            ExpertMode=False,
            Coordinates_box_0=job.doc.path_to_sseq_pdb,
            Structure_box_0=job.doc.path_to_sseq_psf,
            binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
            extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
            binVelocities_box_0=None,
            Coordinates_box_1=job.doc.path_to_sseq_pdb_box_1,
            Structure_box_1=job.doc.path_to_sseq_psf_box_1,
            binCoordinates_box_1=job.doc.path_to_sseq_binCoordinates_box_1,
            extendedSystem_box_1=job.doc.path_to_sseq_extendedSystem_box_1,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": production_pressure_bar,
                "Ewald": True,
                "ElectroStatic": use_ElectroStatics,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": job.doc.Rcut_ang,
                "Exclude": Exclude,
                "VolFreq": VolFreq[-1],
                "MultiParticleFreq": MultiParticleFreq[-1],
                "IntraSwapFreq": IntraSwapFreq[-1],
                "CrankShaftFreq": CrankShaftFreq[-1],
                "SwapFreq": SwapFreq[-1],
                "DisFreq": DisFreq[-1],
                "RotFreq": RotFreq[-1],
                "RegrowthFreq": RegrowthFreq[-1],
                "OutputName": output_name_control_file_calibration_name,
                "EqSteps": Calibration_MC_Eq_Steps,
                "PressureCalc": output_false_list_input,
                "RestartFreq": output_false_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_false_list_input,
                "Potential": cutoff_style,
                "LRC": True,
                "RcutLow": 1.0,
                "RcutCoulomb_box_1" : min((math.floor(job.doc.vap_box_lengths_ang/2.0)-1), 20),
                "CBMC_First": CBMC_First[-1],
                "CBMC_Nth": CBMC_Nth[-1],
                "CBMC_Ang": CBMC_Ang[-1],
                "CBMC_Dih": CBMC_Dih[-1],
            },
        )
        append_wolf_calibration_parameters(job)
        append_checkpoint_line(job, output_name_control_file_calibration_name, job.doc.path_to_sseq_checkpoint)

        ### Need to append Wolf Calibration lines since they aren't in MosDef

        print("#**********************")
        print("Completed: Wolf Calibration GOMC control file writing")
        print("#**********************")


    if (job.sp.solute == "solvent_box"):
        return
    # ******************************************************
    # equilb selected_ensemble, if NVT -> NPT - GOMC control file writing  (start)
    # Note: the control files are written for the max number of gomc_equilb_design_ensemble runs
    # so the Charmm object only needs created 1 time.
    # ******************************************************
    print("#**********************")
    print("Started: equilb NPT or GEMC-NVT GOMC control file writing")
    print("#**********************")
    job.doc.gomc_equilb_design_ensemble_dict = {}
    job.doc.gomc_production_run_ensemble_dict = {}

    for initial_state_sims_i in list(job.doc.InitialState_list):
        output_name_control_file_name = "{}_initial_state_{}".format(
            gomc_equilb_design_ensemble_control_file_name_str, initial_state_sims_i
        )

        job.doc.gomc_equilb_design_ensemble_dict.update(
            {
                initial_state_sims_i: {
                    "restart_control_file_name": restart_control_file_name_str,
                    "output_name_control_file_name": output_name_control_file_name,
                }
            }
        )

        # output all data and calc frequecy
        output_true_list_input = [
            True,
            int(gomc_output_data_every_X_steps),
        ]
        output_false_list_input = [
            False,
            int(gomc_output_data_every_X_steps),
        ]

        if job.doc.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
            useCoul = False
            CBMC_First = (12,)
            CBMC_Nth = (10,)
            CBMC_Ang = (50,)
            CBMC_Dih = (50,)
            if job.doc.equilibration_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.4,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)

            elif job.doc.equilibration_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.39,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)
            elif job.doc.equilibration_ensemble in ["GEMC_NVT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (0.00,)
                IntraSwapFreq = (0.10,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (0.20,)
                DisFreq = (0.29,)
                RotFreq = (0.20,)
                RegrowthFreq = (0.20,)
                CrankShaftFreq = (0.0,)
            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )

        if job.doc.solute in ["ETOH", "ETOH-OPLS", "solvent_box"]:
            useCoul = True
            CBMC_First = (10,)
            CBMC_Nth = (10,)
            CBMC_Ang = (100,)
            CBMC_Dih = (50,)
            if job.doc.equilibration_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.50,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)

            elif job.doc.equilibration_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.49,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)
            elif job.doc.equilibration_ensemble in ["GEMC_NVT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (0.00,)
                IntraSwapFreq = (0.10,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (0.20,)
                DisFreq = (0.29,)
                RotFreq = (0.20,)
                RegrowthFreq = (0.20,)
                CrankShaftFreq = (0.0,)
            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )                
                
        # Only use namd run from ewald to ensure both start at exact same configuration.
        gomc_control.write_gomc_control_file(
            gomc_charmm_object_with_files,
            output_name_control_file_name,
            job.doc.equilibration_ensemble,
            MC_steps,
            production_temperature_K,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(gomc_ff_filename_str),
            Restart= True,
            RestartCheckpoint=True,
            ExpertMode=False,
            #Coordinates_box_0= Coordinates_box_0 if job.sp.electrostatic_method == "Ewald" else job.doc.path_to_ref_pdb,
            #Structure_box_0=Structure_box_0 if job.sp.electrostatic_method == "Ewald" else job.doc.path_to_ref_psf,
            Coordinates_box_0=job.doc.path_to_sseq_pdb,
            Structure_box_0=job.doc.path_to_sseq_psf,
            binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
            extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
            binVelocities_box_0=None,
            Coordinates_box_1=job.doc.path_to_sseq_pdb_box_1,
            Structure_box_1=job.doc.path_to_sseq_psf_box_1,
            binCoordinates_box_1=job.doc.path_to_sseq_binCoordinates_box_1,
            extendedSystem_box_1=job.doc.path_to_sseq_extendedSystem_box_1,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": production_pressure_bar,
                "Ewald": True if job.sp.electrostatic_method == "Ewald" else False,
                "ElectroStatic": use_ElectroStatics,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": job.doc.Rcut_ang,
                "Exclude": Exclude,
                "VolFreq": VolFreq[-1],
                "MultiParticleFreq": MultiParticleFreq[-1],
                "IntraSwapFreq": IntraSwapFreq[-1],
                "CrankShaftFreq": CrankShaftFreq[-1],
                "SwapFreq": SwapFreq[-1],
                "DisFreq": DisFreq[-1],
                "RotFreq": RotFreq[-1],
                "RegrowthFreq": RegrowthFreq[-1],
                "OutputName": output_name_control_file_name,
                "EqSteps": EqSteps,
                "PressureCalc": output_false_list_input,
                "RestartFreq": restart_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_false_list_input,
                "Potential": cutoff_style,
                "LRC": True,
                "RcutLow": 0,
                "CBMC_First": CBMC_First[-1],
                "CBMC_Nth": CBMC_Nth[-1],
                "CBMC_Ang": CBMC_Ang[-1],
                "CBMC_Dih": CBMC_Dih[-1],
                "FreeEnergyCalc": FreeEnergyCalc,
                "MoleculeType": MoleculeType,
                "InitialState": initial_state_sims_i,
                "LambdaVDW": list(job.doc.LambdaVDW_list),
                "LambdaCoulomb":  list(job.doc.LambdaCoul_list) if useCoul else None,
            },
        )
        append_checkpoint_line(job, output_name_control_file_name, job.doc.path_to_sseq_checkpoint)
        print("#**********************")
        print("Completed: equilb NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")

        # ******************************************************
        # equilb selected_ensemble, if NVT -> NPT - GOMC control file writing  (end)
        # Note: the control files are written for the max number of gomc_equilb_design_ensemble runs
        # so the Charmm object only needs created 1 time.
        # ******************************************************

        # ******************************************************
        # production NPT or GEMC-NVT - GOMC control file writing  (start)
        # ******************************************************
        print("#**********************")
        print("Started: production NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")

    for initial_state_sims_i in list(job.doc.InitialState_list):
        output_name_control_file_name = "{}_initial_state_{}".format(
            gomc_production_control_file_name_str, initial_state_sims_i
        )
        restart_control_file_name_str = "{}_initial_state_{}".format(
            gomc_equilb_design_ensemble_control_file_name_str, int(initial_state_sims_i)
        )
        job.doc.gomc_production_run_ensemble_dict.update(
            {
                initial_state_sims_i: {
                    "restart_control_file_name": restart_control_file_name_str,
                    "output_name_control_file_name": output_name_control_file_name,
                }
            }
        )

        # output all data and calc frequecy
        output_true_list_input = [
            True,
            int(gomc_output_data_every_X_steps),
        ]
        output_false_list_input = [
            False,
            int(gomc_output_data_every_X_steps),
        ]

        # calc MC steps
        MC_steps = int(gomc_steps_lamda_production)
        


        # output all data and calc frequecy
        output_true_list_input = [
            True,
            int(gomc_output_data_every_X_steps),
        ]
        output_false_list_input = [
            False,
            int(gomc_output_data_every_X_steps),
        ]
        
        if job.doc.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
            useCoul = False
            CBMC_First = (12,)
            CBMC_Nth = (10,)
            CBMC_Ang = (50,)
            CBMC_Dih = (50,)
            if job.doc.production_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.4,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)

            elif job.doc.production_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (None,)
                SwapFreq = (None,)
                DisFreq = (0.39,)
                RotFreq = (0.3,)
                RegrowthFreq = (0.3,)

            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )

        if job.doc.solute in ["ETOH", "ETOH-OPLS", "solvent_box"]:
            useCoul = True
            CBMC_First = (10,)
            CBMC_Nth = (10,)
            CBMC_Ang = (100,)
            CBMC_Dih = (50,)
            if job.doc.production_ensemble in ["NVT"]:
                VolFreq = (0.00,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.50,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)

            elif job.doc.production_ensemble in ["NPT"]:
                VolFreq = (0.01,)
                MultiParticleFreq = (None,)
                IntraSwapFreq = (0.0,)
                CrankShaftFreq = (0.1,)
                SwapFreq = (None,)
                DisFreq = (0.49,)
                RotFreq = (0.2,)
                RegrowthFreq = (0.20,)

            else:
                raise ValueError(
                    "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                    "in the GOMC control file writer."
                )

        Prod_Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
            restart_control_file_name_str
        )
        Prod_Structure_box_0 = "{}_BOX_0_restart.psf".format(
            restart_control_file_name_str
        )
        Prod_binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
            restart_control_file_name_str
        )
        Prod_extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
            restart_control_file_name_str
        )

        gomc_control.write_gomc_control_file(
            gomc_charmm_object_with_files,
            output_name_control_file_name,
            job.doc.production_ensemble,
            MC_steps,
            production_temperature_K,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(gomc_ff_filename_str),
            Restart=True,
            RestartCheckpoint=True,
            ExpertMode=False,
            Coordinates_box_0=Prod_Coordinates_box_0,
            Structure_box_0=Prod_Structure_box_0,
            binCoordinates_box_0=Prod_binCoordinates_box_0,
            extendedSystem_box_0=Prod_extendedSystem_box_0,
            binVelocities_box_0=None,
            Coordinates_box_1=None,
            Structure_box_1=None,
            binCoordinates_box_1=None,
            extendedSystem_box_1=None,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": production_pressure_bar,
                "Ewald": True if job.sp.electrostatic_method == "Ewald" else False,
                "ElectroStatic": use_ElectroStatics,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": job.doc.Rcut_ang,
                "Exclude": Exclude,
                "VolFreq": VolFreq[-1],
                "MultiParticleFreq": MultiParticleFreq[-1],
                "IntraSwapFreq": IntraSwapFreq[-1],
                "CrankShaftFreq": CrankShaftFreq[-1],
                "SwapFreq": SwapFreq[-1],
                "DisFreq": DisFreq[-1],
                "RotFreq": RotFreq[-1],
                "RegrowthFreq": RegrowthFreq[-1],
                "OutputName": output_name_control_file_name,
                "EqSteps": EqSteps,
                "PressureCalc": output_false_list_input,
                "RestartFreq": restart_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_false_list_input,
                "Potential": cutoff_style,
                "LRC": True,
                "RcutLow": 0,
                "CBMC_First": CBMC_First[-1],
                "CBMC_Nth": CBMC_Nth[-1],
                "CBMC_Ang": CBMC_Ang[-1],
                "CBMC_Dih": CBMC_Dih[-1],
                "FreeEnergyCalc": FreeEnergyCalc,
                "MoleculeType": MoleculeType,
                "InitialState": initial_state_sims_i,
                "LambdaVDW": list(job.doc.LambdaVDW_list),
                "LambdaCoulomb":  list(job.doc.LambdaCoul_list) if useCoul else None,
            },
        )
        append_checkpoint_line(job, output_name_control_file_name, job.fn("{}_restart.chk".format(restart_control_file_name_str)))

        print("#**********************")
        print("Completed: production NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")
        # ******************************************************
        # production NPT or GEMC-NVT - GOMC control file writing  (end)
        # ******************************************************

# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# namd_equilb_NPT -starting the NAMD simulations (start)
# ******************************************************
# ******************************************************
# Only run namd on the Ewald directories, then use the same 
# final trajectory for Wolf.
@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(mosdef_input_written)
@Project.pre(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.pre(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.post(part_3a_output_namd_equilb_NVT_box_0_started)
@Project.post(part_3a_output_namd_equilb_NVT_box_1_started)
@Project.post(part_4a_job_namd_equilb_NVT_box_0_completed_properly)
@Project.post(part_4a_job_namd_equilb_NVT_box_1_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.namd_node_ncpu,
        "ngpu": lambda job: job.doc.namd_node_ngpu,
        "memory": memory_needed,
        "walltime": walltime_namd_hr,
    }
)
@flow.with_job
@flow.cmd
def run_namd_equilb_NVT_box_0_gomc_command(job):
    """Run the namd_equilb_NPT simulation."""
    print("#**********************")
    print("# Started the run_namd_equilb_NVT_box_0_gomc_command.")
    print("#**********************")
    """Run the gomc_calibration_run_ensemble simulation."""
    
    control_file_name_str = namd_equilb_NVT_control_file_box_0_name_str
    
    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(namd_binary_path),
        str(job.doc.namd_equilb_NVT_gomc_binary_file),
        str(job.doc.namd_node_ncpu),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    print('namd run_command = ' + str(run_command))

    return run_command
    """
    run_command = "echo namdcopied"
    print('gomc gomc_sseq_run_ensemble run_command = ' + str(run_command))
    
    return run_command
    """

# ******************************************************
# ******************************************************
# namd_equilb_NPT -starting the NAMD simulations (start)
# ******************************************************
# ******************************************************
# Only run namd on the Ewald directories, then use the same 
# final trajectory for Wolf.
@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(mosdef_input_written)
@Project.pre(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.pre(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.post(part_3a_output_namd_equilb_NVT_box_0_started)
@Project.post(part_3a_output_namd_equilb_NVT_box_1_started)
@Project.post(part_4a_job_namd_equilb_NVT_box_0_completed_properly)
@Project.post(part_4a_job_namd_equilb_NVT_box_1_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.namd_node_ncpu,
        "ngpu": lambda job: job.doc.namd_node_ngpu,
        "memory": memory_needed,
        "walltime": walltime_namd_hr,
    }
)
@flow.with_job
@flow.cmd
def run_namd_equilb_NVT_box_1_gomc_command(job):
    """Run the namd_equilb_NPT simulation."""
    print("#**********************")
    print("# Started the run_namd_equilb_NVT_box_1_gomc_command.")
    print("#**********************")
    """Run the gomc_calibration_run_ensemble simulation."""
    
    control_file_name_str = namd_equilb_NVT_control_file_box_1_name_str
    
    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(namd_binary_path),
        str(job.doc.namd_equilb_NVT_gomc_binary_file),
        str(job.doc.namd_node_ncpu),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    print('namd run_command = ' + str(run_command))

    return run_command
    """
    run_command = "echo namdcopied"
    print('gomc gomc_sseq_run_ensemble run_command = ' + str(run_command))
    
    return run_command
    """
# ******************************************************
# ******************************************************
# namd_equilb_NPT -starting the NAMD simulations (end)
# ******************************************************
# ******************************************************
# ******************************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(part_4a_job_namd_equilb_NVT_box_0_completed_properly)
@Project.pre(part_4a_job_namd_equilb_NVT_box_1_completed_properly)
@Project.pre(mosdef_input_written)
@Project.pre(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.pre(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.post(part_3b_output_gomc_sseq_started)
@Project.post(part_4b_job_gomc_sseq_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.gomc_ncpu,
        "ngpu": lambda job: job.doc.gomc_ngpu,
        "memory": memory_needed,
        "walltime": walltime_gomc_equilbrium_hr,
    }
)
@flow.with_job
@flow.cmd
def run_sseq_run_gomc_command(job):

    Single_state_gomc_eq_control_file_name = "single_state_eq"

    print(f"Running simulation job id {job}")

    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(gomc_binary_path),
        str(job.doc.gomc_equilb_design_ensemble_gomc_binary_file),
        str(job.doc.gomc_ncpu),
        str(Single_state_gomc_eq_control_file_name),
        str(Single_state_gomc_eq_control_file_name),
    )
    print('gomc gomc_sseq_run_ensemble run_command = ' + str(run_command))
    return run_command

#@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_model != "Calibrator")
@Project.pre(part_1a_initial_data_input_to_json)
@Project.pre(mosdef_input_written)
@Project.pre(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.pre(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.pre(part_4b_job_gomc_wolf_parameters_found)
@Project.pre(part_4b_job_gomc_wolf_parameters_appended)
@Project.post(part_3b_output_gomc_wolf_sanity_started)
@Project.post(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.gomc_ncpu if job.sp.electrostatic_method == "Ewald" else 8,
        "ngpu": lambda job: job.doc.gomc_ngpu if (job.sp.electrostatic_method == "Ewald") else 0,
        "memory": memory_needed,
        "walltime": walltime_gomc_equilbrium_hr,
    }
)
@flow.with_job
@flow.cmd
def run_wolf_sanity_run_gomc_command(job):
    """Run the gomc_calibration_run_ensemble simulation."""
    wolf_sanity_control_file_name = "wolf_sanity"

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(gomc_binary_path),
        "GOMC_GPU_GEMC" if (job.sp.electrostatic_method == "Ewald") else "GOMC_CPU_GEMC",
        str(job.doc.gomc_ncpu) if (job.sp.electrostatic_method == "Ewald") else str(8),
        str(wolf_sanity_control_file_name),
        str(wolf_sanity_control_file_name),
    )

    print('gomc gomc_wolf_sanity_run run_command = ' + str(run_command))

    return run_command

# ******************************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Calibrator")
@Project.pre(lambda j: j.sp.wolf_model == "Calibrator")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(mosdef_input_written)
@Project.pre(part_2a_namd_equilb_NVT_box_0_control_file_written)
@Project.pre(part_2a_namd_equilb_NVT_box_1_control_file_written)
@Project.pre(part_2b_gomc_equilb_design_ensemble_control_file_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.pre(part_4a_job_namd_equilb_NVT_box_0_completed_properly)
@Project.pre(part_4a_job_namd_equilb_NVT_box_1_completed_properly)
@Project.post(part_3b_output_gomc_calibration_started)
@Project.post(part_4b_job_gomc_calibration_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.gomc_ncpu,
        "ngpu": lambda job: job.doc.gomc_ngpu,
        "memory": memory_needed,
        "walltime": 26,
    }
)
@flow.with_job
@flow.cmd
def run_calibration_run_gomc_command(job):
    """Run the gomc_calibration_run_ensemble simulation."""
    control_file_name_str = "wolf_calibration"

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(gomc_binary_path),
        str(job.doc.gomc_calibration_gomc_binary_file),
        str(job.doc.gomc_ncpu),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    print('gomc gomc_calibration_run_ensemble run_command = ' + str(run_command))
    return run_command


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Calibrator")
@Project.pre(lambda j: j.sp.wolf_model == "Calibrator")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_4b_job_gomc_calibration_completed_properly)
@Project.post(part_4b_job_gomc_wolf_parameters_found)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def part_4b_job_gomc_calibration_find_minimum(job):

    from src.utils.surface import find_minimum
    import pickle as pickle
    import re
    regex = re.compile("Wolf_Calibration_(\w+?)_(\w+?)_BOX_(\d+)_(\w+?).dat")

    bestValueFileName = "bestWolfParameters"
    if (not job.isfile(bestValueFileName+".pickle")):
        model2BestWolfAlphaRCut = dict()
        for root, dirs, files in os.walk(job.fn("")):
            for file in files:
                if regex.match(file):
                    groups = regex.search(file)
                    wolfKind = groups.group(1)
                    potential = groups.group(2)
                    box = groups.group(3)
                    tupleMin = find_minimum(job.fn(file), job.sp.solute, wolfKind, potential, box, True)
                    # Use smaller error, either BF or Grad Desc
                    model2BestWolfAlphaRCut[(wolfKind, potential, box)] = dict(tupleMin)
        with open(bestValueFileName+".pickle", 'wb') as handle:
            pickle.dump(model2BestWolfAlphaRCut, handle, protocol=pickle.HIGHEST_PROTOCOL)

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Calibrator")
@Project.pre(lambda j: j.sp.wolf_model == "Calibrator")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_4b_job_gomc_calibration_completed_properly)
@Project.post(part_4b_job_gomc_all_surface_plot_created)

@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def part_4b_job_gomc_plot_surfaces(job):

    from src.utils.surface import plot_all_surfaces
    import pickle as pickle
    import re
    regex = re.compile("Wolf_Calibration_(\w+?)_(\w+?)_BOX_(\d+)_(\w+?).dat")
    try:
        model2BestWolfAlphaRCut = dict()
        for root, dirs, files in os.walk(job.fn("")):
            for file in files:
                if regex.match(file):
                    groups = regex.search(file)
                    wolfKind = groups.group(1)
                    potential = groups.group(2)
                    box = groups.group(3)
                    allSurfacesFile = job.sp.solute+"_"+wolfKind+"_"+potential+"_Box_"+box+"_allSurfaces.html"
                    if (not job.isfile(allSurfacesFile+".html")):
                        plot_all_surfaces(pr, job, file, job.sp.solute, wolfKind, potential, box, True)
    except:
        return False
# ******************************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************


@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Calibrator")
@Project.pre(lambda j: j.sp.wolf_model == "Calibrator")
@Project.pre(lambda j: j.sp.solute == "solvent_box")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_4b_wolf_sanity_analysis_completed)
@Project.post(part_4b_wolf_sanity_histograms_created)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def part_4b_create_wolf_sanity_histograms(job):
    df1 = pd.DataFrame()
    ewald_sp = job.statepoint()
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_model']="Calibrator"        
    ewald_sp['wolf_potential']="Calibrator"   
    ewald_sp['solute']="solvent_box"   
    ewald_sp['replica_number_int']=0
    jobs = list(pr.find_jobs(ewald_sp))
    try:
        for ewald_job in jobs:
            if (ewald_job.isfile("wolf_sanity_all_energies.csv")):
                df1 = pd.read_csv (ewald_job.fn('wolf_sanity_all_energies.csv'), sep=',', header=0, na_values='NaN', index_col=0)
            else:
                return False
    except:
        return False

    print(df1)

    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.stats as st

    numBins = 100
    nskip = 100
    ref_ewald = df1["Ewald_Ewald"]

    from pymbar import timeseries
    t0, g, Neff_max = timeseries.detectEquilibration(ref_ewald, nskip=nskip) # compute indices of uncorrelated timeseries
    A_t_equil_ewald = ref_ewald[t0:]
    A_t_equil_steps_ewald = ref_ewald[t0:]

    Col_Dict = {"GROSS_DSF": "Waibel2018a", "VLUGT_DSF": 'Rahbari', "VLUGTWINTRACUTOFF_DSF": 'Waibel2018b',
    "GROSS_DSP": 'Waibel2018a', "VLUGT_DSP": 'Rahbari', "VLUGTWINTRACUTOFF_DSP": 'Waibel2018b'}

    colList = df1.columns.tolist()
    colList.remove("Ewald_Ewald")
    colList.remove("steps")

    figSP, axs = plt.subplots(2, 3)
    counter = 0
    for col, col_i in zip(colList, range(0, len(colList))):

        wolf = df1[col]
        t0, g, Neff_max = timeseries.detectEquilibration(wolf, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil_wolf = wolf[t0:]
        A_t_equil_steps_wolf = wolf[t0:]


        ref_min = min(A_t_equil_ewald)
        ref_max = max(A_t_equil_ewald)

        wolf_min = min(A_t_equil_wolf)
        wolf_max = max(A_t_equil_wolf)

        print(ref_min)
        print(ref_max)

        print(wolf_min)
        print(wolf_max)

        xmin = min(ref_min, wolf_min)
        xmax = min(ref_max, wolf_max)

        binWidth =  (xmax - xmin)/float(numBins)
        binList = np.arange(xmin, xmax+binWidth, binWidth)
        # estimate the line with probability density function (PDF)
        kde1 = st.gaussian_kde(A_t_equil_ewald).pdf(binList)

        #Plot Ewald
        plt.plot(binList, kde1, color="black", linewidth=2, label="Ewald")

        kde2 = st.gaussian_kde(A_t_equil_wolf).pdf(binList)
        #plt.hist(wolf, density=True, bins=binList, alpha=1, label=col)  # density=False would make counts
        plt.plot(binList, kde2, color="red", linewidth=2, label=Col_Dict[col])
        plt.xlim(xmin, xmax)
        plt.ylabel('Probability')
        plt.xlabel('Total Energy (K)')
        #plt.legend()
        plt.rcParams.update({'font.size': 22})
        axs[counter / 3, counter % 3].plot(binList, kde1, color="black", linewidth=2, label="Ewald")
        axs[counter / 3, counter % 3].plot(binList, kde2, color="red", linewidth=2, label=Col_Dict[col])

        plt.savefig("PotentialEnergyDistribution_Ewald_vs_{}".format(col), dpi=300, bbox_inches='tight')
        plt.figure().clear()
    figSP.savefig("PotentialEnergyDistribution_Ewald_vs_All", dpi=300, bbox_inches='tight')

# ******************************************************
# ******************************************************
# data analysis - get the average data from each replicate (start)
# ******************************************************
# ******************************************************


@Project.operation.with_directives(
     {
         "np": 1,
         "ngpu": 0,
         "memory": memory_needed,
         "walltime": walltime_gomc_analysis_hr,
     }
)
@FlowProject.pre(
     lambda * jobs: all(
         part_4b_job_gomc_wolf_sanity_completed_properly(job, gomc_production_control_file_name_str)
         for job in jobs
     )
)
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.post(part_5a_analysis_individual_simulation_averages_completed)
@flow.with_job
def part_5a_analysis_individual_simulation_averages(*jobs):
    # remove the total averged replicate data and all analysis data after this,
    # as it is no longer valid when adding more simulations
    if os.path.isfile(f'../../analysis/{output_avg_std_of_replicates_txt_file_name_liq}'):
        os.remove(f'../../analysis/{output_avg_std_of_replicates_txt_file_name_liq}')
    if os.path.isfile(f'../../analysis/{output_avg_std_of_replicates_txt_file_name_vap}'):
        os.remove(f'../../analysis/{output_avg_std_of_replicates_txt_file_name_vap}')
    if os.path.isfile(f'../../analysis/{output_critical_data_replicate_txt_file_name}'):
        os.remove(f'../../analysis/{output_critical_data_replicate_txt_file_name}')
    if os.path.isfile(f'../../analysis/{output_critical_data_avg_std_of_replicates_txt_file_name}'):
        os.remove(f'../../analysis/{output_critical_data_avg_std_of_replicates_txt_file_name}')
    if os.path.isfile(f'../../analysis/{output_boiling_data_replicate_txt_file_name}'):
        os.remove(f'../../analysis/{output_boiling_data_replicate_txt_file_name}')
    if os.path.isfile(f'../../analysis/{output_boiling_data_avg_std_of_replicates_txt_file_name}'):
        os.remove(f'../../analysis/{output_boiling_data_avg_std_of_replicates_txt_file_name}')


    # this is set to basically use all values.  However, allows the ability to set if needed
    step_start = 0 * 10 ** 6
    step_finish = 1 * 10 ** 12

    # get the averages from each individual simulation and write the csv's.
    for job in jobs:

        reading_file_box_0 = job.fn(f'Blk_{gomc_production_control_file_name_str}_BOX_0.dat')
        reading_file_box_1 = job.fn(f'Blk_{gomc_production_control_file_name_str}_BOX_1.dat')

        output_column_temp_title = 'temp_K'  # column title title for temp
        output_column_no_step_title = 'Step'  # column title title for iter value
        output_column_no_pressure_title = 'P_bar'  # column title title for PRESSURE
        output_column_total_molecules_title = "No_mol"  # column title title for TOT_MOL
        output_column_Rho_title = 'Rho_kg_per_m_cubed'  # column title title for TOT_DENS
        output_column_box_volume_title = 'V_ang_cubed'  # column title title for VOLUME
        output_column_box_length_if_cubed_title = 'L_m_if_cubed'  # column title title for VOLUME
        output_column_box_Hv_title = 'Hv_kJ_per_mol'  # column title title for HEAT_VAP

        blk_file_reading_column_no_step_title = '#STEP'  # column title title for iter value
        blk_file_reading_column_no_pressure_title = 'PRESSURE'  # column title title for PRESSURE
        blk_file_reading_column_total_molecules_title = "TOT_MOL"  # column title title for TOT_MOL
        blk_file_reading_column_Rho_title = 'TOT_DENS'  # column title title for TOT_DENS
        blk_file_reading_column_box_volume_title = 'VOLUME'  # column title title for VOLUME
        blk_file_reading_column_box_Hv_title = 'HEAT_VAP'  # column title title for HEAT_VAP

        # Programmed data
        step_start_string = str(int(step_start))
        step_finish_string = str(int(step_finish))

        # *************************
        # drawing in data from single file and extracting specific rows for the liquid box (start)
        # *************************
        data_box_0 = pd.read_csv(reading_file_box_0, sep='\s+', header=0, na_values='NaN', index_col=False)

        data_box_0 = pd.DataFrame(data_box_0)
        step_no_title_mod = blk_file_reading_column_no_step_title[1:]
        header_list = list(data_box_0.columns)
        header_list[0] = step_no_title_mod
        data_box_0.columns = header_list

        data_box_0 = data_box_0.query(step_start_string + ' <= ' + step_no_title_mod + ' <= ' + step_finish_string)

        iter_no_box_0 = data_box_0.loc[:, step_no_title_mod]
        iter_no_box_0 = list(iter_no_box_0)
        iter_no_box_0 = np.transpose(iter_no_box_0)

        pressure_box_0 = data_box_0.loc[:, blk_file_reading_column_no_pressure_title]
        pressure_box_0 = list(pressure_box_0)
        pressure_box_0 = np.transpose(pressure_box_0)
        pressure_box_0_mean = np.nanmean(pressure_box_0)

        total_molecules_box_0 = data_box_0.loc[:, blk_file_reading_column_total_molecules_title]
        total_molecules_box_0 = list(total_molecules_box_0)
        total_molecules_box_0 = np.transpose(total_molecules_box_0)
        total_molecules_box_0_mean = np.nanmean(total_molecules_box_0)

        Rho_box_0 = data_box_0.loc[:, blk_file_reading_column_Rho_title]
        Rho_box_0 = list(Rho_box_0)
        Rho_box_0 = np.transpose(Rho_box_0)
        Rho_box_0_mean = np.nanmean(Rho_box_0)

        volume_box_0 = data_box_0.loc[:, blk_file_reading_column_box_volume_title]
        volume_box_0 = list(volume_box_0)
        volume_box_0 = np.transpose(volume_box_0)
        volume_box_0_mean = np.nanmean(volume_box_0)
        length_if_cube_box_0_mean = (volume_box_0_mean) ** (1 / 3)

        Hv_box_0 = data_box_0.loc[:, blk_file_reading_column_box_Hv_title]
        Hv_box_0 = list(Hv_box_0)
        Hv_box_0 = np.transpose(Hv_box_0)
        Hv_box_0_mean = np.nanmean(Hv_box_0)

        # *************************
        # drawing in data from single file and extracting specific rows for the liquid box (end)
        # *************************

        # *************************
        # drawing in data from single file and extracting specific rows for the vapor box (start)
        # *************************
        data_box_1 = pd.read_csv(reading_file_box_1, sep='\s+', header=0, na_values='NaN', index_col=False)
        data_box_1 = pd.DataFrame(data_box_1)

        step_no_title_mod = blk_file_reading_column_no_step_title[1:]
        header_list = list(data_box_1.columns)
        header_list[0] = step_no_title_mod
        data_box_1.columns = header_list

        data_box_1 = data_box_1.query(step_start_string + ' <= ' + step_no_title_mod + ' <= ' + step_finish_string)

        iter_no_box_1 = data_box_1.loc[:, step_no_title_mod]
        iter_no_box_1 = list(iter_no_box_1)
        iter_no_box_1 = np.nanmean(iter_no_box_1)

        pressure_box_1 = data_box_1.loc[:, blk_file_reading_column_no_pressure_title]
        pressure_box_1 = list(pressure_box_1)
        pressure_box_1 = np.transpose(pressure_box_1)
        pressure_box_1_mean = np.nanmean(pressure_box_1)

        total_molecules_box_1 = data_box_1.loc[:, blk_file_reading_column_total_molecules_title]
        total_molecules_box_1 = list(total_molecules_box_1)
        total_molecules_box_1 = np.transpose(total_molecules_box_1)
        total_molecules_box_1_mean = np.nanmean(total_molecules_box_1)

        Rho_box_1 = data_box_1.loc[:, blk_file_reading_column_Rho_title]
        Rho_box_1 = list(Rho_box_1)
        Rho_box_1 = np.transpose(Rho_box_1)
        Rho_box_1_mean = np.nanmean(Rho_box_1)

        volume_box_1 = data_box_1.loc[:, blk_file_reading_column_box_volume_title]
        volume_box_1 = list(volume_box_1)
        volume_box_1 = np.transpose(volume_box_1)
        volume_box_1_mean = np.nanmean(volume_box_1)
        length_if_cube_box_1_mean = (volume_box_1_mean) ** (1 / 3)

        Hv_box_1 = data_box_1.loc[:, blk_file_reading_column_box_Hv_title]
        Hv_box_1 = list(Hv_box_1)
        Hv_box_1 = np.transpose(Hv_box_1)
        Hv_box_1_mean = np.nanmean(Hv_box_1)

        # sort boxes based on density to liquid or vapor
        if (Rho_box_0_mean > Rho_box_1_mean) or (Rho_box_0_mean == Rho_box_1_mean):
            pressure_box_liq_mean = pressure_box_0_mean
            total_molecules_box_liq_mean = total_molecules_box_0_mean
            Rho_box_liq_mean = Rho_box_0_mean
            volume_box_liq_mean = volume_box_0_mean
            length_if_cube_box_liq_mean = length_if_cube_box_0_mean
            Hv_box_liq_mean = Hv_box_0_mean

            pressure_box_vap_mean = pressure_box_1_mean
            total_molecules_box_vap_mean = total_molecules_box_1_mean
            Rho_box_vap_mean = Rho_box_1_mean
            volume_box_vap_mean = volume_box_1_mean
            length_if_cube_box_vap_mean = length_if_cube_box_1_mean
            Hv_box_vap_mean = Hv_box_1_mean

        elif Rho_box_0_mean < Rho_box_1_mean:
            pressure_box_liq_mean = pressure_box_1_mean
            total_molecules_box_liq_mean = total_molecules_box_1_mean
            Rho_box_liq_mean = Rho_box_1_mean
            volume_box_liq_mean = volume_box_1_mean
            length_if_cube_box_liq_mean = length_if_cube_box_1_mean
            Hv_box_liq_mean = Hv_box_1_mean

            pressure_box_vap_mean = pressure_box_0_mean
            total_molecules_box_vap_mean = total_molecules_box_0_mean
            Rho_box_vap_mean = Rho_box_0_mean
            volume_box_vap_mean = volume_box_0_mean
            length_if_cube_box_vap_mean = length_if_cube_box_0_mean
            Hv_box_vap_mean = Hv_box_0_mean

        # *************************
        # drawing in data from single file and extracting specific rows for the vapor box (end)
        # *************************

        box_liq_replicate_data_txt_file = open(output_replicate_txt_file_name_liq, "w")
        box_liq_replicate_data_txt_file.write(
            f"{output_column_temp_title: <30} "
            f"{output_column_no_pressure_title: <30} "
            f"{output_column_total_molecules_title: <30} "
            f"{output_column_Rho_title: <30} "
            f"{output_column_box_volume_title: <30} "
            f"{output_column_box_length_if_cubed_title: <30} "
            f"{output_column_box_Hv_title: <30} "
            f" \n"
        )
        box_liq_replicate_data_txt_file.write(
            f"{job.sp.production_temperature_K: <30} "
            f"{pressure_box_liq_mean: <30} "
            f"{total_molecules_box_liq_mean: <30} "
            f"{Rho_box_liq_mean: <30} "
            f"{volume_box_liq_mean: <30} "
            f"{length_if_cube_box_liq_mean: <30} "
            f"{Hv_box_liq_mean: <30} "
            f" \n"
        )

        box_vap_replicate_data_txt_file = open(output_replicate_txt_file_name_vap, "w")
        box_vap_replicate_data_txt_file.write(
            f"{output_column_temp_title: <30} "
            f"{output_column_no_pressure_title: <30} "
            f"{output_column_total_molecules_title: <30} "
            f"{output_column_Rho_title: <30} "
            f"{output_column_box_volume_title: <30} "
            f"{output_column_box_length_if_cubed_title: <30} "
            f"{output_column_box_Hv_title: <30} "
            f" \n"
        )
        box_vap_replicate_data_txt_file.write(
            f"{job.sp.production_temperature_K: <30} "
            f"{pressure_box_vap_mean: <30} "
            f"{total_molecules_box_vap_mean: <30} "
            f"{Rho_box_vap_mean: <30} "
            f"{volume_box_vap_mean: <30} "
            f"{length_if_cube_box_vap_mean: <30} "
            f"{Hv_box_vap_mean: <30} "
            f" \n"
        )

        box_liq_replicate_data_txt_file.close()
        box_vap_replicate_data_txt_file.close()


        # ***********************
        # calc the avg data from the liq and vap boxes (end)
        # ***********************

# ******************************************************
# ******************************************************
# data analysis - get the average data from each replicate (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# data analysis - get the average and std. dev. from/across all the replicate (start)
# ******************************************************
# ******************************************************
#when you add replicates, uncomment this out
"""
@aggregator.groupby(key=statepoint_without_replica, sort_by="production_temperature_K", sort_ascending=False)
@Project.operation.with_directives(
     {
         "np": 1,
         "ngpu": 0,
         "memory": memory_needed,
         "walltime": walltime_gomc_analysis_hr,
     }
)
@Project.pre(lambda *jobs: all(part_5a_analysis_individual_simulation_averages_completed(j)
                               for j in jobs[0]._project))
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.post(part_5b_analysis_replica_averages_completed)
def part_5b_analysis_replica_averages(*jobs):
    # ***************************************************
    #  create the required lists and file labels total averages across the replicates (start)
    # ***************************************************
    # get the list used in this function
    temp_repilcate_list = []

    pressure_repilcate_box_liq_list = []
    total_molecules_repilcate_box_liq_list = []
    Rho_repilcate_box_liq_list = []
    volume_repilcate_box_liq_list = []
    length_if_cube_repilcate_box_liq_list = []
    Hv_repilcate_box_liq_list = []

    pressure_repilcate_box_vap_list = []
    total_molecules_repilcate_box_vap_list = []
    Rho_repilcate_box_vap_list = []
    volume_repilcate_box_vap_list = []
    length_if_cube_repilcate_box_vap_list = []
    Hv_repilcate_box_vap_list = []

    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_no_pressure_title = 'P_bar'  # column title title for PRESSURE
    output_column_total_molecules_title = "No_mol"  # column title title for TOT_MOL
    output_column_Rho_title = 'Rho_kg_per_m_cubed'  # column title title for TOT_DENS
    output_column_box_volume_title = 'V_ang_cubed'  # column title title for VOLUME
    output_column_box_length_if_cubed_title = 'L_m_if_cubed'  # column title title for VOLUME
    output_column_box_Hv_title = 'Hv_kJ_per_mol'  # column title title for HEAT_VAP

    output_column_temp_std_title = 'temp_std_K'  # column title title for temp
    output_column_no_pressure_std_title = 'P_std_bar'  # column title title for PRESSURE
    output_column_total_molecules_std_title = "No_mol_std"  # column title title for TOT_MOL
    output_column_Rho_std_title = 'Rho_std_kg_per_m_cubed'  # column title title for TOT_DENS
    output_column_box_volume_std_title = 'V_std_ang_cubed'  # column title title for VOLUME
    output_column_box_length_if_cubed_std_title = 'L_std_m_if_cubed'  # column title title for VOLUME
    output_column_box_Hv_std_title = 'Hv_std_kJ_per_mol'  # column title title for HEAT_VAP

    output_txt_file_header = f"{output_column_temp_title: <30} "\
                             f"{output_column_temp_std_title: <30} "\
                             f"{output_column_no_pressure_title: <30} "\
                             f"{output_column_no_pressure_std_title: <30} "\
                             f"{output_column_total_molecules_title: <30} "\
                             f"{output_column_total_molecules_std_title: <30} "\
                             f"{output_column_Rho_title: <30} "\
                             f"{output_column_Rho_std_title: <30} "\
                             f"{output_column_box_volume_title: <30} "\
                             f"{output_column_box_volume_std_title: <30} "\
                             f"{output_column_box_length_if_cubed_title: <30} "\
                             f"{output_column_box_length_if_cubed_std_title: <30} "\
                             f"{output_column_box_Hv_title: <30} "\
                             f"{output_column_box_Hv_std_title: <30} "\
                             f"\n"


    write_file_path_and_name_liq = f'analysis/{output_avg_std_of_replicates_txt_file_name_liq}'
    write_file_path_and_name_vap = f'analysis/{output_avg_std_of_replicates_txt_file_name_vap}'
    if os.path.isfile(write_file_path_and_name_liq):
        box_liq_data_txt_file = open(write_file_path_and_name_liq, "a")
    else:
        box_liq_data_txt_file = open(write_file_path_and_name_liq, "w")
        box_liq_data_txt_file.write(output_txt_file_header)

    if os.path.isfile(write_file_path_and_name_vap):
        box_vap_data_txt_file = open(write_file_path_and_name_vap, "a")
    else:
        box_vap_data_txt_file = open(write_file_path_and_name_vap, "w")
        box_vap_data_txt_file.write(output_txt_file_header)

    # ***************************************************
    # create the required lists and file labels total averages across the replicates (end)
    # ***************************************************

    for job in jobs:

        # *************************
        # drawing in data from single simulation file and extracting specific
        # *************************
        reading_file_box_liq = job.fn(output_replicate_txt_file_name_liq)
        reading_file_box_vap = job.fn(output_replicate_txt_file_name_vap)

        data_box_liq = pd.read_csv(reading_file_box_liq, sep='\s+', header=0, na_values='NaN', index_col=False)
        data_box_liq = pd.DataFrame(data_box_liq)

        pressure_repilcate_box_liq = data_box_liq.loc[:, output_column_no_pressure_title][0]
        total_molecules_repilcate_box_liq = data_box_liq.loc[:, output_column_total_molecules_title][0]
        Rho_repilcate_box_liq = data_box_liq.loc[:, output_column_Rho_title][0]
        volume_repilcate_box_liq = data_box_liq.loc[:, output_column_box_volume_title][0]
        length_if_cube_repilcate_box_liq = (volume_repilcate_box_liq) ** (1 / 3)
        Hv_repilcate_box_liq = data_box_liq.loc[:, output_column_box_Hv_title][0]

        # *************************
        # drawing in data from single file and extracting specific rows for the liquid box (end)
        # *************************

        # *************************
        # drawing in data from single file and extracting specific rows for the vapor box (start)
        # *************************
        data_box_vap = pd.read_csv(reading_file_box_vap, sep='\s+', header=0, na_values='NaN', index_col=False)
        data_box_vap = pd.DataFrame(data_box_vap)

        pressure_repilcate_box_vap = data_box_vap.loc[:, output_column_no_pressure_title][0]
        total_molecules_repilcate_box_vap = data_box_vap.loc[:, output_column_total_molecules_title][0]
        Rho_repilcate_box_vap = data_box_vap.loc[:, output_column_Rho_title][0]
        volume_repilcate_box_vap = data_box_vap.loc[:, output_column_box_volume_title][0]
        length_if_cube_repilcate_box_vap = (volume_repilcate_box_vap) ** (1 / 3)
        Hv_repilcate_box_vap = data_box_vap.loc[:, output_column_box_Hv_title][0]


        temp_repilcate_list.append(job.sp.production_temperature_K)

        pressure_repilcate_box_liq_list.append(pressure_repilcate_box_liq)
        total_molecules_repilcate_box_liq_list.append(total_molecules_repilcate_box_liq)
        Rho_repilcate_box_liq_list.append(Rho_repilcate_box_liq)
        volume_repilcate_box_liq_list.append(volume_repilcate_box_liq)
        length_if_cube_repilcate_box_liq_list.append(length_if_cube_repilcate_box_liq)
        Hv_repilcate_box_liq_list.append(Hv_repilcate_box_liq)

        pressure_repilcate_box_vap_list.append(pressure_repilcate_box_vap)
        total_molecules_repilcate_box_vap_list.append(total_molecules_repilcate_box_vap)
        Rho_repilcate_box_vap_list.append(Rho_repilcate_box_vap)
        volume_repilcate_box_vap_list.append(volume_repilcate_box_vap)
        length_if_cube_repilcate_box_vap_list.append(length_if_cube_repilcate_box_vap)
        Hv_repilcate_box_vap_list.append(Hv_repilcate_box_vap)

        # *************************
        # drawing in data from single file and extracting specific rows for the vapor box (end)
        # *************************

    # *************************
    # get the replica means and std.devs (start)
    # *************************
    temp_mean = np.mean(temp_repilcate_list)
    temp_std = np.std(temp_repilcate_list, ddof=1)

    pressure_mean_box_liq = np.mean(pressure_repilcate_box_liq_list)
    total_molecules_mean_box_liq = np.mean(total_molecules_repilcate_box_liq_list)
    Rho_mean_box_liq = np.mean(Rho_repilcate_box_liq_list)
    volume_mean_box_liq = np.mean(volume_repilcate_box_liq_list)
    length_if_cube_mean_box_liq = np.mean(length_if_cube_repilcate_box_liq_list)
    Hv_mean_box_liq = np.mean(Hv_repilcate_box_liq_list)

    pressure_std_box_liq = np.std(pressure_repilcate_box_liq_list, ddof=1)
    total_molecules_std_box_liq = np.std(total_molecules_repilcate_box_liq_list, ddof=1)
    Rho_std_box_liq= np.std(Rho_repilcate_box_liq_list, ddof=1)
    volume_std_box_liq = np.std(volume_repilcate_box_liq_list, ddof=1)
    length_if_cube_std_box_liq = np.std(length_if_cube_repilcate_box_liq_list, ddof=1)
    Hv_std_box_liq = np.std(Hv_repilcate_box_liq_list, ddof=1)

    pressure_mean_box_vap = np.mean(pressure_repilcate_box_vap_list)
    total_molecules_mean_box_vap = np.mean(total_molecules_repilcate_box_vap_list)
    Rho_mean_box_vap = np.mean(Rho_repilcate_box_vap_list)
    volume_mean_box_vap = np.mean(volume_repilcate_box_vap_list)
    length_if_cube_mean_box_vap = np.mean(length_if_cube_repilcate_box_vap_list)
    Hv_mean_box_vap = np.mean(Hv_repilcate_box_vap_list)

    pressure_std_box_vap = np.std(pressure_repilcate_box_vap_list, ddof=1)
    total_molecules_std_box_vap = np.std(total_molecules_repilcate_box_vap_list, ddof=1)
    Rho_std_box_vap = np.std(Rho_repilcate_box_vap_list, ddof=1)
    volume_std_box_vap = np.std(volume_repilcate_box_vap_list, ddof=1)
    length_if_cube_std_box_vap = np.std(length_if_cube_repilcate_box_vap_list, ddof=1)
    Hv_std_box_vap = np.std(Hv_repilcate_box_vap_list, ddof=1)

    # *************************
    # get the replica means and std.devs (end)
    # *************************

    # ************************************
    # write the analysis data files for the liquid and vapor boxes (start)
    # ************************************

    box_liq_data_txt_file.write(
        f"{temp_mean: <30} "
        f"{temp_std: <30} "
        f"{pressure_mean_box_liq: <30} "
        f"{pressure_std_box_liq: <30} "
        f"{total_molecules_mean_box_liq: <30} "
        f"{total_molecules_std_box_liq: <30} "
        f"{Rho_mean_box_liq: <30} "
        f"{Rho_std_box_liq: <30} "
        f"{volume_mean_box_liq: <30} "
        f"{volume_std_box_liq: <30} "
        f"{length_if_cube_mean_box_liq: <30} "
        f"{length_if_cube_std_box_liq: <30} "
        f"{Hv_mean_box_liq: <30} "
        f"{Hv_std_box_liq: <30} "
        f" \n"
    )

    box_vap_data_txt_file.write(
        f"{temp_mean: <30} "
        f"{temp_std: <30} "
        f"{pressure_mean_box_vap: <30} "
        f"{pressure_std_box_vap: <30} "
        f"{total_molecules_mean_box_vap: <30} "
        f"{total_molecules_std_box_vap: <30} "
        f"{Rho_mean_box_vap: <30} "
        f"{Rho_std_box_vap: <30} "
        f"{volume_mean_box_vap: <30} "
        f"{volume_std_box_vap: <30} "
        f"{length_if_cube_mean_box_vap: <30} "
        f"{length_if_cube_std_box_vap: <30} "
        f"{Hv_mean_box_vap: <30} "
        f"{Hv_std_box_vap: <30} "
        f" \n"
    )
    # ************************************
    # write the analysis data files for the liquid and vapor boxes (end)
    # ************************************
"""

# ******************************************************
# ******************************************************
# data analysis - get the average and std. dev. from/across all the replicate (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# data analysis - get the critical and boiling point data for each set of replicates (start)
# ******************************************************
# ******************************************************

@aggregator.groupby(key=statepoint_without_temperature, sort_by="production_temperature_K", sort_ascending=True)
@Project.operation.with_directives(
     {
         "np": 1,
         "ngpu": 0,
         "memory": memory_needed,
         "walltime": walltime_gomc_analysis_hr,
     }
)
@FlowProject.pre(
     lambda * jobs: all(
         gomc_sim_completed_properly(job, gomc_production_control_file_name_str)
         for job in jobs
     )
)
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5b_analysis_replica_averages_completed)
@Project.post(part_5c_analysis_critical_and_boiling_points_replicate_data_completed)
def part_5c_analysis_critical_and_boiling__points_replicate_data(*jobs):
    # ***************************************************
    #  user changable variables (start)
    # ***************************************************
    lowest_Tr_for_critical_calcs = 0.8

    Beta_for_critical_points = 0.326

    # ***************************************************
    #  user changable variables (end)
    # ***************************************************

    # ***************************************************
    #  create the required lists and file labels for the replicates (start)
    # ***************************************************
    input_column_temp_title = 'temp_K'  # column title title for temp
    input_column_no_pressure_title = 'P_bar'  # column title title for PRESSURE
    input_column_Rho_title = 'Rho_kg_per_m_cubed'  # column title title for TOT_DENS
    input_column_box_Hv_kJ_per_mol_title = 'Hv_kJ_per_mol'  # column title title for HEAT_VAP

    output_column_Tc_K_title = 'Tc_K'
    output_column_Rho_c_kg_per_m_cubed_title = 'Rho_c_kg_per_m_cubed'
    output_column_Pc_bar_title = 'Pc_bar'
    output_column_Hv_kJ_per_mol_298_15K_title = 'Hv_kJ_per_mol_298_15K'
    output_column_lowest_T_K_for_Tc = "lowest_T_K_for_Tc"
    output_column_highest_T_K_for_Tc = "highest_T_K_for_Tc"
    output_column_lowest_Tr_K_for_Tc = "lowest_Tr_K_for_Tc"
    output_column_highest_Tr_K_for_Tc = "lowest_Tr_K_for_Tc"
    output_column_No_T_K_for_Tc = "No_T_K_for_Tc"

    # write the critical data file
    output_critical_data_replicate_txt_file_header = f"{output_column_Tc_K_title: <30} " \
                                                     f"{output_column_Rho_c_kg_per_m_cubed_title: <30} " \
                                                     f"{output_column_Pc_bar_title: <30} " \
                                                     f"{output_column_Hv_kJ_per_mol_298_15K_title: <30} " \
                                                     f"{output_column_lowest_T_K_for_Tc: <30} " \
                                                     f"{output_column_highest_T_K_for_Tc: <30} " \
                                                     f"{output_column_lowest_Tr_K_for_Tc: <30} " \
                                                     f"{output_column_highest_Tr_K_for_Tc: <30} " \
                                                     f"{output_column_No_T_K_for_Tc: <30} " \
                                                     f"\n"


    write_critial_file_path_and_name = f'analysis/{output_critical_data_replicate_txt_file_name}'
    if os.path.isfile(write_critial_file_path_and_name):
        critial_data_replicate_txt_file = open(write_critial_file_path_and_name, "a")
    else:
        critial_data_replicate_txt_file = open(write_critial_file_path_and_name, "w")
        critial_data_replicate_txt_file.write(output_critical_data_replicate_txt_file_header)

    # ***************************************************
    #  create the required lists and file labels for the replicates (end)
    # ***************************************************

    # *************************
    # drawing in data from the avg and std dev data file (start)
    # *************************
    temp_avg_box_liq_list = []
    Rho_avg_box_liq_list = []

    temp_K_avg_box_vap_list = []
    pressure_bar_avg_box_vap_list = []
    Rho_avg_box_vap_list = []
    Hv_kJ_per_mol_avg_box_vap_list = []
    for job in jobs:


        reading_file_avg_std_liq = job.fn(f'{output_replicate_txt_file_name_liq}')
        data_avg_std_liq = pd.read_csv(reading_file_avg_std_liq, sep='\s+', header=0, na_values='NaN', index_col=False)

        temp_avg_box_liq_list.append(data_avg_std_liq.loc[:, input_column_temp_title][0])
        Rho_avg_box_liq_list.append( data_avg_std_liq.loc[:, input_column_Rho_title][0])


        reading_file_avg_std_vap = job.fn(f'{output_replicate_txt_file_name_vap}')
        data_avg_std_vap = pd.read_csv(reading_file_avg_std_vap, sep='\s+', header=0, na_values='NaN', index_col=False)

        temp_K_avg_box_vap_list.append(data_avg_std_vap.loc[:, input_column_temp_title][0])
        pressure_bar_avg_box_vap_list.append(data_avg_std_vap.loc[:, input_column_no_pressure_title][0])
        Rho_avg_box_vap_list.append(data_avg_std_vap.loc[:, input_column_Rho_title][0])
        Hv_kJ_per_mol_avg_box_vap_list.append(data_avg_std_vap.loc[:, input_column_box_Hv_kJ_per_mol_title][0])


    # *************************
    # drawing in data from the avg and std dev data file (end)
    # *************************

    Tc_K_replicate_list = []
    Rho_c_kg_per_m_cubed_replicate_list = []
    Pc_bar_replicate_list = []
    Hv_kJ_per_mol_298_15K_replicate_list = []

    lowest_T_K_for_Tc_replicate_list  = []
    highest_T_K_for_Tc_replicate_list  = []
    No_T_K_for_Tc_replicate_list  = []
    lowest_Tr_K_for_Tc_replicate_list  = []
    highest_Tr_K_for_Tc_replicate_list = []

    # check and make sure the liq and vap temps are in order for comparison
    for t_i in range(0, len(temp_avg_box_liq_list)):
        if temp_avg_box_liq_list[t_i] != temp_K_avg_box_vap_list[t_i]:
            raise ValueError("ERROR: the temperatures in the vapor and liquid average "
                             "and std. dev. data do not match in order.")

    # get the full list of temps for the critical point calcs
    starting_temps_K_for_critical_point_list = temp_avg_box_liq_list


    Rho_avg_box_liq_and_1_list = []
    for Rho_i in range(0, len(Rho_avg_box_liq_list)):
        Rho_avg_box_liq_and_1_list.append(
            (Rho_avg_box_liq_list[Rho_i] + Rho_avg_box_vap_list[Rho_i]) / 2
        )



    # ***************************************************
    #  create the required lists and file labels for the replicates (end)
    # ***************************************************

    # ***********************
    # calc the Critical points (start)
    # ***********************
    # loop thru to find lowers temps to use for critical calcs
    # set use_last_no_temps_for_critical_calcs to 1 as adding 1 on first for loop
    use_last_no_temps_for_critical_calcs = 1
    for temp_i in range(0, len(starting_temps_K_for_critical_point_list)):
        use_last_no_temps_for_critical_calcs += 1

        if len(starting_temps_K_for_critical_point_list) >= use_last_no_temps_for_critical_calcs:

            # (Rho_box_liq + Rho_box_liq)/2 = (Rho_critical) + A *(T -Tc)
            # (Rho_box_liq - Rho_box_liq) = B * |T -Tc |^(beta)
            # eq1 : Yr = ar + br * Xr = (Rho_box_liq + Rho_box_liq)/2 = (Rho_critical - A Tc) + A * T
            # eq2 : Ys =  as +bs *Xs = (Rho_box_liq - Rho_box_liq) = B^(1/beta) * Tc + (-B^(1/beta)) * T
            # Yr = (Rho_box_liq + Rho_box_liq)/2
            # Xr = T
            # ar = Rho_critical - A * Tc
            # br = A
            # Ys = (Rho_box_liq - Rho_box_liq) = B^(1/beta)
            # Xs = T
            # as = B^(1/beta)
            # bs = -B^(1/beta)

            temp_for_critical_list = []
            Rho_0_minus_Rho_1_for_critical_list = []
            Rho_0_plus_Rho_1_for_critical_list = []
            Rho_0_plus_Rho_1_div2_for_critical_list = []
            Rho_0_minus_Rho_1_power_1_div_Beta_for_critical_list = []

            ln_P_per_bar_for_critical_list = []
            T_power_neg1_for_critical_list = []
            Hv_kJ_per_mol_box_vap_for_critical_list = []

            # eq3: ln(P1/bar) = -deltaHv/(R*T1) + F
            # eq3: Hv_kJ/mol = W + Q * T

            R_kJ_per_mol_K = 0.008134
            T_298_15K = 298.15


            if use_last_no_temps_for_critical_calcs > len(starting_temps_K_for_critical_point_list):
                print('Error: Larger number of temperatures for the critical '
                      'point calculations than are available.'
                      )

            elif use_last_no_temps_for_critical_calcs <= 0:
                print('Error: The number of temperatures for the critical point '
                      'calculations negative or zero.'
                      )

            else:
                for k in range(len(starting_temps_K_for_critical_point_list) - use_last_no_temps_for_critical_calcs,
                               len(starting_temps_K_for_critical_point_list)
                               ):
                    # concentration 1 files and inputs
                    temp_value_iter = starting_temps_K_for_critical_point_list[k]
                    Rho_avg_box_liq_list_iter = Rho_avg_box_liq_list[k]
                    Rho_avg_box_vap_list_iter = Rho_avg_box_vap_list[k]
                    Rho_0_minus_Rho_1_iter = Rho_avg_box_liq_list_iter - Rho_avg_box_vap_list_iter
                    Rho_0_plus_Rho_1_iter = Rho_avg_box_liq_list_iter + Rho_avg_box_vap_list_iter
                    Rho_0_plus_Rho_1_div2_iter = \
                        (Rho_avg_box_liq_list_iter + Rho_avg_box_vap_list_iter) / 2

                    temp_for_critical_list.append(temp_value_iter)
                    Rho_0_minus_Rho_1_for_critical_list.append(Rho_0_minus_Rho_1_iter)
                    Rho_0_minus_Rho_1_power_1_div_Beta_for_critical_list.append(
                        Rho_0_minus_Rho_1_iter ** (1 / Beta_for_critical_points))
                    Rho_0_plus_Rho_1_for_critical_list.append(Rho_0_plus_Rho_1_iter)
                    Rho_0_plus_Rho_1_div2_for_critical_list.append(Rho_0_plus_Rho_1_div2_iter)

                    Hv_kJ_per_mol_box_vap_critical_iter = Hv_kJ_per_mol_avg_box_vap_list[k]
                    Hv_kJ_per_mol_box_vap_for_critical_list.append(Hv_kJ_per_mol_box_vap_critical_iter)

                    pressure_box_vap_iter = pressure_bar_avg_box_vap_list[k]
                    ln_P_per_bar_for_critical_list.append(np.log(pressure_box_vap_iter))
                    T_power_neg1_for_critical_list.append(1 / (starting_temps_K_for_critical_point_list[k]))

                    Rho_avg_box_liq_and_1_iter = \
                        (Rho_avg_box_liq_list_iter + Rho_avg_box_vap_list_iter) / 2
                    Rho_avg_box_liq_and_1_list.append(Rho_avg_box_liq_and_1_iter)

                eq1_slope, eq1_intercept, eq1_r_value, eq1_p_value, eq1_std_err = \
                    stats.linregress(temp_for_critical_list, Rho_0_plus_Rho_1_div2_for_critical_list)
                eq2_slope, eq2_intercept, eq2_r_value, eq2_p_value, eq2_std_err = \
                    stats.linregress(temp_for_critical_list, Rho_0_minus_Rho_1_power_1_div_Beta_for_critical_list)

                eq3_slope, eq3_intercept, eq3_r_value, eq3_p_value, eq3_std_err = \
                    stats.linregress(T_power_neg1_for_critical_list, ln_P_per_bar_for_critical_list)
                eq4_slope, eq4_intercept, eq4_r_value, eq4_p_value, eq4_std_err = \
                    stats.linregress(temp_for_critical_list, Hv_kJ_per_mol_box_vap_for_critical_list)

                Tc_K_eval = -eq2_intercept / eq2_slope
                if np.isnan(Tc_K_eval) == True:
                    Tc_K = 'NA'
                else:
                    Tc_K = Tc_K_eval

                Rho_c_kg_per_m_cubed_eval = eq1_intercept + eq1_slope * Tc_K
                if np.isnan(Rho_c_kg_per_m_cubed_eval) == True:
                    Rho_c_kg_per_m_cubed = 'NA'
                else:
                    Rho_c_kg_per_m_cubed = Rho_c_kg_per_m_cubed_eval

                Pc_bar_eval = np.exp(eq3_intercept + eq3_slope * (1 / Tc_K))
                if np.isnan(Pc_bar_eval) == True:
                    Pc_bar = 'NA'
                else:
                    Pc_bar = Pc_bar_eval

                Delta_Hv_kJ_per_mol_at_298_15K_eval = eq4_slope * T_298_15K + eq4_intercept
                if np.isnan(Delta_Hv_kJ_per_mol_at_298_15K_eval) == True:
                    Delta_Hv_kJ_per_mol_at_298_15K = 'NA'
                else:
                    Delta_Hv_kJ_per_mol_at_298_15K = Delta_Hv_kJ_per_mol_at_298_15K_eval


                # check if passes the lowest_Tr_for_critical_calcs before writing
                lowest_T_K_for_Tc_iter = temp_for_critical_list[0]
                highest_T_K_for_Tc_iter = temp_for_critical_list[-1]
                No_T_K_for_Tc_iter = len(temp_for_critical_list)
                lowest_Tr_K_for_Tc_iter = float(lowest_T_K_for_Tc_iter / Tc_K)
                highest_Tr_K_for_Tc_iter = float(highest_T_K_for_Tc_iter / Tc_K)


                if lowest_Tr_K_for_Tc_iter >= lowest_Tr_for_critical_calcs:
                    print('utilized temperatures used for critical calcs = ' + str(temp_for_critical_list))
                    Tc_K_replicate_list.append(Tc_K)
                    Rho_c_kg_per_m_cubed_replicate_list.append(Rho_c_kg_per_m_cubed)
                    Pc_bar_replicate_list.append(Pc_bar)
                    Hv_kJ_per_mol_298_15K_replicate_list.append(Delta_Hv_kJ_per_mol_at_298_15K)

                    lowest_T_K_for_Tc_replicate_list.append(lowest_T_K_for_Tc_iter)
                    highest_T_K_for_Tc_replicate_list.append(highest_T_K_for_Tc_iter)
                    No_T_K_for_Tc_replicate_list.append(No_T_K_for_Tc_iter)
                    lowest_Tr_K_for_Tc_replicate_list.append(lowest_Tr_K_for_Tc_iter)
                    highest_Tr_K_for_Tc_replicate_list.append(highest_Tr_K_for_Tc_iter)

    # write the critical data
    critial_data_replicate_txt_file.write(
        f"{Tc_K_replicate_list[-1]: <30} "
        f"{Rho_c_kg_per_m_cubed_replicate_list[-1]: <30} "
        f"{Pc_bar_replicate_list[-1]: <30} "
        f"{Hv_kJ_per_mol_298_15K_replicate_list[-1]: <30} "
        f"{lowest_T_K_for_Tc_replicate_list[-1]: <30} "
        f"{highest_T_K_for_Tc_replicate_list[-1]: <30} "
        f"{lowest_Tr_K_for_Tc_replicate_list[-1]: <30} "
        f"{highest_Tr_K_for_Tc_replicate_list[-1]: <30} "
        f"{No_T_K_for_Tc_replicate_list[-1]: <30} "
        f" \n"
    )
    # ***********************
    # calc the Critical points (end)
    # ***********************

    # ***********************
    # calc the Boiling points (end)
    # ***********************

    # set the standard pressure for the boiling point (bp) via the standard temperature and pressure (STP)
    Pbp_STP_bar = 1.01325
    ln_Pbp_STP_bar = np.log(Pbp_STP_bar)

    # ***************************************************
    #  create the required lists and file labels for the replicates (start)
    # ***************************************************
    output_column_Tbp_K_title = 'Tbp_K'
    output_column_Pbp_bar_title = 'Pbp_bar' # this is atom pressure (1.01325 bar)
    output_column_lowest_T_K_for_Tbp = "lowest_T_K_for_Tbp"
    output_column_highest_T_K_for_Tbp = "highest_T_K_for_Tbp"
    output_column_lowest_Tr_K_for_Tbp = "lowest_Tr_K_for_Tbp"
    output_column_highest_Tr_K_for_Tbp = "lowest_Tr_K_for_Tbp"
    output_column_No_T_K_for_Tbp = "No_T_K_for_Tbp"

    # write the boiling data file
    output_boiling_data_replicate_txt_file_header = f"{output_column_Tbp_K_title: <30} " \
                                                     f"{output_column_Pbp_bar_title: <30} " \
                                                     f"{output_column_lowest_T_K_for_Tbp: <30} " \
                                                     f"{output_column_highest_T_K_for_Tbp: <30} " \
                                                     f"{output_column_lowest_Tr_K_for_Tbp: <30} " \
                                                     f"{output_column_highest_Tr_K_for_Tbp: <30} " \
                                                     f"{output_column_No_T_K_for_Tbp: <30} " \
                                                     f"\n"

    write_boiling_file_path_and_name = f'analysis/{output_boiling_data_replicate_txt_file_name}'
    if os.path.isfile(write_boiling_file_path_and_name):
        boiling_data_replicate_txt_file = open(write_boiling_file_path_and_name, "a")
    else:
        boiling_data_replicate_txt_file = open(write_boiling_file_path_and_name, "w")
        boiling_data_replicate_txt_file.write(output_boiling_data_replicate_txt_file_header)
    # ***************************************************
    #  create the required lists and file labels for the replicates (end)
    # ***************************************************
    # get the full list of temps, pressure, and ln(pressure) for the critical point calcs
    temps_K_for_boiling_point_list = temp_K_avg_box_vap_list
    inverse_temps_K_for_boiling_point_list = [1/ temp_bp_i for temp_bp_i in temps_K_for_boiling_point_list]
    pressures_bar_for_boiling_point_list = pressure_bar_avg_box_vap_list
    ln_pressures_bar_for_boiling_point_list = [np.log(press_bp_i) for press_bp_i in pressure_bar_avg_box_vap_list]
    print('utilized temperatures used for boiling point calcs = ' + str(temps_K_for_boiling_point_list))

    lowest_T_K_for_Tbp_replicate = temps_K_for_boiling_point_list[0]
    highest_T_K_for_Tbp_replicate = temps_K_for_boiling_point_list[-1]
    No_T_K_for_Tbp_replicate = len(temps_K_for_boiling_point_list)
    lowest_Tr_K_for_Tbp_replicate = float(temps_K_for_boiling_point_list[0] / Tc_K_replicate_list[-1])
    highest_Tr_K_for_Tbp_replicate = float(temps_K_for_boiling_point_list[-1] / Tc_K_replicate_list[-1])


    eq_bp_slope, eq_bp_intercept, eq_bp_r_value, eq_bp_p_value, eq_bp_std_err = \
        stats.linregress(ln_pressures_bar_for_boiling_point_list, inverse_temps_K_for_boiling_point_list)

    Tbp_iter = (eq_bp_slope * ln_Pbp_STP_bar + eq_bp_intercept)**(-1)

    # write the critical data
    boiling_data_replicate_txt_file.write(
        f"{Tbp_iter: <30} "
        f"{Pbp_STP_bar: <30} "
        f"{lowest_T_K_for_Tbp_replicate: <30} "
        f"{highest_T_K_for_Tbp_replicate: <30} "
        f"{lowest_Tr_K_for_Tbp_replicate: <30} "
        f"{highest_Tr_K_for_Tbp_replicate: <30} "
        f"{No_T_K_for_Tbp_replicate: <30} "
        f" \n"
    )

    # ***********************
    # calc the Boiling points (end)
    # ***********************


# ******************************************************
# ******************************************************
# data analysis - get the critical and boilin point data for each set of replicates (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# data analysis - get the critical and boilin point data avg and std. dev across the replicates (start)
# ******************************************************
# ******************************************************
@aggregator.groupby(key=statepoint_without_temperature, sort_by="production_temperature_K", sort_ascending=True)
@Project.operation.with_directives(
     {
         "np": 1,
         "ngpu": 0,
         "memory": memory_needed,
         "walltime": walltime_gomc_analysis_hr,
     }
)
@FlowProject.pre(
     lambda * jobs: all(
         gomc_sim_completed_properly(job, gomc_production_control_file_name_str)
         for job in jobs
     )
)
@Project.pre(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.pre(part_5b_analysis_replica_averages_completed)
@Project.pre(part_5c_analysis_critical_and_boiling_points_replicate_data_completed)
@Project.post(part_5d_analysis_critical_and_boiling_points_avg_std_data_completed)
def part_5d_analysis_critical_and_boiling_points_avg_std_data(*jobs):
    # ***********************
    # calc the Critical points (start)
    # ***********************

    output_column_Tc_K_title = 'Tc_K'
    output_column_Rho_c_kg_per_m_cubed_title = 'Rho_c_kg_per_m_cubed'
    output_column_Pc_bar_title = 'Pc_bar'
    output_column_Hv_kJ_per_mol_298_15K_title = 'Hv_kJ_per_mol_298_15K'
    output_column_lowest_T_K_for_Tc = "lowest_T_K_for_Tc"
    output_column_highest_T_K_for_Tc = "highest_T_K_for_Tc"
    output_column_lowest_Tr_K_for_Tc = "lowest_Tr_K_for_Tc"
    output_column_highest_Tr_K_for_Tc = "lowest_Tr_K_for_Tc"
    output_column_No_T_K_for_Tc = "No_T_K_for_Tc"

    output_column_Tc_K_std_title = 'Tc_std_K'
    output_column_Rho_c_kg_per_m_cubed_std_title = 'Rho_c_std_kg_per_m_cubed'
    output_column_Pc_bar_std_title = 'Pc_std_bar'
    output_column_Hv_kJ_per_mol_298_15K_std_title = 'Hv_std_kJ_per_mol_298_15K'

    output_critical_data_avg_std_txt_file_header = f"{output_column_Tc_K_title: <30} " \
                                                   f"{output_column_Tc_K_std_title: <30} " \
                                                   f"{output_column_Rho_c_kg_per_m_cubed_title: <30} " \
                                                   f"{output_column_Rho_c_kg_per_m_cubed_std_title: <30} " \
                                                   f"{output_column_Pc_bar_title: <30} " \
                                                   f"{output_column_Pc_bar_std_title: <30} " \
                                                   f"{output_column_Hv_kJ_per_mol_298_15K_title: <30} " \
                                                   f"{output_column_Hv_kJ_per_mol_298_15K_std_title: <30} " \
                                                   f"{output_column_lowest_T_K_for_Tc: <30} " \
                                                   f"{output_column_highest_T_K_for_Tc: <30} " \
                                                   f"{output_column_lowest_Tr_K_for_Tc: <30} " \
                                                   f"{output_column_highest_Tr_K_for_Tc: <30} " \
                                                   f"{output_column_No_T_K_for_Tc: <30} " \
                                                   f"\n"

    for job in jobs:
        if job.isfile(f'../../analysis/{output_critical_data_replicate_txt_file_name}'):
                reading_file_final_critical_replicate_data = job.fn(
                    f'../../analysis/{output_critical_data_replicate_txt_file_name}')
                data_final_critical_replicate = pd.read_csv(reading_file_final_critical_replicate_data,
                                                            sep='\s+',
                                                            header=0,
                                                            na_values='NaN',
                                                            index_col=False
                                                            )

                Tc_K_final_replicate_data_list = \
                    data_final_critical_replicate.loc[:, output_column_Tc_K_title]
                Rho_avg_box_vap_final_replicate_list = \
                    data_final_critical_replicate.loc[:, output_column_Rho_c_kg_per_m_cubed_title]
                Pc_bar_final_replicate_data_list = \
                    data_final_critical_replicate.loc[:, output_column_Pc_bar_title]
                Hv_kJ_per_mol_298_15K_final_replicate_list = \
                    data_final_critical_replicate.loc[:, output_column_Hv_kJ_per_mol_298_15K_title]

                lowest_T_K_for_Tc_final_replicate_data_list = \
                    data_final_critical_replicate.loc[:, output_column_lowest_T_K_for_Tc]
                highest_T_K_for_Tc_final_replicate_data_list = \
                    data_final_critical_replicate.loc[:, output_column_highest_T_K_for_Tc]
                lowest_Tr_K_for_Tc_final_replicate_list = \
                    data_final_critical_replicate.loc[:, output_column_lowest_Tr_K_for_Tc]
                highest_Tr_K_for_Tc_final_replicate_list = \
                    data_final_critical_replicate.loc[:, output_column_highest_Tr_K_for_Tc]
                No_T_K_for_Tc_final_replicate_list = \
                    data_final_critical_replicate.loc[:, output_column_No_T_K_for_Tc]

                # check that all replicates have same Temps used for the analysis:
                for temp_iter_k in range(0, len(lowest_T_K_for_Tc_final_replicate_data_list)):
                    if len(lowest_T_K_for_Tc_final_replicate_data_list) >= 1 \
                            and len(highest_T_K_for_Tc_final_replicate_data_list) >= 1 \
                            and len(lowest_Tr_K_for_Tc_final_replicate_list) >= 1 \
                            and len(highest_Tr_K_for_Tc_final_replicate_list) >= 1 \
                            and len(No_T_K_for_Tc_final_replicate_list) >= 1:

                        if not math.isclose(lowest_T_K_for_Tc_final_replicate_data_list[0],
                                            lowest_T_K_for_Tc_final_replicate_data_list[temp_iter_k],
                                            rel_tol=0, abs_tol=1e-9
                                            ) \
                                and not math.isclose(highest_T_K_for_Tc_final_replicate_data_list[0],
                                                     highest_T_K_for_Tc_final_replicate_data_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(lowest_Tr_K_for_Tc_final_replicate_list[0],
                                                     lowest_Tr_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(highest_Tr_K_for_Tc_final_replicate_list[0],
                                                     highest_Tr_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(No_T_K_for_Tc_final_replicate_list[0],
                                                     No_T_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ):
                            raise ValueError("ERROR: In the critical point data, "
                                             "the replicate data temperatures or number of points used is "
                                             "different for the replicates, and they need to be the same.  ")

                    else:
                        raise ValueError("ERROR: There is not 1 or any replicate data for the critical point analysis.")

                critial_data_avg_std_txt_file = open(
                    f'analysis/{output_critical_data_avg_std_of_replicates_txt_file_name}', "w"
                )
                critial_data_avg_std_txt_file.write(output_critical_data_avg_std_txt_file_header)

                Tc_K_mean = np.mean(Tc_K_final_replicate_data_list)
                Rho_c_kg_per_m_cubed_mean = np.mean(Rho_avg_box_vap_final_replicate_list)
                Pc_bar_mean = np.mean(Pc_bar_final_replicate_data_list)
                Hv_kJ_per_mol_298_15K_mean = np.mean(Hv_kJ_per_mol_298_15K_final_replicate_list)

                Tc_K_std = np.std(Tc_K_final_replicate_data_list, ddof=1)
                Rho_c_kg_per_m_cubed_std = np.std(Rho_avg_box_vap_final_replicate_list, ddof=1)
                Pc_bar_std = np.std(Pc_bar_final_replicate_data_list, ddof=1)
                Hv_kJ_per_mol_298_15K_std = np.std(Hv_kJ_per_mol_298_15K_final_replicate_list, ddof=1)

                lowest_T_K_for_Tc = lowest_T_K_for_Tc_final_replicate_data_list[0]
                highest_T_K_for_Tc = highest_T_K_for_Tc_final_replicate_data_list[0]
                No_T_K_for_Tc = lowest_Tr_K_for_Tc_final_replicate_list[0]
                lowest_Tr_K_for_Tc = highest_Tr_K_for_Tc_final_replicate_list[0]
                highest_Tr_K_for_Tc = No_T_K_for_Tc_final_replicate_list[0]

                critial_data_avg_std_txt_file.write(
                    f"{Tc_K_mean: <30} "
                    f"{Tc_K_std: <30} "
                    f"{Rho_c_kg_per_m_cubed_mean: <30} "
                    f"{Rho_c_kg_per_m_cubed_std: <30} "
                    f"{Pc_bar_mean: <30} "
                    f"{Pc_bar_std: <30} "
                    f"{Hv_kJ_per_mol_298_15K_mean: <30} "
                    f"{Hv_kJ_per_mol_298_15K_std: <30} "
                    f"{lowest_T_K_for_Tc: <30} "
                    f"{highest_T_K_for_Tc: <30} "
                    f"{lowest_Tr_K_for_Tc: <30} "
                    f"{highest_Tr_K_for_Tc: <30} "
                    f"{No_T_K_for_Tc: <30} "
                    f" \n"
                )
    # ***********************
    # calc the Critical points (end)
    # ***********************

    # ***********************
    # calc the Boiling points (start)
    # ***********************
    output_column_Tbp_K_title = 'Tbp_K'
    output_column_Tbp_K_std_title = 'Tbp_std_K'
    output_column_Pbp_bar_title = 'Pbp_bar'  # this is atom pressure (1.01325 bar)
    output_column_lowest_T_K_for_Tbp = "lowest_T_K_for_Tbp"
    output_column_highest_T_K_for_Tbp = "highest_T_K_for_Tbp"
    output_column_lowest_Tr_K_for_Tbp = "lowest_Tr_K_for_Tbp"
    output_column_highest_Tr_K_for_Tbp = "lowest_Tr_K_for_Tbp"
    output_column_No_T_K_for_Tbp = "No_T_K_for_Tbp"

    # write the boiling data file
    output_boiling_data_avg_std_txt_file_header = f"{output_column_Tbp_K_title: <30} " \
                                                  f"{output_column_Tbp_K_std_title : <30} " \
                                                  f"{output_column_Pbp_bar_title: <30} " \
                                                  f"{output_column_lowest_T_K_for_Tbp: <30} " \
                                                  f"{output_column_highest_T_K_for_Tbp: <30} " \
                                                  f"{output_column_lowest_Tr_K_for_Tbp: <30} " \
                                                  f"{output_column_highest_Tr_K_for_Tbp: <30} " \
                                                  f"{output_column_No_T_K_for_Tbp: <30} " \
                                                  f"\n"

    for job in jobs:
        if job.isfile(f'../../analysis/{output_boiling_data_replicate_txt_file_name}'):
                reading_file_final_boiling_replicate_data = job.fn(
                    f'../../analysis/{output_boiling_data_replicate_txt_file_name}')
                data_final_boiling_replicate = pd.read_csv(reading_file_final_boiling_replicate_data,
                                                            sep='\s+',
                                                            header=0,
                                                            na_values='NaN',
                                                            index_col=False
                                                            )

                Tbp_K_final_replicate_data_list = \
                    data_final_boiling_replicate.loc[:, output_column_Tbp_K_title]
                Pbp_bar_final_replicate_data_list = \
                    data_final_boiling_replicate.loc[:, output_column_Pbp_bar_title]

                lowest_T_K_for_Tbp_final_replicate_data_list = \
                    data_final_boiling_replicate.loc[:, output_column_lowest_T_K_for_Tbp]
                highest_T_K_for_Tbp_final_replicate_data_list = \
                    data_final_boiling_replicate.loc[:, output_column_highest_T_K_for_Tbp]
                lowest_Tr_K_for_Tbp_final_replicate_list = \
                    data_final_boiling_replicate.loc[:, output_column_lowest_Tr_K_for_Tbp]
                highest_Tr_K_for_Tbp_final_replicate_list = \
                    data_final_boiling_replicate.loc[:, output_column_highest_Tr_K_for_Tbp]
                No_T_K_for_Tbp_final_replicate_list = \
                    data_final_boiling_replicate.loc[:, output_column_No_T_K_for_Tbp]

                # check that all replicates have same Temps used for the analysis:
                for temp_iter_k in range(0, len(lowest_T_K_for_Tbp_final_replicate_data_list)):
                    if len(lowest_T_K_for_Tbp_final_replicate_data_list) >= 1 \
                            and len(highest_T_K_for_Tbp_final_replicate_data_list) >= 1 \
                            and len(lowest_Tr_K_for_Tbp_final_replicate_list) >= 1 \
                            and len(highest_Tr_K_for_Tbp_final_replicate_list) >= 1 \
                            and len(No_T_K_for_Tbp_final_replicate_list) >= 1:

                        if not math.isclose(lowest_T_K_for_Tbp_final_replicate_data_list[0],
                                            lowest_T_K_for_Tbp_final_replicate_data_list[temp_iter_k],
                                            rel_tol=0, abs_tol=1e-9
                                            ) \
                                and not math.isclose(Pbp_bar_final_replicate_data_list[0],
                                                     Pbp_bar_final_replicate_data_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(highest_T_K_for_Tbp_final_replicate_data_list[0],
                                                     highest_T_K_for_Tc_final_replicate_data_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(lowest_Tr_K_for_Tbp_final_replicate_list[0],
                                                     lowest_Tr_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(highest_Tr_K_for_Tbp_final_replicate_list[0],
                                                     highest_Tr_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ) \
                                and not math.isclose(No_T_K_for_Tbp_final_replicate_list[0],
                                                     No_T_K_for_Tc_final_replicate_list[temp_iter_k],
                                                     rel_tol=0, abs_tol=1e-9
                                                     ):
                            raise ValueError("ERROR: In the boiling point data, "
                                             "the replicate data temperatures or number of points used is "
                                             "different for the replicates, and they need to be the same.  ")

                    else:
                        raise ValueError("ERROR: There is not 1 or any replicate data for the boiling point analysis.")

                boiling_data_avg_std_txt_file = open(
                    f'analysis/{output_boiling_data_avg_std_of_replicates_txt_file_name}', "w"
                )
                boiling_data_avg_std_txt_file.write(output_boiling_data_avg_std_txt_file_header)

                Tbp_K_mean = np.mean(Tbp_K_final_replicate_data_list)
                Tbp_K_std = np.std(Tbp_K_final_replicate_data_list, ddof=1)

                Pbp_bar = np.mean(Pbp_bar_final_replicate_data_list)

                lowest_T_K_for_Tbp = lowest_T_K_for_Tbp_final_replicate_data_list[0]
                highest_T_K_for_Tbp = highest_T_K_for_Tbp_final_replicate_data_list[0]
                No_T_K_for_Tbp = lowest_Tr_K_for_Tbp_final_replicate_list[0]
                lowest_Tr_K_for_Tbp = highest_Tr_K_for_Tbp_final_replicate_list[0]
                highest_Tr_K_for_Tbp = No_T_K_for_Tbp_final_replicate_list[0]

                boiling_data_avg_std_txt_file.write(
                    f"{Tbp_K_mean: <30} "
                    f"{Tbp_K_std: <30} "
                    f"{Pbp_bar: <30} "
                    f"{lowest_T_K_for_Tbp: <30} "
                    f"{highest_T_K_for_Tbp: <30} "
                    f"{lowest_Tr_K_for_Tbp: <30} "
                    f"{highest_Tr_K_for_Tbp: <30} "
                    f"{No_T_K_for_Tbp: <30} "
                    f" \n"
                )
    # ***********************
    # calc the Boiling points (end)
    # ***********************
# ******************************************************
# ******************************************************
# data analysis - get the critical and boiling point data avg and std. dev across the replicates (end)
# ******************************************************
# ******************************************************





# ******************************************************
# ******************************************************
# signac end code (start)
# ******************************************************
# ******************************************************
if __name__ == "__main__":
    pr = Project()
    pr.main()
# ******************************************************
# ******************************************************
# signac end code (end)
# ******************************************************
# ******************************************************
