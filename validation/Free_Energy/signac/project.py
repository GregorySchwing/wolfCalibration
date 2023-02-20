"""GOMC's setup for signac, signac-flow, signac-dashboard for this study."""
# project.py


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

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


class Grid(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""

    hostname_pattern = r".*\.grid\.wayne\.edu"
    template = "grid.sh"


class Potoff(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""

    hostname_pattern = r".*reslab32ai8111"
    template = "potoff.sh"




# ******************************************************
# users typical variables, but not all (start)
# ******************************************************
# set binary path to gomc binary files (the bin folder).
# If the gomc binary files are callable directly from the terminal without a path,
# please just enter and empty string (i.e., "" or '')

# WSU grid binary paths
gomc_binary_path = "/wsu/home/go/go24/go2432/wolfCalibration/validation/Free_Energy/signac/bin"
namd_binary_path = "/wsu/home/go/go24/go2432/wolfCalibration/validation/Free_Energy/signac/bin"


#gomc_binary_path = "/wsu/home/go/go24/go2432/wolfCalibrationLong/validation/Free_Energy/signac/bin"
#namd_binary_path = "/wsu/home/go/go24/go2432/wolfCalibrationLong/validation/Free_Energy/signac/bin"

# Potoff cluster bin paths
# Potoff cluster bin paths
#gomc_binary_path = "/home6/go2432/wolfCalibration/validation/Free_Energy/signac/bin"
#namd_binary_path = "/home6/go2432/wolfCalibration/validation/Free_Energy/signac/bin"

# local bin paths
#gomc_binary_path = "/home/greg/Documents/wolfCalibration/validation/Free_Energy/signac/bin"
#namd_binary_path = "/home/greg/Documents/wolfCalibration/validation/Free_Energy/signac/bin"

#WSL local bin paths
#gomc_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"
#namd_binary_path = "/mnt/c/Users/grego/OneDrive/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"

# brads workstation binary paths
gomc_binary_path = "/home/greg/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"
namd_binary_path = "/home/greg/Desktop/wolfCalibration/validation/Free_Energy/signac/bin"

# number of simulation steps
#"""
gomc_steps_equilb_design_ensemble = 30 * 10**6 # set value for paper = 10 * 10**6
precal_eq_gomc_steps =  gomc_steps_equilb_design_ensemble # set value for paper = 10 * 10**6

gomc_steps_lamda_production = 50 * 10**6 # set value for paper = 50 * 10**6
gomc_console_output_data_every_X_steps = 5 * 10**2 # set value for paper = 100 * 10**3
gomc_output_data_every_X_steps = 5 * 10**6 # set value for paper = 100 * 10**3
#gomc_free_energy_output_data_every_X_steps = 10 * 10**3 # set value for paper = 10 * 10**3

MC_steps = int(gomc_steps_equilb_design_ensemble)
EqSteps = 5 * 10**6
Calibration_MC_steps = 5 * 10**5
Calibration_MC_Eq_Steps = 5 * 10**4
Wolf_Sanity_MC_steps = 100 * 10**6 # set value for paper = 50 * 10**6
"""
# number of simulation steps
gomc_steps_equilb_design_ensemble = 5 * 10**3 # set value for paper = 10 * 10**6
precal_eq_gomc_steps =  5 * 10**3 # set value for paper = 10 * 10**6

gomc_steps_lamda_production = 5 * 10**3 # set value for paper = 50 * 10**6
gomc_console_output_data_every_X_steps = 5 * 10**2 # set value for paper = 100 * 10**3
gomc_output_data_every_X_steps = 5 * 10**3 # set value for paper = 100 * 10**3
#gomc_free_energy_output_data_every_X_steps = 10 * 10**3 # set value for paper = 10 * 10**3

MC_steps = int(gomc_steps_equilb_design_ensemble)
EqSteps = 1000
Calibration_MC_steps = 5 * 10**3
Calibration_MC_Eq_Steps = 1 * 10**3 
Wolf_Sanity_MC_steps = 1 * 10**4
"""
"""
During the
production run, the change in energy (DeltaU i,j ) between
the current lambda state and all other lambda states,
and the derivative of potential with respect to lambda
(dU Coul /dλ Coul , dU LJ /dλ LJ ), were evaluated and stored
for post-simulation analysis every 5 × 10 3 MCS.
"""
gomc_free_energy_output_data_every_X_steps = 5 * 10**3 # set value for paper = 10 * 10**3

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
namd_equilb_NPT_control_file_name_str = "namd_equilb_NPT"

# The equilb using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
gomc_equilb_design_ensemble_control_file_name_str = "gomc_equilb_design_ensemble"

# The production run using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
gomc_production_control_file_name_str = "gomc_production_run"


preliminary_output_replicate_txt_file_name_box_0 = "preliminary_analysis_avg_data_box_0.txt"
preliminary_uncorrelated_output_replicate_txt_file_name_box_0 = "preliminary_analysis_uncorrelated_data_avg_data_box_0.txt"
# Analysis (each replicates averages):
# Output text (txt) file names for each replicates averages
# directly put in each replicate folder (.txt, .dat, etc)
output_replicate_txt_file_name_box_0 = "analysis_avg_data_box_0.txt"
output_uncorrelated_replicate_txt_file_name_box_0 = "analysis_uncorrelated_avg_data_box_0.txt"

# Analysis (averages and std. devs. of  # all the replcates):
# Output text (txt) file names for the averages and std. devs. of all the replcates,
# including the extention (.txt, .dat, etc)
output_avg_std_of_replicates_txt_file_name_box_0 = "analysis_avg_std_of_replicates_box_0.txt"

preliminary_output_avg_std_of_replicates_txt_file_name_box_0 = "preliminary_analysis_avg_std_of_replicates_box_0.txt"


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
memory_needed = 4



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

def append_wolf_calibration_parameters(job, filename, iter_num):

    WolfCutoffBoxList = [0]

    WolfCutoffCoulombLowerBoundList = [10]
    WolfCutoffCoulombUpperBoundList = [15]
    WolfCutoffCoulombIntervalList = [0.1]

    WolfAlphaLowerBoundList = [0.0]
    WolfAlphabUpperBoundList = [0.5]
    WolfAlphaIntervalList = [0.005]

    wolfCalFreq = 1000

    with open(job.fn(filename), "a") as myfile:
        defPotLine = "InitStep\t{zero}\n".format(zero=0)
        myfile.write(defPotLine)
        defPotLine = "Wolf\t{freq}\n".format(freq=False)
        myfile.write(defPotLine)
        defPotLine = "WolfCalibrationFreq\tTrue\t{freq}\n".format(freq=wolfCalFreq)
        myfile.write(defPotLine)
        for box, wolfCutoffLower, wolfCutoffUpper, wolfCutoffInterval, wolfAlphaLower, wolfAlphaUpper, wolfAlphaInterval \
        in zip(WolfCutoffBoxList, WolfCutoffCoulombLowerBoundList, WolfCutoffCoulombUpperBoundList, WolfCutoffCoulombIntervalList, \
        WolfAlphaLowerBoundList, WolfAlphabUpperBoundList, WolfAlphaIntervalList):
            alphaLine = "WolfAlphaRange\t{box}\t{lb}\t{ub}\t{inter}\n".format(box=box, lb=wolfAlphaLower, ub=wolfAlphaUpper, inter=wolfAlphaInterval)
            myfile.write(alphaLine)


def append_checkpoint_line(job, config_file_name, path_to_previous_checkpoint_file):
    with open(job.fn("{}.conf".format(config_file_name)), "a") as myfile:
        checkpointLine = "Checkpoint\tTrue\t{}\n".format(path_to_previous_checkpoint_file)
        myfile.write(checkpointLine)


def append_default_wolf_parameters_line(job, config_file_name):
    with open(job.fn(config_file_name), "a") as myfile:
        defPotLine = "Wolf\t{freq}\n".format(freq=True)
        myfile.write(defPotLine)
        defPotLine = "WolfKind\t{freq}\n".format(freq=job.sp.wolf_model)
        myfile.write(defPotLine)
        defPotLine = "WolfPotential\t{freq}\n".format(freq=job.sp.wolf_potential)
        myfile.write(defPotLine)   
        defPotLine = "WolfAlpha\t0\t{freq}\n".format(freq=job.sp.alpha)
        myfile.write(defPotLine)   

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
def part_1b_under_equilb_design_ensemble_run_limit(job):
    """Check that the equilbrium design ensemble run is under it's run limit."""
    try:
        if (
            job.doc.calibration_iteration_number
            >= job.doc.calibration_iteration_number_max_number
        ):
            job.doc.calibration_iteration_number_under_limit = False
            return job.doc.calibration_iteration_number_under_limit

        else:
            return True
    except:
        return False


@Project.label
@flow.with_job
def part_4b_append_done(job):
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return job.isfile(job.doc.path_to_wolf_results_repl_0_dir+"append_done.txt")
    except:
        return False


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
    equilibration_ensemble = "NVT"
    #equilibration_ensemble = "NPT"
    production_ensemble = "NVT"

    # set the GOMC production ensemble temp, pressure, molecule, box dimenstion and residue names
    job.doc.equilibration_ensemble = equilibration_ensemble
    job.doc.production_ensemble = production_ensemble
    job.doc.production_pressure_bar = (1 * u.atm).to('bar')
    job.doc.production_temperature_K = job.sp.production_temperature_K

    """
    Liquid phase systems contained one solute in a solvent
    box of 200 1-octanol, 150 n-hexadecane, or 1000 water
    molecules. Initial cubic box sizes were selected to produce
    densities that were close to equilibrium, with a side length
    of 37.6, 41.6, and 31.3 Å for 1-octanol, n-hexadecane,
    and water, respectively.
    """
    """
    angstrom3 = (u.angstrom * u.angstrom * u.angstrom)
    cm3 = (u.cm * u.cm * u.cm)
    job.doc.volume = ((31.3 * u.angstrom) * (31.3 * u.angstrom) * (31.3 * u.angstrom)).to(cm3)
    from scipy import constants
    molar_mass_of_solvent = 18.01528 * u.mol
    job.doc.N_liquid_solvent = int((constants.Avogadro * job.sp.density * job.doc.volume )/ molar_mass_of_solvent)
    """
    job.doc.N_liquid_solvent = 1000

    print(job.doc.N_liquid_solvent)
    if (job.sp.solute == "solvent_box"):
        job.doc.N_liquid_solute = 0
    else:
        job.doc.N_liquid_solute = 1


    job.doc.liq_box_lengths_ang = 31.3 * u.angstrom

    if job.sp.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
        job.doc.Rcut_ang = 15 * u.angstrom  # this is the Rcut for GOMC it is the Rswitch for NAMD
    else:
        job.doc.Rcut_ang = 14 * u.angstrom  # this is the Rcut for GOMC it is the Rswitch for NAMD

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
    g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
    kg_per_m3 = u.kg / (u.m * u.m * u.m)


    #job.doc.density = (job.sp.density * g_per_cm3).to(kg_per_m3)

    job.doc.namd_node_ncpu = 4
    #job.doc.namd_node_ngpu = 1
    job.doc.namd_node_ngpu = 0

    job.doc.gomc_ncpu = 4  # 1 is optimal but I want data quick.  run time is set for 1 cpu
    #job.doc.gomc_ngpu = 1
    job.doc.gomc_ngpu = 0

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

    # set the initial iteration number of the simulation
    if equilibration_ensemble == "NPT":
        job.doc.namd_equilb_NPT_gomc_binary_file = f"namd2"
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NPT"
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_CPU_NPT"
        #job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_NPT"
    elif equilibration_ensemble == "NVT":
        job.doc.namd_equilb_NPT_gomc_binary_file = f"namd2"
        job.doc.gomc_equilb_design_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NVT"
        job.doc.gomc_calibration_gomc_binary_file = f"GOMC_CPU_NVT"
        #job.doc.gomc_calibration_gomc_binary_file = f"GOMC_GPU_NVT"

    else:
        raise ValueError(
            "ERROR: The 'GCMC', 'GEMC_NVT', 'GEMC_NPT' ensembles is not currently available for this project.py "
        )
    

    if production_ensemble == "NPT":
        job.doc.gomc_production_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NPT"
    elif production_ensemble == "NVT":
        job.doc.gomc_production_ensemble_gomc_binary_file = f"GOMC_{job.doc.gomc_cpu_or_gpu}_NVT"
    else:
        raise ValueError(
            "ERROR: The 'GCMC', 'GEMC_NVT', 'GEMC_NPT' ensembles is not currently available for this project.py "
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
    try:
        if (
            job.isfile(job.doc.namd_output_dir+f"{namd_ff_filename_str}.inp")
            and job.isfile(job.doc.namd_output_dir+f"{gomc_ff_filename_str}.inp")
            and job.isfile(
                job.doc.namd_output_dir+f"{mosdef_structure_box_0_name_str}.psf"
            )
            and job.isfile(
                job.doc.namd_output_dir+f"{mosdef_structure_box_0_name_str}.pdb"
            )
        ):
            file_written_bool = True

        return file_written_bool
    except:
        return False

# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and FF files were written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC control file was written (end)
# ******************************************************
# ******************************************************

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

def gomc_cal_completed_properly(job, control_filename_str):
    """General check to see if the gomc simulation was completed properly."""
    job_run_properly_bool = False
    output_log_file = "out_{}.dat".format(control_filename_str)
    if job.isfile(output_log_file):
        with open(job.fn(output_log_file), "r") as fp:
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

# ******************************************************
# ******************************************************
# check if GOMC and NAMD simulation are completed properly (start)
# ******************************************************
# ******************************************************
# function for checking if GOMC simulations are completed properly
def gomc_sseq_completed_properly(job, control_filename_str):
    """General check to see if the gomc simulation was completed properly."""
    job_run_properly_bool = False
    try:
        if job.isfile(job.doc.path_to_gomc_eq_console):
            with open(job.doc.path_to_gomc_eq_console, "r") as fp:
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
    except:
            return False



# function for checking if NAMD simulations are completed properly
def namd_sim_completed_properly(job, control_filename_str):
    """General check to see if the namd simulation was completed properly."""
    job_run_properly_bool = False

    try:
        if job.isfile(job.doc.path_to_namd_console):
            with open(job.doc.path_to_namd_console, "r") as fp:
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
def part_4a_job_namd_equilb_NPT_completed_properly(job):
    """Check to see if the  namd_equilb_NPT_control_file was completed properly
    (high temperature to set temperature NAMD control file)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return namd_sim_completed_properly(job, namd_equilb_NPT_control_file_name_str)
    except:
        return False


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_calibration_converged(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        if (job.doc.current_best_alpha == job.doc.previous_best_alpha):
            job.doc.calibration_converged = True
    except:
        return False

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_all_replicates_best_alpha_found(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return job.isfile(job.doc.path_to_wolf_results_repl_0_dir+"best_alpha_all_replicas.csv")
    except:
        return False

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_check_convergence(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return job.doc.check_convergence
    except:
        return False

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_calibration_needs_to_be_checked(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return not job.doc.check_convergence
    except:
        return True


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_wrote_alpha_csv(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try:
        return job.isfile(job.doc.path_to_wolf_results_my_repl_dir+"full_calibration_replica_{}_BOX_{}.csv".format(job.doc.replica_number_int, 0))
    except:
        return False



# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_sseq_average_obtained(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    try: 
        return job.isfile(job.doc.path_to_ew_results_my_repl_dir+"ewald_average.csv")
    except:
        return False

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_calibration_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    #This will cause Ewald sims to wait for Wolf calibration to complete.
    #if(job.sp.wolf_potential == "Results"):
    #    return True
    try: 
        output_name_control_file_name = "wolf_calibration_{}".format(job.doc.calibration_iteration_number)
        #output_name_control_file_name = "wolf_calibration"
        return gomc_sim_completed_properly(
            job,
            output_name_control_file_name,
        )
    except:
        return False

def part_4b_job_gomc_append_wolf_parameters_to_calibration(job):
    try:
        if (job.doc.calibration_iteration_number >= job.doc.calibration_iteration_number_max_number):
            job.doc.calibration_converged = True
            return True
        cols = ["MODEL", "POT", "ALPHA"]
        colList = ["ALPHA","WAIBEL2018_DSF", "RAHBARI_DSF", "WAIBEL2019_DSF", "WAIBEL2018_DSP", "RAHBARI_DSP", "WAIBEL2019_DSP"]
        print(job.doc.calibration_iteration_number)
        converged = True
        NUMBOXES = 1
        for b in range (NUMBOXES):
            curr = pd.DataFrame()
            prev = pd.DataFrame()
            if (job.isfile("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(job.doc.calibration_iteration_number, b))):
                curr = pd.read_csv (job.fn("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(job.doc.calibration_iteration_number, b)), delim_whitespace=True, header=None)
            else:
                return False
            if (job.doc.calibration_iteration_number != 0):
                if (job.isfile("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(job.doc.calibration_iteration_number-1, b))):
                    prev = pd.read_csv (job.fn("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(job.doc.calibration_iteration_number-1, b)), delim_whitespace=True, header=None)
                else:
                    return False

            if (job.isfile("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}.dat".format(job.doc.calibration_iteration_number, b))):
                modelRelErr = pd.read_csv (job.fn("wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}.dat".format(job.doc.calibration_iteration_number, b)), delim_whitespace=True, header=None, names=colList, index_col=0)
                print(modelRelErr)
                print(modelRelErr.abs().min().idxmin().split("_"))
                print(modelRelErr.abs().min().min())
            else:
                return False 

            print(prev)
            print(curr)
            print("Change in best alpha")
            changeBwRuns = pd.concat([curr,prev]).drop_duplicates(keep=False)
            print(changeBwRuns)
            print(changeBwRuns.empty)
            converged = converged and changeBwRuns.empty
            WolfDefaultKind = modelRelErr.abs().min().idxmin().split("_")[0]
            WolfDefaultPotential = modelRelErr.abs().min().idxmin().split("_")[1]
            mask = curr['MODEL'] == WolfDefaultKind
            mask2 = curr['POT'] == WolfDefaultPotential
            c = np.logical_and(mask, mask2)
            nextAlpha = curr[c]["ALPHA"].values[0]
            print(nextAlpha)
            print("Run next iteration {} with alpha {} ({} {})".format(job.doc.calibration_iteration_number+1,nextAlpha,WolfDefaultKind,WolfDefaultPotential))
            output_name_control_file_name = "wolf_calibration_{}.conf".format(
                job.doc.calibration_iteration_number+1
            )
            with open(job.fn(output_name_control_file_name), "a") as myfile:
                #defAlphaLine = "WolfKind\t{val}\n".format(val=WolfDefaultKind)
                #myfile.write(defAlphaLine)
                #defAlphaLine = "WolfPotential\t{val}\n".format(val=WolfDefaultPotential)
                #myfile.write(defAlphaLine)
                defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=b, val=nextAlpha)
                myfile.write(defAlphaLine)
        job.doc.calibration_converged = converged
        return converged   
    except:
        print(repr(e))
        return False


def part_4b_extract_best_initial_guess_from_ewald_calibration(job):
    try:
        cols = ["MODEL", "POT", "ALPHA"]
        colList = ["ALPHA","WAIBEL2018_DSF", "RAHBARI_DSF", "WAIBEL2019_DSF", "WAIBEL2018_DSP", "RAHBARI_DSP", "WAIBEL2019_DSP"]
        print(job.doc.calibration_iteration_number)
        bestA = []
        NUMBOXES = 1
        for b in range (NUMBOXES):
            curr = pd.DataFrame()
            if (job.isfile(job.doc.path_to_ew_results_repl_0_dir+"wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(0, b))):
                curr = pd.read_csv (job.doc.path_to_ew_results_repl_0_dir+"wolf_calibration_{}_WOLF_CALIBRATION_BOX_{}_BEST_ALPHAS.csv".format(0, b), delim_whitespace=True, header=None)
                curr.columns = cols
            else:
                return False

            print(curr)

            WolfDefaultKind = job.sp.wolf_model
            WolfDefaultPotential = job.sp.wolf_potential
            mask = curr['MODEL'] == WolfDefaultKind
            mask2 = curr['POT'] == WolfDefaultPotential
            c = np.logical_and(mask, mask2)
            nextAlpha = curr[c]["ALPHA"].values[0]
            print(nextAlpha)
            print("Run next iteration {} with alpha {} ({} {})".format(job.doc.calibration_iteration_number,nextAlpha,WolfDefaultKind,WolfDefaultPotential))
            output_name_control_file_name = "wolf_calibration_{}.conf".format(
                job.doc.calibration_iteration_number
            )
            with open(job.fn(output_name_control_file_name), "a") as myfile:
                defPotLine = "InitStep\t{zero}\n".format(zero=0)
                myfile.write(defPotLine)
                defPotLine = "Wolf\t{freq}\n".format(freq=True)
                myfile.write(defPotLine)
                defPotLine = "WolfKind\t{freq}\n".format(freq=job.sp.wolf_model)
                myfile.write(defPotLine)
                defPotLine = "WolfPotential\t{freq}\n".format(freq=job.sp.wolf_potential)
                myfile.write(defPotLine)   
                defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=b, val=nextAlpha)
                myfile.write(defAlphaLine)
                bestA.append(nextAlpha)
            return bestA
    except:
        print(repr(e))
        return False


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_sseq_completed_properly(job):
    """Check to see if the sseq simulation is finished (set temperature)."""
    try: 
        output_name_control_file_name = "single_state_eq"

        return gomc_sseq_completed_properly(
            job,
            output_name_control_file_name,
        )
    except:
        return False
# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_initial_guesses_generated(job):
#This will cause Ewald sims to wait for Wolf calibration to complete.
    #if(job.sp.wolf_potential == "Results"):
    #    return True
    try: 
        output_name_control_file_name = "wolf_calibration_{}".format(0)
        job_run_properly_bool = False
        output_log_file = job.doc.path_to_ew_results_repl_0_dir+"out_{}.dat".format(output_name_control_file_name)
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
    except:
        return False
# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file

@Project.label
@flow.with_job
def part_4b_job_gomc_wolf_sanity_completed_properly(job):
    """Check to see if the wolf_sanity simulation is finished (set temperature)."""
    try: 
        if (job.doc.winning_alpha):
            output_name_control_file_name = "wolf_sanity"
            return gomc_sim_completed_properly(
                job,
                output_name_control_file_name,
            )
        else:
            return False
    except:
        return False


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_wolf_sanity_individual_simulation_averages_completed(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    
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
@Project.pre(mosdef_input_written)
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
    dict_of_equilibrated_energies_stats = {}
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
                    #energies.append(float(line.split()[2]))
                    # Use Total_Elec to avoid underreporting error.
                    energies.append(float(line.split()[7]))

                except:
                    print(line)
                    print("An exception occurred") 
            if DensRegex.match(line):
                #print('\n'.join(line.split()[1] for line in f))
                try:
                    if (job.sp.solute in "solvent_box"):
                        if (job.doc.production_ensemble in ["NVT"]):
                            densities.append(float(line.split()[3]))
                        elif (job.doc.production_ensemble in ["NPT"]):
                            densities.append(float(line.split()[4]))      
                    elif(job.sp.solute in "ETOH"):    
                        if (job.doc.production_ensemble in ["NVT"]):
                            densities.append(float(line.split()[7]))
                        elif (job.doc.production_ensemble in ["NPT"]):
                            densities.append(float(line.split()[8]))      
                        elif (job.doc.production_ensemble in ["GEMC_NVT"]):
                            densities.append(float(line.split()[4])) 
                except:
                    print("An exception occurred") 
    steps_np = np.array(steps)
    energies_np = np.array(energies)
    densities_np = np.array(densities)

    nskip = 100

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

    print("Num correlated equilibrated energy samples",np.shape(energies_np)[0])
    dict_of_equilibrated_energies_stats[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_mean'] = [A_t_equil.mean()]
    dict_of_equilibrated_energies_stats[f'{job.sp.wolf_model}_{job.sp.wolf_potential}_std'] = [A_t_equil.std()]
    dfUC3 = pd.DataFrame.from_dict(dict_of_equilibrated_energies_stats)
    dfUC3.to_csv('wolf_sanity_equilibrated_energies_avg{}.csv'.format(job.id))


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
    try:
        if (job.isfile(job.doc.path_to_wolf_results_my_repl_dir+"/wolf_statistics.csv")):
            return True
    except:
        return False

@Project.label
@flow.with_job
def part_4b_wolf_sanity_histograms_created(job):
    try:
        df1 = pd.DataFrame()
        if (job.isfile(job.doc.path_to_wolf_results_repl_0_dir+"/wolf_sanity_all_energies.csv") and job.isfile(job.doc.path_to_wolf_results_repl_0_dir+"/PotentialEnergyDistribution_Ewald_vs_All.png")):
            return True
        else:
            return False
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
@Project.pre(lambda j: j.sp.wolf_potential == "Results")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_append_done)
@Project.pre(lambda *jobs: not any(part_4b_job_gomc_winning_alpha(j) and not\
     part_4b_wolf_sanity_individual_simulation_averages_completed(j) for j in jobs[0]._project))
@Project.post(part_4b_wolf_sanity_analysis_completed)
@flow.with_job
def part_4b_wolf_sanity_analysis(job):

    df1 = pd.DataFrame()
    df3 = pd.DataFrame()
    df5 = pd.DataFrame()

    # All different wolf models and ewald within a replicas
    jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "solute": job.sp.solute}))
    print(jobs)
    for other_job in jobs:
            if (other_job.sp.electrostatic_method == "Wolf" and \
                other_job.sp.wolf_potential == "Results"):
                continue
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
    ref_mean = df5["Results_Results"].mean()
    for method in listOfWolfMethods:
        print("Comparing statistical identicallness of Ewald and", method)
        welchs_output = scipy.stats.ttest_ind(df5["Results_Results"], df5[method], equal_var=False, nan_policy='omit')
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
    ref_mean = df1["Results_Results"].mean()
    for method in listOfWolfMethods:
        print("Comparing statistical identicallness of Ewald and", method)
        welchs_output = scipy.stats.ttest_ind(df1["Results_Results"], df1[method], equal_var=False, nan_policy='omit')
        statistics[method] = [df1[method].mean(), df1[method].std(),(df1[method].mean()-ref_mean)/ref_mean, welchs_output[0], welchs_output[1]]

    # Change the row indexes
    statistics.index = ['mean', 'std', 'relative_error', 't-statistic', 'p-value']   
    statistics = statistics.T.sort_values('p-value', ascending=False).T
    statistics.to_csv('wolf_statistics.csv', sep = ' ', )

    #job.doc.winningWolfModel = (statistics.columns[1]).split("_")[0]
    #job.doc.winningWolfPotential = (statistics.columns[1]).split("_")[1]
    print(statistics)


# ******************************************************
# ******************************************************
# data analysis - get the average data from each individual simulation (start)
# ******************************************************
# ******************************************************




# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
# For some reason this is failing on all but replica 0..
@Project.label
@flow.with_job
def part_4b_job_gomc_winning_alpha(job):
    try:
        if (job.doc.winning_alpha and not part_4b_job_gomc_wolf_parameters_appended(job)):
            print(job.fn(""), job.sp.wolf_model, job.sp.wolf_potential)
        return job.doc.winning_alpha
    except:
        return False


@Project.label
@flow.with_job
def part_4b_job_gomc_wolf_parameters_appended(job):

    if (job.sp.electrostatic_method == "Ewald"):
        return True

    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    import re
    regex = re.compile("(\w+?)_initial_state_(\w+?).conf")
    success = True
    atLeastOneMatchExists = False
    if (job.sp.electrostatic_method == "Ewald"):
        return True

    if(not job.isfile("wolf_sanity.conf")):
        return False
    regex = re.compile("wolf_sanity.conf")
    for root, dirs, files in os.walk(job.fn("")):
        for file in files:
            if regex.match(file):
                atLeastOneMatchExists = True
                with open(file, "r") as openedfile:
                    last_line = openedfile.readlines()[-1]
                if ("WolfAlpha" in last_line):
                    continue
                else:
                    success = success and False

    return success and atLeastOneMatchExists

# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(lambda j: j.sp.wolf_potential == "Results")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_all_replicates_best_alpha_found)
@Project.post(part_4b_append_done)
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
    import shutil

    cols = ["MODEL", "POT", "ALPHA"]
    bestAlphas = pd.read_csv("best_alpha_all_replicas.csv")
    print(bestAlphas)
    if (job.doc.equilibration_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]):  
        box_list = [0, 1]
    else:
        box_list = [0]
    for index, row in bestAlphas.iterrows():
        print(row)
        jobs = list(pr.find_jobs({"electrostatic_method": "Wolf", \
            "wolf_model": row['WOLF_KIND'], \
            "wolf_potential" : row['COUL_KIND']}))
        for ref_job in jobs:
            if (ref_job.sp.alpha == row['ALPHA']):
                print(ref_job.fn(""))

                #shutil.copyfile(job.doc.path_to_gomc_sseq_dir+wolf_sanity_control_file_name, wolf_sanity_control_file_name)
                #append_default_wolf_parameters_line(job,wolf_sanity_control_file_name)

                import re
                regex = re.compile("(\w+?)_initial_state_(\w+?).conf")
                for root, dirs, files in os.walk(ref_job.doc.path_to_gomc_sseq_dir):
                    for file in files:
                        if regex.match(file):
                            print(os.path.join(ref_job.doc.path_to_gomc_sseq_dir, file), "copied to", ref_job.fn(file))
                            shutil.copyfile(os.path.join(ref_job.doc.path_to_wolf_template_dir, file), ref_job.fn(file))
                            append_default_wolf_parameters_line(ref_job,file)
                
                regex = re.compile("wolf_sanity.conf")
                for root, dirs, files in os.walk(ref_job.doc.path_to_gomc_sseq_dir):
                    for file in files:
                        if regex.match(file):
                            print(os.path.join(ref_job.doc.path_to_gomc_sseq_dir, file), "copied to", ref_job.fn(file))
                            shutil.copyfile(os.path.join(ref_job.doc.path_to_wolf_template_dir, file), ref_job.fn(file))
                            append_default_wolf_parameters_line(ref_job,file)
                ref_job.doc.winning_alpha = True
            else:
                ref_job.doc.winning_alpha = False
    f = open("append_done.txt", "a")
    f.write("Now the file has more content!")
    f.close()
    job.doc.append_done = True
        




# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4b_job_gomc_equilb_design_ensemble_completed_properly(job):
    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    try:
        statesDone = True
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                filename_4b_iter = job.doc.gomc_equilb_design_ensemble_dict[
                    str(initial_state_i)
                ]["output_name_control_file_name"]

                statesDone = statesDone and gomc_sim_completed_properly(
                    job,
                    filename_4b_iter,
                )
            except:
                return False
        return statesDone
    except:
        return False

# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@flow.with_job
def part_4c_job_production_run_completed_properly(job):

    """Check to see if the gomc_equilb_design_ensemble simulation was completed properly (set temperature)."""
    try:
        statesDone = True
        for initial_state_i in list(job.doc.InitialState_list):
            try:
                filename_4b_iter = job.doc.gomc_production_run_ensemble_dict[
                    str(initial_state_i)
                ]["output_name_control_file_name"]

                statesDone = statesDone and gomc_sim_completed_properly(
                    job,
                    filename_4b_iter,
                )
            except:
                return False
        return statesDone
    except:
        return False

# ******************************************************
# ******************************************************
# check if GOMC and NAMD simulation are completed properly (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC anaylsis is completed properly (start)
# ******************************************************
# ******************************************************
# check if analysis is done for the individual replicates wrote the gomc files
@Project.pre(part_4b_job_gomc_equilb_design_ensemble_completed_properly)
@Project.label
@flow.with_job
def part_5a_preliminary_analysis_individual_simulation_averages_completed(job):
    """Check that the individual simulation averages files are written ."""
    file_written_bool = False
    if (
        job.isfile(
            f"{preliminary_output_replicate_txt_file_name_box_0}"
        )
    ):
        file_written_bool = True

    return file_written_bool



# check if analysis is done for the individual replicates wrote the gomc files
@Project.pre(part_4c_job_production_run_completed_properly)
@Project.label
@flow.with_job
def part_5a_analysis_individual_simulation_averages_completed(job):
    """Check that the individual simulation averages files are written ."""
    file_written_bool = False
    if (
        job.isfile(
            f"{output_replicate_txt_file_name_box_0}"
        )
    ):
        file_written_bool = True

    return file_written_bool

# check if analysis for averages of all the replicates is completed
@Project.pre(part_5a_preliminary_analysis_individual_simulation_averages_completed)
@Project.label
def part_5b_preliminary_analysis_replica_averages_completed(*jobs):
    """Check that the simulation replicate average and std. dev. files are written."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{preliminary_output_avg_std_of_replicates_txt_file_name_box_0}"
            )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return all_file_written_bool_pass


# check if analysis for averages of all the replicates is completed
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.label
def part_5b_analysis_replica_averages_completed(*jobs):
    """Check that the simulation replicate average and std. dev. files are written."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{output_avg_std_of_replicates_txt_file_name_box_0}"
            )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return all_file_written_bool_pass

# check if analysis for averages of all the replicates is completed
@Project.pre(part_5a_preliminary_analysis_individual_simulation_averages_completed)
@Project.label
def part_5b_preliminary_analysis_replica_averages_completed(*jobs):
    """Check that the simulation replicate average and std. dev. files are written."""
    file_written_bool_list = []
    all_file_written_bool_pass = False
    for job in jobs:
        file_written_bool = False

        if (
            job.isfile(
                f"../../analysis/{output_avg_std_of_replicates_txt_file_name_box_0}"
            )
        ):
            file_written_bool = True

        file_written_bool_list.append(file_written_bool)

    if False not in file_written_bool_list:
        all_file_written_bool_pass = True

    return all_file_written_bool_pass
# ******************************************************
# ******************************************************
# check if GOMC anaylsis is completed properly (end)
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

    #if job.doc.solvent not in ["TIP4"]:
        #solvent.energy_minimize(forcefield=forcefield_dict[job.doc.solvent], steps=10 ** 5)
    if (job.doc.N_liquid_solute > 0):
        smiles_or_mol2_solute = get_molecule_path(smiles_or_mol2_name_to_value_dict[job.sp.solute][job.sp.forcefield])

        if job.sp.solute in ["He", "Ne", "Kr", "Ar", "Xe", "Rn"]:
            solute = mb.Compound(name=job.doc.solute)
        else:
            solute = mb.load(smiles_or_mol2_solute[1],
                            smiles=smiles_or_mol2_solute[0]
                            )
        solute.name = job.sp.solute

        solute_ff = get_ff_path(forcefield_residue_to_ff_filename_dict[job.sp.solute][job.sp.forcefield])
        solvent_ff = get_ff_path(forcefield_residue_to_ff_filename_dict[job.sp.solvent][job.sp.forcefield])

        # only put the FF molecules in the simulation in the dictionaly input into the Chamm object.
        minimal_forcefield_dict = {solute.name: solute_ff,
                                solvent.name: solvent_ff
                                }

        #solute.energy_minimize(forcefield=forcefield_dict[job.sp.solute], steps=10 ** 5)

        # for trappe, currently unused'
        if (job.sp.forcefield == "TRAPPE"):
            bead_to_atom_name_dict = { '_CH3':'C', '_CH2':'C',  'O':'O', 'H':'H'}
        else:
            bead_to_atom_name_dict = None

        residues_list = [solute.name, solvent.name]
        print("residues_list  = " +str(residues_list ))

        #if job.doc.solvent in ["TIP4", "TIP3"]:
        # Just to prove the point of why Waibel2018 DSP is failing.
        #gomc_fix_bonds_angles_residues_list = [solvent.name]

        gomc_fix_bonds_angles_residues_list = [solvent.name, solute.name]
        #else:
        #    gomc_fix_bonds_angles_residues_list  = None
        print('Running: filling liquid box')
        box_0 = mb.fill_box(compound=[solute, solvent],
                            n_compounds=[job.doc.N_liquid_solute, job.doc.N_liquid_solvent],
                            box=[u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                                u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                                u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("nm"),
                                ],
                            seed=mbuild_box_seed_no
                            )
        print('Completed: filling liquid box')

    else:
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

    print('Running: GOMC FF file, and the psf and pdb files')
    if job.doc.production_ensemble in ["NVT", "NPT"]:
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

        print('Running: gomc_charmm')
        gomc_charmm = mf_charmm.Charmm(
            box_0,
            mosdef_structure_box_0_name_str,
            structure_box_1=None,
            filename_box_1=None,
            ff_filename=  gomc_ff_filename_str,
            forcefield_selection=minimal_forcefield_dict,
            residues=residues_list,
            bead_to_atom_name_dict=bead_to_atom_name_dict,
            gomc_fix_bonds_angles=gomc_fix_bonds_angles_residues_list,
        )

    else:
        raise ValueError("ERROR: The GCMC and GEMC ensembles are not supported in this script.")

    if write_files == True:
        gomc_charmm.write_inp()

        namd_charmm.write_inp()

        namd_charmm.write_psf()

        namd_charmm.write_pdb()

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
    if (job.sp.wolf_model == "Results"):
        [namd_charmm_object_with_files, gomc_charmm_object_with_files] = build_charmm(job, write_files=True)

    namd_restart_pdb_psf_file_name_str = mosdef_structure_box_0_name_str
    restart_control_file_name_str = namd_equilb_NPT_control_file_name_str

    # Get prefix to namd runs.
    namd_output_prefix = ""
    jobs = list(pr.find_jobs({"replica_number_int": 0, "electrostatic_method": "Ewald", "wolf_model": "Results"}))

    for ref_job in jobs:
        namd_output_prefix = ref_job.fn("")



    # All wolf jobs will just copy the output of namd and gomc precal eq from here to local dir
    job.doc.namd_output_dir = namd_output_prefix

    Coordinates_box_0 = namd_output_prefix+"{}.pdb".format(
        mosdef_structure_box_0_name_str
    )
    Structure_box_0 = namd_output_prefix+"{}.psf".format(
        mosdef_structure_box_0_name_str
    )
    Coordinates_box_1 = namd_output_prefix+"{}.pdb".format(
        mosdef_structure_box_1_name_str
    )
    Structure_box_1 = namd_output_prefix+"{}.psf".format(
        mosdef_structure_box_1_name_str
    )
    binCoordinates_box_0 = namd_output_prefix+"{}.restart.coor".format(
        restart_control_file_name_str
    )
    extendedSystem_box_0 = namd_output_prefix+"{}.restart.xsc".format(
        restart_control_file_name_str
    )
    binCoordinates_box_1 = namd_output_prefix+"{}.coor".format(
        mosdef_structure_box_1_name_str
    )
    extendedSystem_box_1 = namd_output_prefix+"{}.xsc".format(
        mosdef_structure_box_1_name_str
    )

    # Path to namd output
    gomc_eq_control_file_name_str = "single_state_eq"
    job.doc.path_to_namd_console =  namd_output_prefix+f"out_{namd_equilb_NPT_control_file_name_str}.dat"

    job.doc.path_to_ref_pdb =  Coordinates_box_0
    job.doc.path_to_ref_psf =  Structure_box_0
    job.doc.path_to_ref_binCoordinates =  binCoordinates_box_0
    job.doc.path_to_ref_extendedSystem =  extendedSystem_box_0

    if (job.doc.equilibration_ensemble in ["GCMC"]):  
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


    Single_state_gomc_eq_Coordinates_box_0 = namd_output_prefix+"{}_BOX_0_restart.pdb".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Structure_box_0 = namd_output_prefix+"{}_BOX_0_restart.psf".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_binCoordinates_box_0 = namd_output_prefix+"{}_BOX_0_restart.coor".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_extendedSystem_box_0 = namd_output_prefix+"{}_BOX_0_restart.xsc".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Coordinates_box_1 = namd_output_prefix+"{}_BOX_1_restart.pdb".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_Structure_box_1 = namd_output_prefix+"{}_BOX_1_restart.psf".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_binCoordinates_box_1 = namd_output_prefix+"{}_BOX_1_restart.coor".format(
        Single_state_gomc_eq_control_file_name
    )
    Single_state_gomc_eq_extendedSystem_box_1 = namd_output_prefix+"{}_BOX_1_restart.xsc".format(
        Single_state_gomc_eq_control_file_name
    )


    jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "electrostatic_method": "Ewald", "wolf_model": "Results"}))
    for ref_job in jobs:
        #if (ref_job.isfile(f"{Coordinates_box_0}")):
        job.doc.path_to_gomc_sseq_dir =  ref_job.fn("")
        job.doc.path_to_gomc_eq_console =  ref_job.fn(f"out_{gomc_eq_control_file_name_str}.dat")

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

    jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "electrostatic_method": "Wolf", "wolf_model": "Results"}))
    for ref_job in jobs:
        #if (ref_job.isfile(f"{Coordinates_box_0}")):
        job.doc.path_to_wolf_template_dir =  ref_job.fn("")

    jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "electrostatic_method": "Ewald", "wolf_model": "Results", "wolf_potential" : "Results"}))
    for results_job in jobs:
        job.doc.path_to_ew_results_my_repl_dir =  results_job.fn("")

    jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "electrostatic_method": "Wolf", "wolf_model": "Results", "wolf_potential" : "Results"}))
    for results_job in jobs:
        job.doc.path_to_wolf_results_my_repl_dir =  results_job.fn("")

    jobs = list(pr.find_jobs({"replica_number_int": 0, "electrostatic_method": "Ewald", "wolf_model": "Results", "wolf_potential" : "Results"}))
    for results_job in jobs:
        job.doc.path_to_ew_results_repl_0_dir =  results_job.fn("")

    jobs = list(pr.find_jobs({"replica_number_int": 0, "electrostatic_method": "Wolf", "wolf_model": "Results", "wolf_potential" : "Results"}))
    for results_job in jobs:
        job.doc.path_to_wolf_results_repl_0_dir =  results_job.fn("")


    FreeEnergyCalc = [True, int(gomc_free_energy_output_data_every_X_steps)]
    # This has to be off during calibration
    NoFreeEnergyCalc = [False, int(gomc_free_energy_output_data_every_X_steps)]

    MoleculeType = [job.sp.solute, 1]

    use_ElectroStatics = True

    if (job.sp.forcefield in ["OPLS"]):
        VDWGeometricSigma = True
    else:
        VDWGeometricSigma = False

    # max number of equilibrium selected runs
    calibration_iteration_number_max_number = int(10)

    # set the initial iteration number of the calibration simulation
    job.doc.equilb_design_ensemble_dict = {}
    job.doc.calibration_iteration_number = int(0)
    job.doc.calibration_iteration_number_max_number = (
        calibration_iteration_number_max_number
    )
    job.doc.calibration_iteration_number_under_limit = True
    job.doc.calibration_converged = False    
    job.doc.check_convergence = True

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

    box_lengths_ang = [u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       u.unyt_quantity(job.doc.liq_box_lengths_ang, 'angstrom').to_value("angstrom"),
                       ]

    seed_no = job.doc.replica_number_int

    if job.doc.equilibration_ensemble in ["NVT"]:
        namd_template_path_str = os.path.join(project_directory_path, "templates/NAMD_NVT_conf_template.conf")
    elif job.doc.equilibration_ensemble in ["NPT"]:
        namd_template_path_str = os.path.join(project_directory_path, "templates/NAMD_conf_template.conf")

    if job.doc.solvent in ["TIP3", "SPC", "SPCE", "MSPCE"]:
        namd_uses_water = True
        namd_water_model = 'tip3'
    elif job.doc.solvent in ["TIP4"]:
        namd_uses_water = True
        namd_water_model = 'tip4'
    else:
        namd_uses_water = False
        namd_water_model= None

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

    job.doc.winning_alpha = job.sp.electrostatic_method == "Ewald" and job.sp.wolf_model == "Results"

    if (job.sp.wolf_model != "Results"):
        return
    # generate the namd file
    # NOTE: the production and melt temps are converted to intergers so they can be ramped down
    # from hot to cool to equilibrate the system.
    generate_namd_equilb_control_file(template_path_filename=namd_template_path_str,
                                      namd_path_conf_filename=namd_equilb_NPT_control_file_name_str,
                                      namd_path_file_output_names=namd_equilb_NPT_control_file_name_str,
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
                                      box_lengths=box_lengths_ang,
                                      )

    print("#**********************")
    print("Completed: namd_equilb_NPT GOMC control file writing")
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
    output_true_list_input = [
        True,
        int(gomc_output_data_every_X_steps),
    ]
    output_false_list_input = [
        False,
        int(gomc_output_data_every_X_steps),
    ]
    # namd and gomc integrators are off by ~2%, which may cause drift
    # if you calibrate wolf using the restart files from namd
    # therefore a single state npt equilibration is performed
    # in gomc before wolf calibration.
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

        else:
            raise ValueError(
                "Moleules MC move ratios not listed for this solvent and solute or ensemble "
                "in the GOMC control file writer."
            )
    print("#**********************")
    print("Started: equilb NPT NAMD -> NPT GOMC control file writing")
    print("#**********************")

    gomc_control.write_gomc_control_file(
        gomc_charmm_object_with_files,
        Single_state_gomc_eq_control_file_name,
        job.doc.equilibration_ensemble,
        gomc_steps_equilb_design_ensemble,
        510,
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
        Coordinates_box_1=None,
        Structure_box_1=None,
        binCoordinates_box_1=None,
        extendedSystem_box_1=None,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": production_pressure_bar,
            "Ewald": True,
            "RcutCoulomb_box_0" : 14,
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
            "RestartFreq": output_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": console_output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_false_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_true_list_input,
            "Potential": cutoff_style,
            "LRC": True,
            "RcutLow": 1.0,
            "CBMC_First": CBMC_First[-1],
            "CBMC_Nth": CBMC_Nth[-1],
            "CBMC_Ang": CBMC_Ang[-1],
            "CBMC_Dih": CBMC_Dih[-1],
            #"FreeEnergyCalc": NoFreeEnergyCalc,
            #"MoleculeType": MoleculeType,
            #"InitialState": initial_state_sims_i,
            #"LambdaVDW": list(job.doc.LambdaVDW_list),
            #"LambdaCoulomb":  list(job.doc.LambdaCoul_list) if useCoul else None,
        },
    )

    print("#**********************")
    print("Started: Wolf Calibration GOMC control file writing")
    print("#**********************")
    #for cal_run in range(job.doc.calibration_iteration_number_max_number):
    for cal_run in range(2):

        #
        output_name_control_file_calibration_name = "wolf_calibration_{}".format(
            cal_run
        )


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
            510,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(gomc_ff_filename_str),
            Restart=True,
            RestartCheckpoint=True,
            ExpertMode=False,
            Coordinates_box_0=job.doc.path_to_ref_pdb,
            Structure_box_0=job.doc.path_to_ref_psf,
            binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
            extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
            binVelocities_box_0=None,
            Coordinates_box_1=None,
            Structure_box_1=None,
            binCoordinates_box_1=None,
            extendedSystem_box_1=None,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": production_pressure_bar,
                "Ewald": cal_run == 0,
                "RcutCoulomb_box_0" : 14,
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
                "RestartFreq": output_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_true_list_input,
                "Potential": cutoff_style,
                "LRC": True,
                "RcutLow": 1.0,
                "CBMC_First": CBMC_First[-1],
                "CBMC_Nth": CBMC_Nth[-1],
                "CBMC_Ang": CBMC_Ang[-1],
                "CBMC_Dih": CBMC_Dih[-1],
            },
        )
        if(cal_run == 0):
            append_wolf_calibration_parameters(job, output_name_control_file_calibration_name+".conf", cal_run)
        append_checkpoint_line(job, output_name_control_file_calibration_name, job.doc.path_to_sseq_checkpoint)

    ### Need to append Wolf Calibration lines since they aren't in MosDef

    print("#**********************")
    print("Completed: Wolf Calibration GOMC control file writing")
    print("#**********************")


    print("#**********************")
    print("Started:  Wolf Sanity GOMC control file writing")
    print("#**********************")
    wolf_sanity_control_file_name = "wolf_sanity"
    gomc_control.write_gomc_control_file(
        gomc_charmm_object_with_files,
        wolf_sanity_control_file_name,
        job.doc.production_ensemble,
        Wolf_Sanity_MC_steps,
        510,
        ff_psf_pdb_file_directory=None,
        check_input_files_exist=False,
        Parameters="{}.inp".format(gomc_ff_filename_str),
        Restart=True,
        RestartCheckpoint=True,
        ExpertMode=False,
        Coordinates_box_0=job.doc.path_to_ref_pdb,
        Structure_box_0=job.doc.path_to_ref_psf,
        binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
        extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
        binVelocities_box_0=None,
        Coordinates_box_1=None,
        Structure_box_1=None,
        binCoordinates_box_1=None,
        extendedSystem_box_1=None,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": production_pressure_bar,
            "Ewald": job.sp.electrostatic_method == "Ewald",
            "ElectroStatic": use_ElectroStatics,
            "RcutCoulomb_box_0" : 14,
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
            "PressureCalc": output_false_list_input,
            "RestartFreq": output_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": console_output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_false_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_true_list_input,
            "Potential": cutoff_style,
            "LRC": True,
            "RcutLow": 1.0,
            "CBMC_First": CBMC_First[-1],
            "CBMC_Nth": CBMC_Nth[-1],
            "CBMC_Ang": CBMC_Ang[-1],
            "CBMC_Dih": CBMC_Dih[-1],
            #"FreeEnergyCalc": NoFreeEnergyCalc,
            #"MoleculeType": MoleculeType,
            #"InitialState": initial_state_sims_i,
            #"LambdaVDW": list(job.doc.LambdaVDW_list),
            #"LambdaCoulomb":  list(job.doc.LambdaCoul_list) if useCoul else None,
        },
    )
    append_checkpoint_line(job, wolf_sanity_control_file_name, job.doc.path_to_sseq_checkpoint)

    print("#**********************")
    print("Finished: Wolf Sanity GOMC control file writing")
    print("#**********************")

    # ******************************************************
    # equilb selected_ensemble, if NVT -> NPT - GOMC control file writing  (start)
    # Note: the control files are written for the max number of gomc_equilb_design_ensemble runs
    # so the Charmm object only needs created 1 time.
    # ******************************************************
    print("#**********************")
    print("Started: equilb NPT or GEMC-NVT GOMC control file writing")
    print("#**********************")

    for initial_state_sims_i in list(job.doc.InitialState_list):
        output_name_control_file_name = "{}_initial_state_{}".format(
            gomc_equilb_design_ensemble_control_file_name_str, initial_state_sims_i
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
            Coordinates_box_0= Coordinates_box_0 if job.sp.electrostatic_method == "Ewald" else job.doc.path_to_ref_pdb,
            Structure_box_0=Structure_box_0 if job.sp.electrostatic_method == "Ewald" else job.doc.path_to_ref_psf,
            binCoordinates_box_0=job.doc.path_to_sseq_binCoordinates,
            extendedSystem_box_0=job.doc.path_to_sseq_extendedSystem,
            binVelocities_box_0=None,
            Coordinates_box_1=None,
            Structure_box_1=None,
            binCoordinates_box_1=None,
            extendedSystem_box_1=None,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": production_pressure_bar,
                "Ewald": job.sp.electrostatic_method == "Ewald",
                "ElectroStatic": use_ElectroStatics,
                "RcutCoulomb_box_0" : 14,
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
                "RestartFreq": output_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_true_list_input,
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
                "Ewald": job.sp.electrostatic_method == "Ewald",
                "RcutCoulomb_box_0" : 14,
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
                "RestartFreq": output_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": console_output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_true_list_input,
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
        append_checkpoint_line(job, output_name_control_file_name, "{}_restart.chk".format(restart_control_file_name_str))

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
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(mosdef_input_written)
@Project.post(part_4a_job_namd_equilb_NPT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.namd_node_ncpu,
        "ngpu": lambda job: 0,
        "memory": memory_needed,
        "walltime": walltime_namd_hr,
    }
)
@flow.with_job
@flow.cmd
def run_namd_equilb_NPT_gomc_command(job):
    """Run the namd_equilb_NPT simulation."""
    print("#**********************")
    print("# Started the run_namd_equilb_NPT_gomc_command.")
    print("#**********************")
    """Run the gomc_calibration_run_ensemble simulation."""
    
    control_file_name_str = namd_equilb_NPT_control_file_name_str
    
    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(namd_binary_path),
        str(job.doc.namd_equilb_NPT_gomc_binary_file),
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
# Note replica isnt required to be zero.  replicas start branching at this point
# Since GOMC is semi-reproducible (CPU binaries) with a seed.
@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(part_4a_job_namd_equilb_NPT_completed_properly)
@Project.pre(mosdef_input_written)
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

@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(part_4a_job_namd_equilb_NPT_completed_properly)
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.post(part_4b_job_gomc_sseq_average_obtained)
@Project.operation.with_directives(
    {
        "np": lambda job: job.doc.gomc_ncpu,
        "ngpu": lambda job: job.doc.gomc_ngpu,
        "memory": memory_needed,
        "walltime": walltime_gomc_equilbrium_hr,
    }
)
@flow.with_job
def average_sseq_electrostatic_energy(job):
    import re
    EnRegex = re.compile("ENER_0")
    
    blk_file = f'out_single_state_eq.dat'
    steps = []
    energies = []
    densities = []
    with open(blk_file, 'r', encoding='utf8') as f:
        for line in f:
            if EnRegex.match(line):
                try:
                    steps.append(float(line.split()[1]))
                    #energies.append(float(line.split()[2]))
                    # Use Total_Elec to avoid underreporting error.
                    energies.append(float(line.split()[7]))

                except:
                    print(line)
                    print("An exception occurred") 
    nskip = 20
    steps_np = np.array(steps)
    energies_np = np.array(energies)
    from pymbar import timeseries
    t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
    A_t_equil = energies_np[t0:]
    A_t_equil_steps = steps_np[t0:]

    dict_of_equilibrated_energies = {}
    dict_of_equilibrated_energies_stats = {}

    dict_of_equilibrated_energies["steps"] = A_t_equil_steps
    dict_of_equilibrated_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = A_t_equil

    dict_of_energy_mean = {"EWALD_MEAN":A_t_equil.mean()}
    print(dict_of_energy_mean)
    d = {'EWALD_MEAN': [A_t_equil.mean()]}
    df = pd.DataFrame(data=d)
    #df = pd.DataFrame.from_dict(dict_of_energy_mean)
    df.to_csv("ewald_average.csv", header=True)


def extract_electrostatic_energy(job, filename):
    import re
    EnRegex = re.compile("ENER_0")
    
    steps = []
    energies = []
    densities = []
    with open(filename, 'r', encoding='utf8') as f:
        for line in f:
            if EnRegex.match(line):
                try:
                    steps.append(float(line.split()[1]))
                    #energies.append(float(line.split()[2]))
                    # Use Total_Elec to avoid underreporting error.
                    energies.append(float(line.split()[7]))

                except:
                    print(line)
                    print("An exception occurred") 
    nskip = 20
    steps_np = np.array(steps)
    energies_np = np.array(energies)
    from pymbar import timeseries
    t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
    A_t_equil = energies_np[t0:]
    A_t_equil_steps = steps_np[t0:]


    """

    dict_of_equilibrated_energies = {}
    dict_of_equilibrated_energies_stats = {}

    dict_of_equilibrated_energies["steps"] = A_t_equil_steps
    dict_of_equilibrated_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = A_t_equil

    dict_of_energy_mean = {"EWALD_MEAN":A_t_equil.mean()}
    print(dict_of_energy_mean)
    d = {'EWALD_MEAN': [A_t_equil.mean()]}
    df = pd.DataFrame(data=d)
    #df = pd.DataFrame.from_dict(dict_of_energy_mean)
    df.to_csv("ewald_average.csv", header=True)
    """

    return A_t_equil.mean()
    

@Project.pre(part_1a_initial_data_input_to_json)
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.pre(part_4b_job_gomc_wolf_parameters_appended)
@Project.pre(part_4b_job_gomc_winning_alpha)
@Project.post(part_4b_job_gomc_wolf_sanity_completed_properly)
@Project.operation.with_directives(
    {
        "np": 8,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_gomc_production_hr,
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
        str(job.doc.gomc_production_ensemble_gomc_binary_file),
        str(job.doc.gomc_ncpu),
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
#@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.electrostatic_method == "Ewald")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_1b_under_equilb_design_ensemble_run_limit)
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.post(part_4b_job_gomc_calibration_completed_properly)
@Project.post(part_4b_initial_guesses_generated)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": 7,
    }
)
@flow.with_job
def generate_initial_guesses_for_calibration_run_gomc_command(job):
    """Run the gomc_calibration_run_ensemble simulation."""
    #control_file_name_str = "wolf_calibration"
    control_file_name_str = "wolf_calibration_{}".format(job.doc.calibration_iteration_number)

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(gomc_binary_path),
        str("GOMC_GPU_NVT"),
        str(1),
        str(control_file_name_str),
        str(control_file_name_str),
    )
    print('gomc generate_initial_guesses_for_calibration_run_gomc_command run_command = ' + str(run_command))
    import subprocess

    exec_run_command = subprocess.Popen(
        run_command, shell=True, stderr=subprocess.STDOUT
    )
    os.waitpid(exec_run_command.pid, 0)  # os.WSTOPPED) # 0)


# ******************************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
#@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_model != "Results")
#@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_1b_under_equilb_design_ensemble_run_limit)
@Project.pre(mosdef_input_written)
@Project.pre(part_4b_job_gomc_sseq_completed_properly)
@Project.pre(part_4b_initial_guesses_generated) 
@Project.post(part_4b_job_gomc_calibration_completed_properly)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": 7,
    }
)
@flow.with_job
def run_calibration_run_gomc_command(job):
    ew_ref = pd.read_csv (job.doc.path_to_ew_results_my_repl_dir+"ewald_average.csv", header=0)
    ew_mean = ew_ref['EWALD_MEAN'].iloc[0]
    print("Ew mean", ew_mean)
    from skopt import Optimizer
    import matplotlib.pyplot as plt
    from skopt.plots import plot_convergence, plot_evaluations, plot_objective, plot_regret
    opt = Optimizer([(0, 0.5)], "GP", acq_func="EI",
                acq_optimizer="sampling",
                initial_point_generator="lhs",
                n_initial_points=5)
    #for cal_run in range(job.doc.calibration_iteration_number_max_number):
    for cal_run in range(10):
        print("Calibration_iteration_number", cal_run)
        control_file_name_str = "wolf_calibration_{}".format(cal_run)
        """Run the gomc_calibration_run_ensemble simulation."""
        #control_file_name_str = "wolf_calibration"
        template_control_file_name_str = "wolf_calibration_{}".format(1)

        #job.doc.calibration_iteration_number = job.doc.calibration_iteration_number+1

        import shutil

        shutil.copyfile(job.doc.path_to_ew_results_repl_0_dir+template_control_file_name_str+".conf", control_file_name_str+".conf")
        shutil.copyfile(job.doc.path_to_ew_results_repl_0_dir+"in_gomc_FF.inp", "in_gomc_FF.inp")
        suggested = []
        if (cal_run == 0):
            suggested = part_4b_extract_best_initial_guess_from_ewald_calibration(job)
        else:
            suggested = opt.ask()
        
        with open(job.fn(control_file_name_str), "a") as myfile:
            defPotLine = "InitStep\t{zero}\n".format(zero=0)
            myfile.write(defPotLine)
            defPotLine = "Wolf\t{freq}\n".format(freq=True)
            myfile.write(defPotLine)
            defPotLine = "WolfKind\t{freq}\n".format(freq=job.sp.wolf_model)
            myfile.write(defPotLine)
            defPotLine = "WolfPotential\t{freq}\n".format(freq=job.sp.wolf_potential)
            myfile.write(defPotLine)  
            defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=0, val=suggested)
            myfile.write(defAlphaLine)

        print(f"Running simulation job id {job}")
        run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
            str(gomc_binary_path),
            str(job.doc.gomc_calibration_gomc_binary_file),
            str(1),
            str(control_file_name_str),
            str(control_file_name_str),
        )
        print('gomc gomc_calibration_run_ensemble run_command = ' + str(run_command))
        import subprocess

        exec_run_command = subprocess.Popen(
            run_command, shell=True, stderr=subprocess.STDOUT
        )
        os.waitpid(exec_run_command.pid, 0)  # os.WSTOPPED) # 0)
        # test if the simulation actualy finished before checkin and adding 1 to the equilb counter
        # test if the simulation actualy finished before checkin and adding 1 to the equilb counter

        if gomc_sim_completed_properly(
            job,
            control_file_name_str
        ):
            y = extract_electrostatic_energy(job, "out_{}.dat".format(control_file_name_str))
            print(suggested, y-ew_mean)
            res = opt.tell(suggested, y-ew_mean)
            #print("Incrementing calibration_iteration_number", job.doc.calibration_iteration_number)
            #job.doc.calibration_iteration_number = job.doc.calibration_iteration_number+1
            #print("Incrementing calibration_iteration_number", job.doc.calibration_iteration_number)
    print("x*=%.2f f(x*)=%.2f" % (res.x[0], res.fun))
    job.doc.best_alpha = res.x[0]
    job.doc.best_alpha_elec_mean = res.fun
    _ = plot_objective(res)
    plt.savefig('plot_objective.png', bbox_inches='tight')
    plt.show()
    _ = plot_convergence(res)
    plt.savefig('plot_convergence.png', bbox_inches='tight')
    plt.show()
    _ = plot_evaluations(res)
    plt.savefig('plot_evaluations.png', bbox_inches='tight')
    plt.show()
    _ = plot_regret(res)
    plt.savefig('plot_regret.png', bbox_inches='tight')
    plt.show()

# ******************************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************

@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Results")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
#@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(lambda *jobs: all(part_4b_job_gomc_calibration_completed_properly(j)
                               for j in jobs[0]._project))
@Project.pre(lambda *jobs: all(part_4b_job_gomc_sseq_average_obtained(j)
                               for j in jobs[0]._project))
#@Project.pre(lambda *jobs: all(part_4b_job_gomc_calibration_needs_to_be_checked(j)
#                               for j in jobs[0]._project))
@Project.post(part_4b_wrote_alpha_csv)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def write_replicate_alpha_csv(job):
    try:
        ew_ref = pd.read_csv (job.doc.path_to_ew_results_my_repl_dir+"ewald_average.csv", header=0)
        ew_mean = ew_ref['EWALD_MEAN'].iloc[0]
        print("Ew mean", ew_mean)
        jobs = list(pr.find_jobs({"replica_number_int": job.sp.replica_number_int, "electrostatic_method": job.sp.electrostatic_method, "solute": job.sp.solute}))
        try:
            NUMBOXES = 1
            for b in range (NUMBOXES):
                master = pd.DataFrame()
                columnLabels = ["WOLF_KIND","COUL_KIND","ALPHA","MEAN_REL_ERR"]
                cols = ["MODEL", "POT", "ALPHA"]
                colList = ["ALPHA","WAIBEL2018_DSF", "RAHBARI_DSF", "WAIBEL2019_DSF", "WAIBEL2018_DSP", "RAHBARI_DSP", "WAIBEL2019_DSP"]
                for ewald_job in jobs:
                    if (ewald_job.sp.wolf_model == "Results"):
                        continue

                    control_file_name = "wolf_calibration"

                    output_file = f'out_{control_file_name}.dat'

                    import re
                    EnRegex = re.compile("ENER_0")
                    
                    steps = []
                    energies = []
                    densities = []
                    with open(ewald_job.fn(output_file), 'r', encoding='utf8') as f:
                        for line in f:
                            if EnRegex.match(line):
                                try:
                                    steps.append(float(line.split()[1]))
                                    #energies.append(float(line.split()[2]))
                                    # Use Total_Elec to avoid underreporting error.
                                    energies.append(float(line.split()[7]))

                                except:
                                    print(line)
                                    print("An exception occurred") 
                    nskip = 1
                    steps_np = np.array(steps)
                    energies_np = np.array(energies)
                    from pymbar import timeseries
                    t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
                    A_t_equil = energies_np[t0:]
                    A_t_equil_steps = steps_np[t0:]

                    dict_of_equilibrated_energies = {}
                    dict_of_equilibrated_energies_stats = {}

                    dict_of_equilibrated_energies["steps"] = A_t_equil_steps
                    dict_of_equilibrated_energies[f'{job.sp.wolf_model}_{job.sp.wolf_potential}'] = A_t_equil


                    data = {'WOLF_KIND' : ewald_job.sp.wolf_model, 'COUL_KIND' : ewald_job.sp.wolf_potential,\
                         'ALPHA' : ewald_job.sp.alpha, 'MEAN_REL_ERR' : (A_t_equil.mean()-ew_mean)/ew_mean}
                    curr = pd.DataFrame(data, index=[0])  # the `index` argument is important 
                    master = master.append(curr, ignore_index=True)
                    #master = master.sort_values(['WOLF_KIND', 'COUL_KIND'],ascending = [True, False])
                master = master.sort_values(['WOLF_KIND', 'COUL_KIND', 'ALPHA'],ascending = [True, False, True])
                master.to_csv('full_calibration_replica_{}_BOX_{}.csv'.format(job.doc.replica_number_int, b), header=True, index=False, sep=',')

                curr = master.loc[master.groupby(['WOLF_KIND','COUL_KIND']).MEAN_REL_ERR.idxmin()].reset_index(drop=True)
                curr = curr.sort_values(['WOLF_KIND', 'COUL_KIND'],ascending = [True, False])
                print(curr)
                print(curr.columns)
                prev = pd.DataFrame()
                curr.to_csv('best_alphas_replica_{}_BOX_{}.csv'.format(job.doc.replica_number_int, b), header=True, index=False, sep=',')
        except:
            print(repr(e))
            return False
    except:
        print(repr(e))
        return False


@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Results")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(lambda *jobs: all(part_4b_job_gomc_calibration_completed_properly(j)
                               for j in jobs[0]._project))
@Project.pre(lambda *jobs: all(part_4b_wrote_alpha_csv(j)
                               for j in jobs[0]._project))
#@Project.pre(lambda *jobs: all(part_4b_job_gomc_calibration_needs_to_be_checked(j)
#                               for j in jobs[0]._project))
@Project.post(part_4b_all_replicates_best_alpha_found)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def get_minimum_alpha_across_replicas(job):
    try:
        jobs = list(pr.find_jobs({"electrostatic_method": job.sp.electrostatic_method, "solute": job.sp.solute,\
            "wolf_potential": job.sp.wolf_potential,"wolf_model": job.sp.wolf_model,}))
        try:
            master = pd.DataFrame()
            NUMBOXES = 1
            for b in range (NUMBOXES):

                cols = ["MODEL", "POT", "ALPHA"]
                colList = ["ALPHA","WAIBEL2018_DSF", "RAHBARI_DSF", "WAIBEL2019_DSF", "WAIBEL2018_DSP", "RAHBARI_DSP", "WAIBEL2019_DSP"]
                for ewald_job in jobs:
                    try:
                        curr = pd.DataFrame()
                        prev = pd.DataFrame()
                        if (ewald_job.isfile('full_calibration_replica_{}_BOX_{}.csv'.format(ewald_job.doc.replica_number_int, b))):
                            curr = pd.read_csv (ewald_job.fn('full_calibration_replica_{}_BOX_{}.csv'.format(ewald_job.doc.replica_number_int, b)), sep=",", header=0)
                        else:
                            print(ewald_job.fn('full_calibration_replica_{}_BOX_{}.csv'.format(ewald_job.doc.replica_number_int, b)), "DNE")
                            #return False
                        print(curr)
                        master = master.append(curr, ignore_index=False)
                        print("master")
                        print(master)
                    except:
                        print(ewald_job.fn(""), "exc")
                        # Make this job eligibile for the next cal run
            master = master.sort_values(['WOLF_KIND', 'COUL_KIND', 'ALPHA'],ascending = [True, False, True])
            print("master sorted")
            print(master)
            master.to_csv('full_calibration_all_replicas.csv', header=True, index=False, sep=',')
            curr = master.groupby(['WOLF_KIND','COUL_KIND','ALPHA']).mean()
            print("mean")

            print(curr)
            print("best")
            curr.reset_index(inplace=True)
            curr['MEAN_REL_ERR'] = curr['MEAN_REL_ERR'].abs()
            #bestA = curr.loc[curr.groupby(level=['WOLF_KIND','COUL_KIND']).MEAN_REL_ERR.idxmin()]
            bestA = curr.loc[curr.groupby(['WOLF_KIND','COUL_KIND']).MEAN_REL_ERR.idxmin()]
            #bestA = curr.min(level=['WOLF_KIND','COUL_KIND'])
            bestA.to_csv('best_alpha_all_replicas.csv', header=True, index=False, sep=',')
        except:
            print()
            print(repr(e))
            return False
    except:
        print(repr(e))
        return False

@Project.pre(lambda j: j.sp.electrostatic_method == "Wolf")
@Project.pre(lambda j: j.sp.wolf_potential == "Results")
@Project.pre(lambda j: j.sp.wolf_model == "Results")
#@Project.pre(lambda j: j.sp.solute == "solvent_box")
#@Project.pre(lambda j: j.sp.replica_number_int == 0)
@Project.pre(part_4b_append_done)
@Project.pre(lambda *jobs: all(part_4b_wolf_sanity_analysis_completed(j)
                               for j in jobs[0]._project))

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
    df_equilibrated_all = pd.DataFrame()

    jobs = list(pr.find_jobs({"electrostatic_method": job.sp.electrostatic_method, "solute": job.sp.solute, "wolf_model": "Results", "wolf_potential": "Results"}))
    print(jobs)
    try:
        for ewald_job in jobs:
            print(ewald_job.sp.wolf_model,ewald_job.sp.wolf_potential,ewald_job.sp.replica_number_int)
            if (ewald_job.isfile("wolf_sanity_equilibrated_energies.csv")):
                df1 = pd.read_csv (ewald_job.fn('wolf_sanity_equilibrated_energies.csv'), sep=',', header=0, na_values='NaN', index_col=0)
                df_equilibrated_all = df_equilibrated_all.append(df1, ignore_index=True)
            else:
                print(ewald_job.fn("wolf_sanity_equilibrated_energies.csv"), "DNE\n")
    except:
        return False
    
    try:
        for ewald_job in jobs:
            if (ewald_job.isfile("wolf_statistics_equilibrated.csv")):
                df2 = pd.read_csv (ewald_job.fn('wolf_statistics_equilibrated.csv'), sep=' ', header=0, na_values='NaN', index_col=0)
    except:
        return False

    print(df1)
    print(df2)
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.stats as st

    numBins = 100
    nskip = 100
    ref_ewald = df_equilibrated_all["Results_Results"]
    """
    from pymbar import timeseries
    t0, g, Neff_max = timeseries.detectEquilibration(ref_ewald, nskip=nskip) # compute indices of uncorrelated timeseries
    A_t_equil_ewald = ref_ewald[t0:].dropna()
    A_t_equil_steps_ewald = ref_ewald[t0:].dropna()
    """
    A_t_equil_ewald = ref_ewald.dropna()
    print(A_t_equil_ewald)

    #Col_Dict = {"GROSS_DSF": "Waibel2018a (DSF)", "VLUGT_DSF": 'Rahbari (DSF)', "VLUGTWINTRACUTOFF_DSF": 'Waibel2018b (DSF)',
    #"GROSS_DSP": 'Waibel2018a (DSP)', "VLUGT_DSP": 'Rahbari (DSP)', "VLUGTWINTRACUTOFF_DSP": 'Waibel2018b (DSP)'}
    Col_Dict = {"WAIBEL2018_DSF": "Waibel2018a (DSF)", "RAHBARI_DSF": 'Rahbari (DSF)', "WAIBEL2019_DSF": 'Waibel2018b (DSF)',
    "WAIBEL2018_DSP": 'Waibel2018a (DSP)', "RAHBARI_DSP": 'Rahbari (DSP)', "WAIBEL2019_DSP": 'Waibel2018b (DSP)'}

    colList = df_equilibrated_all.columns.tolist()
    colList.remove("Results_Results")
    colList.remove("steps")

    colList = ["WAIBEL2018_DSF", "RAHBARI_DSF", "WAIBEL2019_DSF", "WAIBEL2018_DSP", "RAHBARI_DSP", "WAIBEL2019_DSP"]

    ewald_mean = np.mean(A_t_equil_ewald)
    wolf_mean = 0.0

    for col, col_i in zip(colList, range(0, len(colList))):

        wolf = df_equilibrated_all[col]
        """
        t0, g, Neff_max = timeseries.detectEquilibration(wolf, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil_wolf = wolf[t0:].dropna()
        A_t_equil_steps_wolf = wolf[t0:].dropna()
        """
        A_t_equil_wolf = wolf.dropna()
        print(A_t_equil_wolf)

        ref_min = min(A_t_equil_ewald)
        ref_max = max(A_t_equil_ewald)

        wolf_min = min(A_t_equil_wolf)
        wolf_max = max(A_t_equil_wolf)

        print("ref_min",ref_min)
        print("ref_max",ref_max)

        print("wolf_min",wolf_min)
        print("wolf_max",wolf_max)

        xmin = min(float(ref_min), float(wolf_min))
        xmax = max(float(ref_max), float(wolf_max))

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
        plt.ylabel('Probability', fontsize=22)
        plt.xlabel('Total Electrostatic (K)', fontsize=22)
        plt.legend()
        plt.savefig("PotentialEnergyDistribution_Ewald_vs_{}".format(col), dpi=300, bbox_inches='tight')
        plt.figure().clear()   

    from matplotlib.figure import figaspect
    w, h = figaspect(9/16)
    figSP, axs = plt.subplots(2, 3, figsize=(w,h))
    counter = 0
    
    for col, col_i in zip(colList, range(0, len(colList))):
        print(colList)
        wolf = df_equilibrated_all[col]
        """
        t0, g, Neff_max = timeseries.detectEquilibration(wolf, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil_wolf = wolf[t0:].dropna()
        A_t_equil_steps_wolf = wolf[t0:].dropna()
        """
        A_t_equil_wolf = wolf.dropna()
        print(A_t_equil_wolf)

        wolf_mean = np.mean(A_t_equil_wolf)
        
        ref_min = min(A_t_equil_ewald)
        ref_max = max(A_t_equil_ewald)

        wolf_min = min(A_t_equil_wolf)
        wolf_max = max(A_t_equil_wolf)

        print("ref_min",ref_min)
        print("ref_max",ref_max)

        print("wolf_min",wolf_min)
        print("wolf_max",wolf_max)

        xmin = min(float(ref_min), float(wolf_min))
        xmax = max(float(ref_max), float(wolf_max))

        binWidth =  (xmax - xmin)/float(numBins)
        binList = np.arange(xmin, xmax+binWidth, binWidth)
        # estimate the line with probability density function (PDF)
        kde1 = st.gaussian_kde(A_t_equil_ewald).pdf(binList)
        kde2 = st.gaussian_kde(A_t_equil_wolf).pdf(binList)
        #plt.hist(wolf, density=True, bins=binList, alpha=1, label=col)  # density=False would make counts
        axs[counter // 3, counter % 3].plot(binList, kde1, color="black", linewidth=2, label="Ewald")
        axs[counter // 3, counter % 3].plot(binList, kde2, color="red", linewidth=2, label=Col_Dict[col])
        axs[counter // 3, counter % 3].title.set_text(Col_Dict[col])
        axs[counter // 3, counter % 3].set_xlabel('Total Electrostatic (K)', labelpad=20)
        axs[counter // 3, counter % 3].set_ylabel('Probability')
        from matplotlib.offsetbox import AnchoredText
        anchored_text = AnchoredText("{}%".format(round(((ewald_mean-wolf_mean)/ewald_mean)*100,3)), loc=2)
        axs[counter // 3, counter % 3].add_artist(anchored_text)
        counter = counter + 1

    #figSP.legend()
    #figSP.ylabel('Probability', fontsize=22)
    #figSP.xlabel('Total Electrostatic (K)', fontsize=22)
    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color='red', label='Wolf')
    black_patch = mpatches.Patch(color='black', label='Ewald')
    error = mpatches.Rectangle((1,1), 1,1, color='black', label='Relative Error of Mean', fill = False)
    figSP.legend(handles=[red_patch, black_patch, error], ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.05),)
    figSP.tight_layout(pad=1.5)
    figSP.savefig("PotentialEnergyDistribution_Ewald_vs_All", dpi=300, bbox_inches='tight')
  

    import scipy
    from scipy.stats import ttest_ind
    from scipy.spatial.distance import jensenshannon
    statistics = pd.DataFrame()
    listOfWolfMethods = list(df_equilibrated_all.columns.values.tolist())
    print(listOfWolfMethods)
    listOfWolfMethods.remove("steps")
    print(listOfWolfMethods)
    ref_mean = df_equilibrated_all["Results_Results"].mean()
    for method in listOfWolfMethods:
        print("Comparing statistical identicallness of Ewald and", method)
        welchs_output = scipy.stats.ttest_ind(df_equilibrated_all["Results_Results"], df_equilibrated_all[method], equal_var=False, nan_policy='omit')
        statistics[method] = [df_equilibrated_all[method].mean(), df_equilibrated_all[method].std(),(df_equilibrated_all[method].mean()-ref_mean)/ref_mean, welchs_output[0], welchs_output[1]]

    # Change the row indexes
    statistics.index = ['mean', 'std', 'relative_error', 't-statistic', 'p-value']   
    statistics = statistics.T.sort_values('p-value', ascending=False).T
    statistics.to_csv('wolf_statistics_all_replicates.csv', sep = ' ', )

    job.doc.winningWolfModel = (statistics.columns[1]).split("_")[0]
    job.doc.winningWolfPotential = (statistics.columns[1]).split("_")[1]
    print(statistics)

for initial_state_j in range(0, number_of_lambda_spacing_including_zero_int):
    
    @Project.pre(part_4a_job_namd_equilb_NPT_completed_properly)
    @Project.pre(part_4b_job_gomc_sseq_completed_properly)
    @Project.pre(part_4b_job_gomc_wolf_parameters_appended) 
    @Project.pre(part_4b_wolf_sanity_histograms_created)  
    @Project.pre(part_4b_wolf_sanity_analysis_completed)  
    #@Project.pre(part_4b_is_winning_wolf_model_or_ewald)
    #@Project.post(part_3b_output_gomc_equilb_design_ensemble_started)
    @Project.post(part_4b_job_gomc_equilb_design_ensemble_completed_properly)
    @Project.operation.with_directives(
        {
            "np": lambda job: job.doc.gomc_ncpu,
            "ngpu": lambda job: job.doc.gomc_ngpu,
            "memory": memory_needed,
            "walltime": walltime_gomc_equilbrium_hr,
        },
        name = f"gomc_equilb_design_ensemble_initial_state_{initial_state_j}"
    )
    @flow.with_job
    @flow.cmd
    def run_equilb_run_gomc_command(job, *, initial_state_j=initial_state_j):
        """Run the gomc_equilb_run_ensemble simulation."""
        control_file_name_str = job.doc.gomc_equilb_design_ensemble_dict[
            str(initial_state_j)
        ]["output_name_control_file_name"]

        print(f"Running simulation job id {job}")
        run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
            str(gomc_binary_path),
            str(job.doc.gomc_equilb_design_ensemble_gomc_binary_file),
            str(job.doc.gomc_ncpu),
            str(control_file_name_str),
            str(control_file_name_str),
        )

        print('gomc equilbrium_run run_command = ' + str(run_command))

        return run_command
# *****************************************
# ******************************************************
# equilb NPT - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
for initial_state_i in range(0, number_of_lambda_spacing_including_zero_int):
    @Project.pre(part_4b_job_gomc_equilb_design_ensemble_completed_properly)
    @Project.pre(part_4b_job_gomc_wolf_parameters_appended)
    #@Project.pre(part_5a_preliminary_analysis_individual_simulation_averages_completed)

    @Project.post(part_4c_job_production_run_completed_properly)
    @Project.operation.with_directives(
        {
            "np": lambda job: job.doc.gomc_ncpu,
            "ngpu": lambda job: job.doc.gomc_ngpu,
            "memory": memory_needed,
            "walltime": walltime_gomc_production_hr,
        },
        name = f"gomc_production_ensemble_initial_state_{initial_state_i}"
    )
    @flow.with_job
    @flow.cmd
    def run_production_run_gomc_command(job, *, initial_state_i=initial_state_i):
        """Run the gomc_production_ensemble simulation."""

        control_file_name_str = job.doc.gomc_production_run_ensemble_dict[
            str(initial_state_i)
        ]["output_name_control_file_name"]

        print(f"Running simulation job id {job}")
        run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
            str(gomc_binary_path),
            str(job.doc.gomc_production_ensemble_gomc_binary_file),
            str(job.doc.gomc_ncpu),
            str(control_file_name_str),
            str(control_file_name_str),
        )

        print('gomc production run_command = ' + str(run_command))

        return run_command

# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (end)
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
@Project.pre(part_4b_job_gomc_equilb_design_ensemble_completed_properly)
@Project.post(part_5a_preliminary_analysis_individual_simulation_averages_completed)
@flow.with_job
def part_5a_preliminary_analysis_individual_simulation_averages(job):

    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_solute_title = 'solute'  # column title title for temp
    output_column_dFE_MBAR_title = 'dFE_MBAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_MBAR_std_title = 'dFE_MBAR_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_TI_title = 'dFE_TI_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_TI_std_title = 'dFE_TI_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_BAR_title = 'dFE_BAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_BAR_std_title = 'dFE_BAR_std_kcal_per_mol'  # column title title for ds_MBAR


    files = []
    blk_files = []
    k_b = 1.9872036E-3  # kcal/mol/K
    temperature = job.sp.production_temperature_K
    k_b_T = temperature * k_b
    for initial_state_iter in range(0, number_of_lambda_spacing_including_zero_int):
        reading_filename_box_0_iter = f'Free_Energy_BOX_0_{gomc_equilb_design_ensemble_control_file_name_str}_' \
                                        f'initial_state_{initial_state_iter}.dat'
        files.append(reading_filename_box_0_iter)

    #All samples
    # for TI estimator
    dHdl = pd.concat([extract_dHdl(job.fn(f), T=temperature) for f in files])
    ti = TI().fit(dHdl)
    delta_ti, delta_std_ti = get_delta_TI_or_MBAR(ti, k_b_T)

    # for MBAR estimator
    u_nk = pd.concat([extract_u_nk(job.fn(f), T=temperature) for f in files])
    mbar = MBAR().fit(u_nk)
    delta_mbar, delta_std_mbar = get_delta_TI_or_MBAR(mbar, k_b_T)

    # for BAR estimator
    bar = BAR().fit(u_nk)
    delta_bar, delta_std_bar = get_delta_BAR(bar, k_b_T)


    cols = [output_column_temp_title,output_column_solute_title, output_column_dFE_MBAR_title, \
                output_column_dFE_MBAR_std_title, output_column_dFE_TI_title, output_column_dFE_TI_std_title, \
                output_column_dFE_BAR_title,output_column_dFE_BAR_std_title]

    allData = [job.sp.production_temperature_K,job.sp.solute,delta_mbar,delta_std_mbar,delta_ti,delta_std_ti,delta_bar,delta_std_bar]
    allDataDF = pd.DataFrame(np.array(allData).reshape(-1,len(allData)))
    allDataDF.columns = cols    #allDataDF.columns = cols
    print(allDataDF)
    allDataDF.to_csv(preliminary_output_replicate_txt_file_name_box_0, header=True, index=False, sep=',')
    """
    from pymbar import timeseries
    nskip = 100
    # Read the data for TI estimator and BAR or MBAR estimators.
    list_data_TI = []
    list_data_BAR = []
    for f in files:
        dHdl = extract_dHdl(f, T=temperature)
        u_nkr = extract_u_nk(f, T=temperature)
        #Detect uncorrelated samples using VDW+Coulomb term in derivative 
        # of energy time series (calculated for TI)
        srs = dHdl['VDW'] + dHdl['Coulomb'] 
        t0, g, Neff_max = timeseries.detectEquilibration(srs, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil = srs[t0:]
        list_data_TI.append(ss.statistical_inefficiency(dHdl, series=srs, conservative=False))
        list_data_BAR.append(ss.statistical_inefficiency(u_nkr, series=srs, conservative=False))

    # Correlated samples
    #for TI estimator
    print("Working on TI method ...")
    dHdl = pd.concat([ld for ld in list_data_TI])
    ti = TI().fit(dHdl)
    delta_ti, delta_std_ti = get_delta_TI_or_MBAR(ti, k_b_T)

    #for MBAR estimator
    print("Working on MBAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    mbar = MBAR().fit(u_nk)
    delta_mbar, delta_std_mbar = get_delta_TI_or_MBAR(mbar, k_b_T)

    #for BAR estimator
    print("Working on BAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    bar = BAR().fit(u_nk)
    delta_bar, delta_std_bar = get_delta_BAR(bar, k_b_T)

    deCorrData = [job.sp.production_temperature_K,job.sp.solute,delta_mbar,delta_std_mbar,delta_ti,delta_std_ti,delta_bar,delta_std_bar]
    deCorrDataDF = pd.DataFrame(np.array(deCorrData).reshape(-1,len(deCorrData)))
    deCorrDataDF.columns = cols
    print(deCorrDataDF)
    """

# ******************************************************
# ******************************************************
# data analysis - get the average data from each individual simulation (start)
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
@Project.pre(part_4c_job_production_run_completed_properly)
@Project.post(part_5a_analysis_individual_simulation_averages_completed)
@flow.with_job
def part_5a_analysis_individual_simulation_averages(job):

    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_solute_title = 'solute'  # column title title for temp
    output_column_dFE_MBAR_title = 'dFE_MBAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_MBAR_std_title = 'dFE_MBAR_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_TI_title = 'dFE_TI_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_TI_std_title = 'dFE_TI_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_BAR_title = 'dFE_BAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_BAR_std_title = 'dFE_BAR_std_kcal_per_mol'  # column title title for ds_MBAR


    files = []
    blk_files = []
    k_b = 1.9872036E-3  # kcal/mol/K
    temperature = job.sp.production_temperature_K
    k_b_T = temperature * k_b
    for initial_state_iter in range(0, number_of_lambda_spacing_including_zero_int):
        reading_filename_box_0_iter = f'Free_Energy_BOX_0_{gomc_production_control_file_name_str}_' \
                                        f'initial_state_{initial_state_iter}.dat'
        files.append(reading_filename_box_0_iter)

    #All samples
    # for TI estimator
    dHdl = pd.concat([extract_dHdl(job.fn(f), T=temperature) for f in files])
    ti = TI().fit(dHdl)
    delta_ti, delta_std_ti = get_delta_TI_or_MBAR(ti, k_b_T)

    # for MBAR estimator
    u_nk = pd.concat([extract_u_nk(job.fn(f), T=temperature) for f in files])
    mbar = MBAR().fit(u_nk)
    delta_mbar, delta_std_mbar = get_delta_TI_or_MBAR(mbar, k_b_T)

    # for BAR estimator
    bar = BAR().fit(u_nk)
    delta_bar, delta_std_bar = get_delta_BAR(bar, k_b_T)


    cols = [output_column_temp_title,output_column_solute_title, output_column_dFE_MBAR_title, \
                output_column_dFE_MBAR_std_title, output_column_dFE_TI_title, output_column_dFE_TI_std_title, \
                output_column_dFE_BAR_title,output_column_dFE_BAR_std_title]

    allData = [job.sp.production_temperature_K,job.sp.solute,delta_mbar,delta_std_mbar,delta_ti,delta_std_ti,delta_bar,delta_std_bar]
    allDataDF = pd.DataFrame(np.array(allData).reshape(-1,len(allData)))
    allDataDF.columns = cols    #allDataDF.columns = cols
    print(allDataDF)
    allDataDF.to_csv(output_replicate_txt_file_name_box_0, header=True, index=False, sep=',')
    """
    from pymbar import timeseries
    nskip = 100
    # Read the data for TI estimator and BAR or MBAR estimators.
    list_data_TI = []
    list_data_BAR = []
    for f in files:
        dHdl = extract_dHdl(f, T=temperature)
        u_nkr = extract_u_nk(f, T=temperature)
        #Detect uncorrelated samples using VDW+Coulomb term in derivative 
        # of energy time series (calculated for TI)
        srs = dHdl['VDW'] + dHdl['Coulomb'] 
        t0, g, Neff_max = timeseries.detectEquilibration(srs, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil = srs[t0:]
        list_data_TI.append(ss.statistical_inefficiency(dHdl, series=srs, conservative=False))
        list_data_BAR.append(ss.statistical_inefficiency(u_nkr, series=srs, conservative=False))

    # Correlated samples
    #for TI estimator
    print("Working on TI method ...")
    dHdl = pd.concat([ld for ld in list_data_TI])
    ti = TI().fit(dHdl)
    delta_ti, delta_std_ti = get_delta_TI_or_MBAR(ti, k_b_T)

    #for MBAR estimator
    print("Working on MBAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    mbar = MBAR().fit(u_nk)
    delta_mbar, delta_std_mbar = get_delta_TI_or_MBAR(mbar, k_b_T)

    #for BAR estimator
    print("Working on BAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    bar = BAR().fit(u_nk)
    delta_bar, delta_std_bar = get_delta_BAR(bar, k_b_T)

    deCorrData = [job.sp.production_temperature_K,job.sp.solute,delta_mbar,delta_std_mbar,delta_ti,delta_std_ti,delta_bar,delta_std_bar]
    deCorrDataDF = pd.DataFrame(np.array(deCorrData).reshape(-1,len(deCorrData)))
    deCorrDataDF.columns = cols
    print(deCorrDataDF)
    """


# ******************************************************
# ******************************************************
# data analysis - get the average data from each individual simulation (end)
# ******************************************************
# ******************************************************


# ******************************************************
# ******************************************************
# data analysis - get the average and std. dev. from/across all the replicates (start)
# ******************************************************
# ******************************************************

#@aggregator.groupby(key=statepoint_without_replica,
#                    sort_by="production_temperature_K",
#                    sort_ascending=True
#)
#@Project.operation.with_directives(
#     {
#         "np": 1,
#         "ngpu": 0,
#         "memory": memory_needed,
#         "walltime": walltime_gomc_analysis_hr,
#     }
#)

@Project.pre(part_4b_job_gomc_equilb_design_ensemble_completed_properly)
@Project.pre(part_5a_preliminary_analysis_individual_simulation_averages_completed)
@Project.post(part_5b_analysis_replica_averages_completed)
def part_5b_preliminary_analysis_replica_averages(*jobs):
    # ***************************************************
    #  create the required lists and file labels for the replicates (start)
    # ***************************************************
    # output and labels
    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_temp_std_title = 'temp_std_K'  # column title title for temp
    output_column_solute_title = 'solute'  # column title title for temp
    output_column_dFE_MBAR_title = 'dFE_MBAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_MBAR_std_title = 'dFE_MBAR_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_TI_title = 'dFE_TI_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_TI_std_title = 'dFE_TI_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_BAR_title = 'dFE_BAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_BAR_std_title = 'dFE_BAR_std_kcal_per_mol'  # column title title for ds_MBAR

    # get the list used in this function
    temp_repilcate_list = []
    solute_repilcate_list = []

    delta_MBAR_repilcate_box_0_list = []
    delta_TI_repilcate_box_0_list = []
    delta_BAR_repilcate_box_0_list = []


    output_txt_file_header = f"{output_column_temp_title: <30} " \
                             f"{output_column_temp_std_title: <30} " \
                             f"{output_column_solute_title: <30} "\
                             f"{output_column_dFE_MBAR_title: <30} "\
                             f"{output_column_dFE_MBAR_std_title: <30} "\
                             f"{output_column_dFE_TI_title: <3    0} "\
                             f"{output_column_dFE_TI_std_title: <30} "\
                             f"{output_column_dFE_BAR_title: <30} "\
                             f"{output_column_dFE_BAR_std_title: <30} "\
                             f"\n"


    write_file_path_and_name_box_0 = f'analysis/{preliminary_output_avg_std_of_replicates_txt_file_name_box_0}'
    if os.path.isfile(write_file_path_and_name_box_0):
        box_box_0_data_txt_file = open(write_file_path_and_name_box_0, "a")
    else:
        box_box_0_data_txt_file = open(write_file_path_and_name_box_0, "w")
        box_box_0_data_txt_file.write(output_txt_file_header)


    # ***************************************************
    #  create the required lists and file labels for the replicates (end)
    # ***************************************************

    for job in jobs:

        # *************************
        # drawing in data from single file and extracting specific rows from box 0 (start)
        # *************************
        reading_file_box_box_0 = job.fn(preliminary_output_replicate_txt_file_name_box_0)

        data_box_box_0 = pd.read_csv(reading_file_box_box_0, sep='\s+', header=0, na_values='NaN', index_col=False)
        data_box_box_0 = pd.DataFrame(data_box_box_0)

        temp_repilcate_list.append(data_box_box_0.loc[:, output_column_temp_title][0])
        solute_repilcate_list.append(data_box_box_0.loc[:, output_column_solute_title][0])

        delta_MBAR_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_MBAR_title][0])
        delta_TI_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_TI_title][0])
        delta_BAR_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_BAR_title][0])

        # *************************
        # drawing in data from single file and extracting specific rows from box 0 (end)
        # *************************


    # *************************
    # get the replica means and std.devs (start)
    # *************************
    temp_mean = np.mean(temp_repilcate_list)
    temp_std = np.std(temp_repilcate_list, ddof=1)

    solute_iter = solute_repilcate_list[0]

    delta_MBAR_mean_box_box_0 = np.mean(delta_MBAR_repilcate_box_0_list)
    delta_TI_mean_box_box_0 = np.mean(delta_TI_repilcate_box_0_list)
    delta_BAR_mean_box_box_0 = np.mean(delta_BAR_repilcate_box_0_list)

    delta_std_MBAR_mean_box_box_0 = np.std(delta_MBAR_repilcate_box_0_list, ddof=1)
    delta_std_TI_mean_box_box_0 = np.std(delta_TI_repilcate_box_0_list, ddof=1)
    delta_std_BAR_mean_box_box_0 = np.std(delta_BAR_repilcate_box_0_list, ddof=1)

    # *************************
    # get the replica means and std.devs (end)
    # *************************

    # ************************************
    # write the analysis data files for the liquid and vapor boxes (start)
    # ************************************

    box_box_0_data_txt_file.write(
        f"{temp_mean: <30} "
        f"{temp_std: <30} "
        f"{solute_iter: <30} "
        f"{delta_MBAR_mean_box_box_0: <30} "
        f"{delta_std_MBAR_mean_box_box_0: <30} "
        f"{delta_TI_mean_box_box_0: <30} "
        f"{delta_std_TI_mean_box_box_0: <30} "
        f"{delta_BAR_mean_box_box_0: <30} "
        f"{delta_std_BAR_mean_box_box_0: <30} "
        f" \n"
    )

    # ************************************
    # write the analysis data files for the liquid and vapor boxes (end)
    # ************************************


@Project.pre(part_4c_job_production_run_completed_properly)
@Project.pre(part_5a_analysis_individual_simulation_averages_completed)
@Project.post(part_5b_analysis_replica_averages_completed)
def part_5b_analysis_replica_averages(*jobs):
    # ***************************************************
    #  create the required lists and file labels for the replicates (start)
    # ***************************************************
    # output and labels
    output_column_temp_title = 'temp_K'  # column title title for temp
    output_column_temp_std_title = 'temp_std_K'  # column title title for temp
    output_column_solute_title = 'solute'  # column title title for temp
    output_column_dFE_MBAR_title = 'dFE_MBAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_MBAR_std_title = 'dFE_MBAR_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_TI_title = 'dFE_TI_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_TI_std_title = 'dFE_TI_std_kcal_per_mol'  # column title title for ds_MBAR
    output_column_dFE_BAR_title = 'dFE_BAR_kcal_per_mol'  # column title title for delta_MBAR
    output_column_dFE_BAR_std_title = 'dFE_BAR_std_kcal_per_mol'  # column title title for ds_MBAR

    # get the list used in this function
    temp_repilcate_list = []
    solute_repilcate_list = []

    delta_MBAR_repilcate_box_0_list = []
    delta_TI_repilcate_box_0_list = []
    delta_BAR_repilcate_box_0_list = []


    output_txt_file_header = f"{output_column_temp_title: <30} " \
                             f"{output_column_temp_std_title: <30} " \
                             f"{output_column_solute_title: <30} "\
                             f"{output_column_dFE_MBAR_title: <30} "\
                             f"{output_column_dFE_MBAR_std_title: <30} "\
                             f"{output_column_dFE_TI_title: <3    0} "\
                             f"{output_column_dFE_TI_std_title: <30} "\
                             f"{output_column_dFE_BAR_title: <30} "\
                             f"{output_column_dFE_BAR_std_title: <30} "\
                             f"\n"


    write_file_path_and_name_box_0 = f'analysis/{output_avg_std_of_replicates_txt_file_name_box_0}'
    if os.path.isfile(write_file_path_and_name_box_0):
        box_box_0_data_txt_file = open(write_file_path_and_name_box_0, "a")
    else:
        box_box_0_data_txt_file = open(write_file_path_and_name_box_0, "w")
        box_box_0_data_txt_file.write(output_txt_file_header)


    # ***************************************************
    #  create the required lists and file labels for the replicates (end)
    # ***************************************************

    for job in jobs:

        # *************************
        # drawing in data from single file and extracting specific rows from box 0 (start)
        # *************************
        reading_file_box_box_0 = job.fn(output_replicate_txt_file_name_box_0)

        data_box_box_0 = pd.read_csv(reading_file_box_box_0, sep='\s+', header=0, na_values='NaN', index_col=False)
        data_box_box_0 = pd.DataFrame(data_box_box_0)

        temp_repilcate_list.append(data_box_box_0.loc[:, output_column_temp_title][0])
        solute_repilcate_list.append(data_box_box_0.loc[:, output_column_solute_title][0])

        delta_MBAR_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_MBAR_title][0])
        delta_TI_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_TI_title][0])
        delta_BAR_repilcate_box_0_list.append(data_box_box_0.loc[:, output_column_dFE_BAR_title][0])

        # *************************
        # drawing in data from single file and extracting specific rows from box 0 (end)
        # *************************


    # *************************
    # get the replica means and std.devs (start)
    # *************************
    temp_mean = np.mean(temp_repilcate_list)
    temp_std = np.std(temp_repilcate_list, ddof=1)

    solute_iter = solute_repilcate_list[0]

    delta_MBAR_mean_box_box_0 = np.mean(delta_MBAR_repilcate_box_0_list)
    delta_TI_mean_box_box_0 = np.mean(delta_TI_repilcate_box_0_list)
    delta_BAR_mean_box_box_0 = np.mean(delta_BAR_repilcate_box_0_list)

    delta_std_MBAR_mean_box_box_0 = np.std(delta_MBAR_repilcate_box_0_list, ddof=1)
    delta_std_TI_mean_box_box_0 = np.std(delta_TI_repilcate_box_0_list, ddof=1)
    delta_std_BAR_mean_box_box_0 = np.std(delta_BAR_repilcate_box_0_list, ddof=1)

    # *************************
    # get the replica means and std.devs (end)
    # *************************

    # ************************************
    # write the analysis data files for the liquid and vapor boxes (start)
    # ************************************

    box_box_0_data_txt_file.write(
        f"{temp_mean: <30} "
        f"{temp_std: <30} "
        f"{solute_iter: <30} "
        f"{delta_MBAR_mean_box_box_0: <30} "
        f"{delta_std_MBAR_mean_box_box_0: <30} "
        f"{delta_TI_mean_box_box_0: <30} "
        f"{delta_std_TI_mean_box_box_0: <30} "
        f"{delta_BAR_mean_box_box_0: <30} "
        f"{delta_std_BAR_mean_box_box_0: <30} "
        f" \n"
    )

    # ************************************
    # write the analysis data files for the liquid and vapor boxes (end)
    # ************************************

# ******************************************************
# ******************************************************
# data analysis - get the average and std. dev. from/across all the replicates (end)
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
