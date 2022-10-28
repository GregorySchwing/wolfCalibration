"""Initialize signac statepoints."""

import os
import numpy as np
import signac
import unyt as u

# *******************************************
# the main user varying state points (start)
# *******************************************

project=signac.init_project('density')

forcefield = ["TRAPPE"] # ["Ne", "Rn"]
solute = ["solvent_box"] # ["Ne", "Rn"]
#solute = ["Ne", "ETOH"] # ["Ne", "Rn"]
solvent = ["SPCE"] # ["Ne", "Rn"]
#solvent = ["SPC", "MSPCE"] # ["Ne", "Rn"]
electrostatic_method = ["Wolf", "Ewald"] # ["Ne", "Rn"]

#replicas = [0, 1, 2, 3, 4] # [0, 1, 2, 3, 4]
replicas = [0]# [0, 1, 2, 3, 4]
g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
#densities_to_1 = np.arange (0.2, 0.3, 0.1)
densities_to_1 = np.arange (0.2, 1.1, 0.1)

startingDensities = np.array([0.001, 0.01, 0.1])
densities = np.concatenate((startingDensities, densities_to_1), axis=0)* g_per_cm3
print(densities)
#production_temperatures = [275, 295, 315, 335, 355, 375] * u.K # [275, 295, 315, 335, 355, 375] * u.K
#production_temperatures = [275] * u.K # [275, 295, 315, 335, 355, 375] * u.K
production_temperatures = [298] * u.K # [275, 295, 315, 335, 355, 375] * u.K

 
# *******************************************
# the main user varying state points (end)
# *******************************************



print("os.getcwd() = " +str(os.getcwd()))

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

wolfPotential = ["DSP","DSF"] # ["Ne", "Rn"]
wolfModel = ["VLUGT","VLUGTWINTRACUTOFF","GROSS"] # ["Ne", "Rn"]

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()

for replica_i in replicas:
    for solvent_i in solvent:
        for solute_i in solute:
            for prod_temp_i in production_temperatures:
                for density_i in densities:
                    for e_method in electrostatic_method:
                        if (e_method == "Wolf"):
                            for wolfM in wolfModel:
                                for wolfP in wolfPotential:
                                    statepoint = {
                                        "replica_number_int": replica_i,
                                        "solvent": solvent_i,
                                        "solute": solute_i,
                                        "density" : density_i,
                                        "forcefield" : forcefield[0],
                                        "wolf_model": wolfM,
                                        "wolf_potential": wolfP,
                                        "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                        "electrostatic_method": e_method,
                                    }
                                    total_statepoints.append(statepoint)
                        else:
                            statepoint = {
                                    "replica_number_int": 0,
                                    "solvent": solvent_i,
                                    "solute": solute_i,
                                    "density" : density_i,
                                    "forcefield" : forcefield[0],
                                    "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                    "electrostatic_method": "Ewald",
                                    "wolf_model": "Ewald",
                                    "wolf_potential": "Ewald",
                            }
                            total_statepoints.append(statepoint)           
                    # The calibration statepoints
                    statepoint = {
                                    "replica_number_int": 0,
                                    "solvent": solvent_i,
                                    "solute": solute_i,
                                    "density" : density_i,
                                    "forcefield" : forcefield[0],
                                    "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                    "electrostatic_method": "Wolf",
                                    "wolf_model": "Calibrator",
                                    "wolf_potential": "Calibrator",
                                }
                    total_statepoints.append(statepoint)           
                    

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()


# *******************************************
# the other state points (start)
# *******************************************
