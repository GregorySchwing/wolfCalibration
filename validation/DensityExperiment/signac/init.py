"""Initialize signac statepoints."""

import os
import numpy as np
import signac
import unyt as u

# *******************************************
# the main user varying state points (start)
# *******************************************

project=signac.init_project('density')

forcefield = {"SPCE" : "TRAPPE", "MTOH" : "OPLS"}
#solute = ["ETOH", "solvent_box"] # ["Ne", "Rn"]
solute = ["solvent_box"] # ["Ne", "Rn"]

#solute = ["Ne", "ETOH"] # ["Ne", "Rn"]
solvent = ["SPCE", "MTOH"] # ["Ne", "Rn"]
#solvent = ["SPC", "MSPCE"] # ["Ne", "Rn"]
electrostatic_method = ["Wolf", "Ewald"] # ["Ne", "Rn"]
replicas = [0, 1, 2, 3, 4]
#replicas = [0, 1]

g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
#densities_to_1 = np.arange (0.2, 0.3, 0.1)
densities_to_1 = np.arange (0.2, 1.1, 0.1)

densities_to_pt_7 = np.arange (0.2, 0.8, 0.1)

startingDensities = np.array([0.001, 0.01, 0.1])
densities_spce = np.concatenate((startingDensities, densities_to_1), axis=0)* g_per_cm3
print(densities_spce)

startingDensities_mtoh = np.array([0.00266, 0.01, 0.1])
densities_mtoh_to_7 = np.concatenate((startingDensities_mtoh, densities_to_pt_7), axis=0)
relevantLiquidDensities = np.array([0.692, 0.748, 0.7863])

densities_mtoh = np.concatenate((densities_mtoh_to_7, relevantLiquidDensities), axis=0)* g_per_cm3

#production_temperatures = [275, 295, 315, 335, 355, 375] * u.K # [275, 295, 315, 335, 355, 375] * u.K
#production_temperatures = [275] * u.K # [275, 295, 315, 335, 355, 375] * u.K
production_temperatures = [510] * u.K # [275, 295, 315, 335, 355, 375] * u.K

solvent_to_densities = {"SPCE": densities_spce, "MTOH": densities_mtoh}
# *******************************************
# the main user varying state points (end)
# *******************************************
#alpha_range = np.arange (0.0, 0.505, 0.005)
alpha_range = np.arange (0.0, 0.505, 0.005)
#alpha_range = np.linspace(0.0, 0.5, num=6)
print(alpha_range)
# *******************************************
# the main user varying state points (end)
# *******************************************



print("os.getcwd() = " +str(os.getcwd()))

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

wolfPotential = ["DSP","DSF"] # ["Ne", "Rn"]
wolfModel = ["RAHBARI","WAIBEL2018","WAIBEL2019"] # ["Ne", "Rn"]

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()

for replica_i in replicas:
    for solute_i in solute:
        for solvent_i in solvent:
            for ff_i in forcefield:
                for prod_temp_i in production_temperatures:
                    for density_i in solvent_to_densities[solvent_i]:
                        for e_method in electrostatic_method:
                            if (e_method == "Wolf"):
                                for wolfM in wolfModel:
                                    for wolfP in wolfPotential:
                                        for alpha in alpha_range:
                                            statepoint = {
                                                "replica_number_int": replica_i,
                                                "solvent": solvent_i,
                                                "solute": solute_i,
                                                "density" : np.round(density_i, 3),
                                                "forcefield" : forcefield[solvent_i],
                                                "wolf_model": wolfM,
                                                "wolf_potential": wolfP,
                                                "alpha": np.round(alpha, 3),
                                                "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                                "electrostatic_method": e_method,
                                            }
                                            total_statepoints.append(statepoint)
                                                            # The calibration statepoints
                                statepoint = {
                                                "replica_number_int": replica_i,
                                                "solute": solute_i,
                                                "solvent": solvent_i,
                                                "density" : np.round(density_i, 3),
                                                "forcefield" : forcefield[solvent_i],
                                                "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                                "electrostatic_method": "Wolf",
                                                "wolf_model": "Results",
                                                "wolf_potential": "Results",
                                            }
                                total_statepoints.append(statepoint)                     

                                # The calibration statepoints
                                statepoint = {
                                                "replica_number_int": replica_i,
                                                "solute": solute_i,
                                                "density" : np.round(density_i, 3),
                                                "solvent": solvent_i,
                                                "forcefield" : forcefield[solvent_i],
                                                "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
                                                "electrostatic_method": "Ewald",
                                                "wolf_model": "Results",
                                                "wolf_potential": "Results",
                                            }
                                total_statepoints.append(statepoint)
                        

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()


# *******************************************
# the other state points (start)
# *******************************************
