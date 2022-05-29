"""Initialize signac statepoints."""

import os
import numpy as np
import signac
import unyt as u

# *******************************************
# the main user varying state points (start)
# *******************************************

project=signac.init_project('nobel_gas_in_water')

#solute = ["Ne", "Rn"] # ["Ne", "Rn"]
#solute = ["ETOH"] # ["Ne", "Rn"]
solute = ["Ne"] # ["Ne", "Rn"]
#replicas = [0, 1, 2, 3, 4] # [0, 1, 2, 3, 4]
replicas = [0]# [0, 1, 2, 3, 4]
#production_temperatures = [275, 295, 315, 335, 355, 375] * u.K # [275, 295, 315, 335, 355, 375] * u.K
production_temperatures = [275] * u.K # [275, 295, 315, 335, 355, 375] * u.K


# *******************************************
# the main user varying state points (end)
# *******************************************



print("os.getcwd() = " +str(os.getcwd()))

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()

for replica_i in replicas:
    for solute_i in solute:
        for prod_temp_i in production_temperatures:
            statepoint = {
                "replica_number_int": replica_i,
                "solute": solute_i,
                "production_temperature_K": np.round(prod_temp_i.to_value("K"), 4),
            }

            total_statepoints.append(statepoint)

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()


# *******************************************
# the other state points (start)
# *******************************************
