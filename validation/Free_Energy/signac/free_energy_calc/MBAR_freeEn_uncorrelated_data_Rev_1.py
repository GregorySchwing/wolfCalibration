from alchemlyb.parsing.gomc import  extract_dHdl,  extract_u_nk
from alchemlyb.estimators import MBAR, BAR, TI 
import alchemlyb.preprocessing.subsampling as ss
import pandas as pd
import numpy as np
import os

base_directory = os.path.dirname(os.path.realpath(__file__))
file_main_folder = 'workspace_renamed_solute_temps'


def get_delta(est, k_b_T):
    """ Return the change in free energy and standard deviation for TI and MBAR estimators.
    
    """
    delta = est.delta_f_.iloc[0, -1] * k_b_T
    d_delta = est.d_delta_f_.iloc[0, -1] * k_b_T
    return delta, d_delta


def get_delta2(est, k_b_T):
    """ Return the change in free energy and standard deviation for BAR estimator.
    
    """
    ee = 0.0

    for i in range(len(est.d_delta_f_) - 1):
        ee += est.d_delta_f_.values[i][i+1]**2
    
    delta = est.delta_f_.iloc[0, -1] * k_b_T
    d_delta = k_b_T * ee**0.5
    return delta, d_delta

##################################################

numFile = 17
temps_list = [275, 295, 315, 335, 355, 375] # temps_list = [275, 295, 315, 335, 355, 375] #

solvent_list = ["waterTIP4P"]
solute_list = ["Ne", "Rn"] #["Ne", "Rn"]
fname = "Free_Energy_BOX_0_gomc_production_run_initial_state_"
ext = ".dat"
replicate_list = 3

solvent_list = ["waterTIP4P"]

base_free_energy_file_name_str = "free_energy_data"

# loop through solvent list
for solvent_i in solvent_list:

    # loop through all solute_list
    for solute_i in solute_list:
        avg_output_data_csv = open(f'avg_std_free_energy_data_{solvent_i}_{solute_i}.csv', "w")
        avg_output_data_csv.write(
            "%25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, \n" % ("solvent",
                                                                          "solute",
                                                                          "temp_K",
                                                                          "avg_MBAR_kcal_per_mol",
                                                                          "std_dev_MBAR_kcal_per_mol",
                                                                          "avg_TI_kcal_per_mol",
                                                                          "std_dev_TI_kcal_per_mol",
                                                                          "avg_BAR_kcal_per_mol",
                                                                          "std_dev_BAR_kcal_per_mol",
                                                                          )
        )

        avg_output_data_txt = open(f'avg_std_free_energy_data_{solvent_i}_{solute_i}.txt', "w")
        avg_output_data_txt.write(
            "%25s %25s %25s %25s %25s %25s %25s %25s %25s \n" % ("solvent",
                                                                 "solute",
                                                                 "temp_K",
                                                                 "avg_MBAR_kcal_per_mol",
                                                                 "std_dev_MBAR_kcal_per_mol",
                                                                 "avg_TI_kcal_per_mol",
                                                                 "std_dev_TI_kcal_per_mol",
                                                                 "avg_BAR_kcal_per_mol",
                                                                 "std_dev_BAR_kcal_per_mol",
                                                                 )
        )

        replicate_output_data_csv = open(f'replicate_free_energy_data__{solvent_i}_{solute_i}.csv', "w")
        replicate_output_data_csv.write(
            "%25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, \n" % ("replicate_number",
                                                                                "solvent",
                                                                                "solute",
                                                                                "temp_K",
                                                                                "sum_MBAR_kcal_per_mol",
                                                                                "sum_ds_MBAR_kcal_per_mol",
                                                                                "sum_TI_kcal_per_mol",
                                                                                "sum_ds_TI_kcal_per_mol",
                                                                                "sum_BAR_kcal_per_mol",
                                                                                "sum_ds_BAR_kcal_per_mol",
                                                                                )
        )

        replicate_output_data_txt = open(f'replicate_free_energy_data_{solvent_i}_{solute_i}.txt', "w")
        replicate_output_data_txt.write(
            "%25s %25s %25s %25s %25s %25s %25s %25s %25s %25s\n" % ("replicate_number",
                                                                     "solvent",
                                                                     "solute",
                                                                     "temp_K",
                                                                     "sum_MBAR_kcal_per_mol",
                                                                     "sum_ds_MBAR_kcal_per_mol",
                                                                     "sum_TI_kcal_per_mol",
                                                                     "sum_ds_TI_kcal_per_mol",
                                                                     "sum_BAR_kcal_per_mol",
                                                                     "sum_ds_BAR_kcal_per_mol",
                                                                     )
        )



        for temp_i in temps_list:
            tis = []
            mbars = []
            bars = []
            temprature = temp_i  # temperature (K)
            k_b = 1.9872036E-3  # kcal/mol/K
            k_b_T = temprature * k_b

            for nr in range(replicate_list):
                # read the free energy files
                data_loc = f"{base_directory}/../{file_main_folder}/set_{nr}/{solute_i}/{temprature}K"
                files = []
                for i in range(numFile):
                    freeEn_file = fname + str(i) + ext
                    file_path = data_loc + "/" + freeEn_file
                    files.append(file_path)

                # for TI estimator
                dHdl = pd.concat([extract_dHdl(f, T=temprature) for f in files])
                ti = TI().fit(dHdl)
                sum_ti, sum_ds_ti = get_delta(ti, k_b_T)
                tis.append(sum_ti)

                # for MBAR estimator
                u_nk = pd.concat([extract_u_nk(f, T=temprature) for f in files])
                mbar = MBAR().fit(u_nk)
                sum_mbar, sum_ds_mbar = get_delta(mbar, k_b_T)
                mbars.append(sum_mbar)

                # for BAR estimator
                bar = BAR().fit(u_nk)
                sum_bar, sum_ds_bar = get_delta2(bar, k_b_T)
                bars.append(sum_bar)

                print(f"Replicate %d data in kcal/mol for the solvent = %s and solute = %s at %s K \n"
                      "sum_MBAR = %7.4f, sum_ds_MBAR = %7.4f, \n"
                      "  sum_TI = %7.4f,   sum_ds_TI = %7.4f, \n"
                      " sum_BAR = %7.4f,  sum_ds_BAR = %7.4f" % (
                    nr,
                    f"{solvent_i}",
                    f"{solute_i}",
                    f"{temp_i}",
                    sum_mbar,
                    sum_ds_mbar,
                    sum_ti,
                    sum_ds_ti,
                    sum_bar,
                    sum_ds_bar
                )
                      )

                replicate_output_data_csv.write(
                    "%25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s,\n" % (str(nr),
                                                                                       str(solvent_i),
                                                                                       str(solute_i),
                                                                                       str(temp_i),
                                                                                       str(sum_mbar),
                                                                                       str(sum_ds_mbar),
                                                                                       str(sum_ti),
                                                                                       str(sum_ds_ti),
                                                                                       str(sum_bar),
                                                                                       str(sum_ds_mbar),
                                                                                       )
                )

                replicate_output_data_txt.write(
                    "%25s %25s %25s %25s %25s %25s %25s %25s %25s %25s \n" % (str(nr),
                                                                              str(solvent_i),
                                                                              str(solute_i),
                                                                              str(temp_i),
                                                                              str(sum_mbar),
                                                                              str(sum_ds_mbar),
                                                                              str(sum_ti),
                                                                              str(sum_ds_ti),
                                                                              str(sum_bar),
                                                                              str(sum_ds_mbar),
                                                                              )
                )

            tis = np.array(tis)
            mbars = np.array(mbars)
            bars = np.array(bars)

            print("Average and Standard Deviations in kcal/mol for the solvent = %s and solute = %s at %s K \n"
                  "avg_MBAR = %7.4f, std_MBAR = %7.4f, \n"
                  "  avg_TI = %7.4f,   std_TI = %7.4f, \n"
                  " avg_BAR = %7.4f,  std_BAR = %7.4f" % (
                f"{solvent_i}",
                f"{solute_i}",
                f"{temp_i}",
                np.average(mbars),
                np.std(mbars),
                np.average(tis),
                np.std(tis),
                np.average(bars),
                np.std(bars)
            )
                  )

            avg_output_data_csv.write(
                "%25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, %25s, \n" % (str(solvent_i),
                                                                              str(solute_i),
                                                                              str(temp_i),
                                                                              str(np.average(mbars)),
                                                                              str(np.std(mbars)),
                                                                              str(np.average(tis)),
                                                                              str(np.std(tis)),
                                                                              str(np.average(bars)),
                                                                              str(np.std(bars)),
                                                                              )
            )

            avg_output_data_txt.write(
                "%25s %25s %25s %25s %25s %25s %25s %25s %25s \n" % (str(solvent_i),
                                                                     str(solute_i),
                                                                     str(temp_i),
                                                                     str(np.average(mbars)),
                                                                     str(np.std(mbars)),
                                                                     str(np.average(tis)),
                                                                     str(np.std(tis)),
                                                                     str(np.average(bars)),
                                                                     str(np.std(bars)),
                                                                     )
            )


        avg_output_data_csv.close()
        avg_output_data_txt.close()
        replicate_output_data_csv.close()
        replicate_output_data_csv.close()
