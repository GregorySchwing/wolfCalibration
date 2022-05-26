import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import csv as csv
import pandas as pd
import itertools as it
from collections import defaultdict
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib.axis as axis
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import math
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import pylab as plot

#********************************
#Notes
#********************************

#Plots 2 item:
#1) a pairwise energy for indiviual peeling strands, and non-peeling strands
#2)
#********************************
#End Notes
#********************************
# only change the below reading_file
axis_Label_font_size = 24
legend_font_size = 18
axis_number_font_size = 20

PointSizes = 8
MarkeredgeWidth = 2
ConnectionLineSizes = 2
Calc_Linestyle = 'none'  #'-' =  solid, '--' = dashed, ':'= dotted
Exp_Linestyle = '-'

#**************
# T vs density plot ranges (start)
#*************
Ne_Free_energy_cal_per_mol_min = 0
Ne_Free_energy_cal_per_mol_max = 4

Rn_Free_energy_cal_per_mol_min = 0
Rn_Free_energy_cal_per_mol_max = 5

T_min = 270
T_max = 390

Free_energy_cal_per_mol_major = 1
Free_energy_cal_per_mol_minor = 0.2

T_ticks_major = 20
T_ticks_minor = 5

Inverse_T_scalar_10_to_Exp = 4
#**************
# P vs T plot ranges (end)
#*************

solute_name_list = ["Ne", "Rn"] # ["Ne", "Rn"]

for solute_i in solute_name_list:

    figure_name = f'free_energy_{solute_i}.pdf'

    Free_energy_calc_file_reading_name = f"analysis_avg_std_of_replicates_box_0.txt"
    Free_energy_calc_file_Vrabed_reading_name = f"Vrabec_et_al_Rn_Mick_et_al/" \
                                                f"avg_std_free_energy_data_waterTIP4P_Vrabec_et_al_Rn_Mick_et_al_" \
                                                f"{solute_i}.csv"
    Free_energy_exp_file_NIST_reading_name = f"NIST_website/avg_std_free_energy_data_NIST_website_{solute_i}.csv"


    Color_for_Plot_Data_calc = 'r'
    Color_for_Plot_Vrabec_calc =  'k'
    Color_for_Plot_NIST_exp =  'b'

    marker_gomc = 's'
    marker_Vrabec = 'o'

    #********************************
    #  File importing
    #********************************

    #data import Our data calcd
    Free_energy_calc_data = pd.read_csv(Free_energy_calc_file_reading_name,  sep='\s+',
                                        header=0,
                                        #index_col=0,
                                        #usecols=[2,3,4]
                                        )
    Free_energy_calc_data_df = pd.DataFrame(Free_energy_calc_data)
    Free_energy_calc_data_df = Free_energy_calc_data_df.loc[lambda df: df['solute'] == solute_i]
    print('Free_energy_calc_data_df')
    print(Free_energy_calc_data_df)

    Temp_calc_data = Free_energy_calc_data_df.loc[:, "temp_K"].values.tolist()
    print(Temp_calc_data)

    avg_free_energy_calc_free = \
        Free_energy_calc_data_df.loc[:, "dFE_MBAR_kcal_per_mol"].values.tolist()

    std_dev_free_energy_calc_free = \
        Free_energy_calc_data_df.loc[:, "dFE_MBAR_std_kcal_per_mol"].values.tolist()

    # data import Vrabec data calcd
    Free_energy_Vrabec_calc_data = pd.read_csv(Free_energy_calc_file_Vrabed_reading_name, sep=',')
    Free_energy_Vrabec_calc_data_df = pd.DataFrame(Free_energy_Vrabec_calc_data)

    Temp_Vrabec_calc_data = Free_energy_Vrabec_calc_data_df.iloc[:, 2].values.tolist()
    avg_free_energy_Vrabec_calc_free = Free_energy_Vrabec_calc_data_df.iloc[:, 3].values.tolist()
    std_dev_free_energy_Vrabec_calc_free = Free_energy_Vrabec_calc_data_df.iloc[:, 4].values.tolist()


    # data import NIST experimental
    Free_energy_NIST_exp_data = pd.read_csv(Free_energy_exp_file_NIST_reading_name, sep=',')
    Free_energy_NIST_exp_data_df = pd.DataFrame(Free_energy_NIST_exp_data)

    Temp_NIST_exp_data = Free_energy_NIST_exp_data_df.iloc[:, 2].values.tolist()
    avg_free_energy_NIST_exp_free = Free_energy_NIST_exp_data_df.iloc[:, 3].values.tolist()
    #std_dev_free_energy_NIST_exp_free = Free_energy_NIST_exp_data_df.loc[:, 'std_dev_MBAR_kcal_per_mol'].values.tolist()

    # ********************************
    #  End File importing
    # ********************************




    #****************************************
    #Plot Number 1  (temp vs density) (start)
    #****************************************

    # Plotting curve data below

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1, 1, 1)


    for tick in ax1.xaxis.get_ticklabels():
        tick.set_fontname('Arial')
    for tick in ax1.yaxis.get_ticklabels():
        tick.set_fontname('Arial')

    plt.xlabel('Temperature (K)', fontname="Arial", fontsize=axis_Label_font_size)
    plt.ylabel('Free Energy (kcal/mol)', fontname="Arial", fontsize=axis_Label_font_size)

    avg_free_energy_calc_free_label = 'This Work'
    avg_free_energy_Vrabec_calc_free_label = "Calculated: Linnemann et. al \nSimulations"
    avg_free_energy_NIST_exp_free_label = "Calculated: NIST Website \nExperimental Data"

    # Matplotlib colors b : blue, g : green, r : red, c : cyan, m : magenta, y : yellow, k : black, w : white, gray='x' (x=0 to 1)
    print(f"Temp_calc_data= {Temp_calc_data}")
    print(f"avg_free_energy_calc_free= {avg_free_energy_calc_free}")
    print(f"std_dev_free_energy_calc_free= {std_dev_free_energy_calc_free}")

    print(f"Temp_Vrabec_calc_data= {Temp_Vrabec_calc_data}")
    print(f"avg_free_energy_Vrabec_calc_free= {avg_free_energy_Vrabec_calc_free}")
    print(f"std_dev_free_energy_Vrabec_calc_free= {std_dev_free_energy_Vrabec_calc_free}")

    print(f"Temp_NIST_exp_data= {avg_free_energy_calc_free}")
    print(f"avg_free_energy_NIST_exp_free= {avg_free_energy_NIST_exp_free}")


    plt.errorbar(Temp_calc_data,
                 avg_free_energy_calc_free,
                 yerr=std_dev_free_energy_calc_free,
                 color=Color_for_Plot_Data_calc,
                 marker=marker_gomc,
                 linestyle=Calc_Linestyle,
                 markersize=PointSizes,
                 markeredgewidth=MarkeredgeWidth,
                 linewidth=ConnectionLineSizes,
                 fillstyle='none',
                 capsize=10,
                 label=avg_free_energy_calc_free_label,
                 )
    plt.errorbar(Temp_Vrabec_calc_data,
                 avg_free_energy_Vrabec_calc_free ,
                 yerr=std_dev_free_energy_Vrabec_calc_free,
                 color=Color_for_Plot_Vrabec_calc,
                 marker=marker_Vrabec,
                 linestyle=Calc_Linestyle,
                 markersize=PointSizes,
                 markeredgewidth=MarkeredgeWidth,
                 linewidth=ConnectionLineSizes,
                 fillstyle='none',
                 capsize=10,
                 label=avg_free_energy_Vrabec_calc_free_label
                 )
    plt.plot(Temp_NIST_exp_data,
             avg_free_energy_NIST_exp_free,
             color=Color_for_Plot_NIST_exp,
             marker=None,
             linestyle=Exp_Linestyle,
             markersize=PointSizes,
             markeredgewidth=MarkeredgeWidth,
             linewidth=ConnectionLineSizes,
             fillstyle='none',
             label=avg_free_energy_NIST_exp_free_label,
             )

    major_yticks = np.arange(0, 20+0.001, Free_energy_cal_per_mol_major)
    major_xticks = np.arange(T_min, T_max+0.001, T_ticks_major)

    minor_yticks = np.arange(0, 20+0.001, Free_energy_cal_per_mol_minor )
    minor_xticks = np.arange(T_min, T_max+0.001, T_ticks_minor)


    #plt.gca().set_xlim(left=2, right=105)

    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True)
    ax1.set_yticks(major_yticks)
    ax1.set_yticks(minor_yticks, minor=True)

    ax1.tick_params(axis='both', which='major', length=4, width=2, labelsize=axis_number_font_size, top=True, right=True)
    ax1.tick_params(axis='both', which='minor', length=4, width=1, labelsize=axis_number_font_size, top=True, right=True)

    legend1 = ax1.legend(loc='upper center', shadow=True, fontsize=legend_font_size )

    frame1 = legend1.get_frame()
    frame1.set_facecolor('0.90')
    plt.tight_layout()  # centers layout nice for final paper
    # plt.gcf().subplots_adjust(bottom=0.15) # moves plot up so x label not cutoff

    plt.xlim(T_min, T_max+0.001)  # set plot range on x axis
    plt.gcf().subplots_adjust(left=0.12, bottom=None, right=0.95, top=None, wspace=None, hspace=None) # moves plot  so x label not cutoff

    if solute_i in ["Ne"]:
        plt.ylim(Ne_Free_energy_cal_per_mol_min, Ne_Free_energy_cal_per_mol_max)  # set plot range on x axis
        plt.legend(ncol=1,loc='lower right', fontsize=legend_font_size, prop={'family':'Arial','size': legend_font_size})
    elif solute_i in ["Rn"]:
        plt.ylim(Rn_Free_energy_cal_per_mol_min, Rn_Free_energy_cal_per_mol_max)  # set plot range on x axis
        plt.legend(ncol=1, loc='upper left', fontsize=legend_font_size,
                   prop={'family': 'Arial', 'size': legend_font_size})

    else:
        raise ValueError("ERROR: The solute is not setup for plotting in this script")
    plt.show()
    fig1.savefig(figure_name)

    plt.close()

    #****************************************
    #Plot Number 1  (temp vs density) (end)
    #****************************************




