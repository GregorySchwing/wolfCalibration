#!/usr/bin/python
from alchemlyb.parsing.gomc import  extract_dHdl,  extract_u_nk
#from alchemlyb.estimators import MBAR, BAR, TI 
from alchemlyb.estimators import BAR, TI 
from alchemlyb.estimators import AutoMBAR as MBAR
import alchemlyb.preprocessing.subsampling as ss
import pandas as pd
import numpy as np
import os

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import re
import glob
from matplotlib import cm
from scipy.optimize import minimize, brute
from scipy import interpolate, optimize
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle, Ellipse

temprature = 298           #temperature (K)
k_b = 1.9872036E-3         #kcal/mol/K
k_b_T = temprature * k_b

def get_delta(est):
    """ Return the change in free energy and standard deviation for TI and MBAR estimators.

    """
    delta = est.delta_f_.iloc[0, -1] * k_b_T
    d_delta = est.d_delta_f_.iloc[0, -1] * k_b_T
    return delta, d_delta


def get_delta2(est):
    """ Return the change in free energy and standard deviation for BAR estimator.

    """
    ee = 0.0

    for i in range(len(est.d_delta_f_) - 1):
        ee += est.d_delta_f_.values[i][i+1]**2

    delta = est.delta_f_.iloc[0, -1] * k_b_T
    d_delta = k_b_T * ee**0.5
    return delta, d_delta

def add_point(ax, x, y, z, fc = None, ec = None, radius = 0.005, labelArg = None):
	xy_len, z_len = ax.get_figure().get_size_inches()
	axis_length = [x[1] - x[0] for x in [ax.get_xbound(), ax.get_ybound(), ax.get_zbound()]]
	axis_rotation =  {'z': ((x, y, z), axis_length[1]/axis_length[0]),
			'y': ((x, z, y), axis_length[2]/axis_length[0]*xy_len/z_len),
			'x': ((y, z, x), axis_length[2]/axis_length[1]*xy_len/z_len)}
	i = 0
	for a, ((x0, y0, z0), ratio) in axis_rotation.items():
		p = Ellipse((x0, y0), width = radius, height = radius*ratio, fc=fc, ec=ec, label = labelArg if i == 0 else "")
		ax.add_patch(p)
		art3d.pathpatch_2d_to_3d(p, z=z0, zdir=a)
		i = i + 1

##################################################

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('FreeEnergyAnalysis.py -i <intputfile.p> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('FreeEnergyAnalysis.py -i <intputfile.p> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
            parsingInputs = False
    print('Input file path is', inputfile)
    print('Output file is ', outputfile)
    p = re.compile("Free_Energy_BOX_(\d+)_(\w+?).dat")    
    
    list_data_TI = []
    list_data_BAR = []
    tis = []
    mbars = []
    bars = []

    for root, dirs, files in os.walk(inputfile):
        path = root.split(os.sep)
        for file in files:
            matched = p.match(file)
            if (bool(matched)):
                print("Reading File: %s " % os.path.join(root, file))
                            
                # Read the data for TI estimator and BAR or MBAR estimators.
                dHdl = extract_dHdl(os.path.join(root, file), T=temprature)
                u_nkr = extract_u_nk(os.path.join(root, file), T=temprature)
                #Detect uncorrelated samples using VDW+Coulomb term in derivative 
                # of energy time series (calculated for TI)
                srs = dHdl['VDW'] + dHdl['Coulomb'] 
                list_data_TI.append(ss.statistical_inefficiency(dHdl, series=srs, conservative=False))
                list_data_BAR.append(ss.statistical_inefficiency(u_nkr, series=srs, conservative=False))
                #print(dHdl)
                #print(u_nkr)
                #print(srs)


    #for TI estimator
    #print("Working on TI method ...")
    dhdl = pd.concat([ld for ld in list_data_TI])
    print(dhdl)
    ti = TI().fit(dhdl)
    sum_ti, sum_ds_ti = get_delta(ti)
    tis.append(sum_ti)


    #for MBAR estimator
    #print("Working on MBAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    print(u_nk)
    mbar = MBAR().fit(u_nk)
    sum_mbar, sum_ds_mbar = get_delta(mbar)
    mbars.append(sum_mbar)

    #for BAR estimator
    #print("Working on BAR method ...")
    u_nk = pd.concat([ld for ld in list_data_BAR])
    bar = BAR().fit(u_nk)
    sum_bar, sum_ds_bar = get_delta2(bar)
    bars.append(sum_bar)

    nr = 1
    com = "SPCE"
    print("Replica-%d, %4s,  %7.4f,  %7.4f,  %7.4f,  %7.4f,  %7.4f,  %7.4f" % (nr, com, sum_ti, sum_ds_ti, sum_mbar, sum_ds_mbar, sum_bar, sum_ds_bar))

        
    tis = np.array(tis)
    mbars = np.array(mbars)
    bars = np.array(bars)
    print("Average, %4s,  %7.4f,  %7.4f,  %7.4f,  %7.4f,  %7.4f,  %7.4f" % (com, np.average(tis), np.std(tis), np.average(mbars), np.std(mbars), np.average(bars), np.std(bars)))        

if __name__ == "__main__":
    main(sys.argv[1:])
