#!/usr/bin/python

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import glob
from matplotlib import cm
from scipy.optimize import minimize, brute, shgo, dual_annealing, differential_evolution
from scipy import interpolate, optimize
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle, Ellipse
import pickle
import plotly.io as pio 
import plotly.graph_objects as go
from scipy.interpolate import griddata
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

def find_minimum(path, model, wolfKind, potential, box, plotSuface=False):
    df = pd.read_csv(path,sep='\t',index_col=0)
    df = df.iloc[: , :-1]
    dfMean = df.mean()

    points = dfMean.index.map(lambda x: x.strip('('))
    points = points.map(lambda x: x.strip(')'))
    pointsSplit = points.str.split(pat=", ", expand=False)

    df3 = pd.DataFrame(pointsSplit.tolist(), columns=['rcut','alpha'], dtype=np.float64)
    df4 = pd.DataFrame(dfMean.values, columns=['err'], dtype=np.float64)

    x = df3.iloc[:,0].to_numpy()
    y = df3.iloc[:,1].to_numpy()
    #z = np.abs(df4.iloc[:,0].to_numpy())
    # I wonder if interpolation has problem with abs value
    z = df4.iloc[:,0].to_numpy()

    xi = np.linspace(x.min(), x.max(), 6500)
    yi = np.linspace(y.min(), y.max(), 6500)
    zi = np.linspace(z.min(), z.max(), 6500)

    rranges = slice(x.min(), x.max(), (x.max() - x.min())/650), slice(y.min(), y.max(), (y.max() - y.min())/650)
    print(rranges)

    #F2 = interpolate.interp2d(x, y, z, kind='linear')
    #F2 = interpolate.interp2d(x, y, z, kind='cubic')
    F2 = interpolate.interp2d(x, y, z, kind='quintic')

    X,Y = np.meshgrid(xi,yi)

    Z2 = F2(xi, yi)

    #f = lambda x: np.abs(F2(*x))
    f = lambda x: np.sum(np.abs(F2(np.linspace(*x[0]-0.05, *x[0]+0.05, 10),np.linspace(*x[1]-0.05, *x[1]+0.05, 10))))


    bounds = [(x.min(), x.max()),(y.min(), y.max())]
    bf = brute(f, rranges, full_output=True, finish=optimize.fmin)
    bfXY = np.array(bf[0])
    print(bfXY[0])
    print(bfXY[1])
    x0 = (bfXY[0], bfXY[1])
    gd = minimize(f, x0, method='SLSQP', bounds=bounds)
    print(gd)
    gdXY = np.array(gd.x)
    print(gdXY[0])
    print(gdXY[1])
    gdJacXY = np.array(gd.jac)
    print(gdJacXY[0])
    print(gdJacXY[1])

    print("Calling shgo")
    # Default method is SLSQP
    #shgoOut = shgo(f, x0, method='SLSQP', bounds=bounds)
    # Doesnt work for shgo
    #bounds = [(x.min(), x.max()),(y.min(), y.max())]
    # TypeError: shgo() got multiple values for argument 'bounds'
    # Derivative-free, so can't use gradient here to rank.
    shgoOut = shgo(f, bounds=bounds)
    print(shgoOut)

    print("Calling dual_annealing")
    dual_annealingOut = dual_annealing(f, bounds=bounds, x0=x0)
    print(dual_annealingOut)    

    print("Calling differential_evolution")
    differential_evolutionOut = differential_evolution(f, bounds=bounds, x0=x0)
    print(differential_evolutionOut)  
    print(differential_evolutionOut.keys())   
    quit() 
    """
    differential_evolutionOutXY = np.array(differential_evolutionOut.x)
    print(differential_evolutionOutXY[0])
    print(differential_evolutionOutXY[1])
    differential_evolutionOutJacXY = np.array(differential_evolutionOut.jac)
    print(differential_evolutionOutJacXY[0])
    print(differential_evolutionOutJacXY[1])
    """
    
    ZBF = F2(bfXY[0], bfXY[1])
    ZGD = F2(gdXY[0], gdXY[1])

    d = {'x': [gdXY[0]], 'y': [gdXY[1]], 'z':[ZGD]}

    dfGD = pd.DataFrame(data=d)


    print("ZBF : ", ZBF)
    print("ZGD : ", ZGD)

    if(plotSuface):
        title = model+"_"+wolfKind+"_"+potential+"_Box_"+box
        xi_forplotting = np.linspace(x.min(), x.max(), 1000)
        yi_forplotting = np.linspace(y.min(), y.max(), 1000)

        Z2_forplotting = F2(xi_forplotting, yi_forplotting)

        prefix = os.path.split(path)
        plotPath = os.path.join(prefix[0], title)

        #fig.savefig(fname=plotPath+".png")
        iteractivefig = go.Figure()
        iteractivefig.add_surface(x=xi_forplotting,y=yi_forplotting,z=Z2_forplotting)
        layout = go.Layout(title=title,autosize=True, 
        margin=dict(l=65, r=65, b=65, t=65))
        iteractivefig.update_layout(layout)
        iteractivefig.update_layout(scene = dict(
                    xaxis_title='RCut',
                    yaxis_title='Alpha',
                    zaxis_title='Relative Error'),
                    width=700,
                    margin=dict(r=20, b=10, l=10, t=10))
        iteractivefig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))

        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
