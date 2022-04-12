#!/usr/bin/python

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import glob
from matplotlib import cm
from scipy.optimize import minimize, brute
from scipy import interpolate, optimize
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle, Ellipse
from shapely.geometry import Point, Polygon
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
    z = np.abs(df4.iloc[:,0].to_numpy())

    rranges = slice(x.min(), x.max(), (x.max() - x.min())/650), slice(y.min(), y.max(), (y.max() - y.min())/650)

    F2 = interpolate.interp2d(x, y, z, kind='quintic')

    xi = np.linspace(x.min(), x.max(), 6500)
    yi = np.linspace(y.min(), y.max(), 6500)


    X,Y = np.meshgrid(xi,yi)

    Z2 = F2(xi, yi)

    f = lambda x: np.abs(F2(*x))


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
        layout = go.Layout(title=title,autosize=False, width=500, height=500, 
        margin=dict(l=65, r=65, b=65, t=65))
        iteractivefig.update_layout(layout)
        iteractivefig.update_layout(scene = dict(
                    xaxis_title='Alpha',
                    yaxis_title='RCut',
                    zaxis_title='Relative Error'),
                    width=700,
                    margin=dict(r=20, b=10, l=10, t=10))
        iteractivefig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))

        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    return (bfXY[0], bfXY[1], ZBF, gdXY[0], gdXY[1], ZGD, gdJacXY[0], gdJacXY[1])


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('3DSurface.py -i <intputfile.p> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('3DSurface.py -i <intputfile.p> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
            parsingInputs = False
    print('Input file path is', inputfile)
    print('Output file is ', outputfile)
    p = re.compile("Wolf_Calibration_(\w+?)_(\w+?)_BOX_(\d+)_(\w+?).dat")

    calibrationFiles = sorted(glob.glob(os.path.join(inputfile,'Wolf_Calibration_*.dat')), key=os.path.getmtime)
    print(calibrationFiles)
    for calFile in calibrationFiles:
        justFileName = os.path.basename(calFile)
        print(justFileName)
        groups = p.search(justFileName)
        wolfKind = groups.group(1)
        potential = groups.group(2)
        box = groups.group(3)
        print ("wolf Kind" , wolfKind)
        print ("potential Kind" , potential)
        print ("box" , box)

        df = pd.read_csv(calFile,sep='\t',index_col=0)
        df = df.iloc[: , :-1]
        dfMean = df.mean()

        points = dfMean.index.map(lambda x: x.strip('('))
        points = points.map(lambda x: x.strip(')'))
        pointsSplit = points.str.split(pat=", ", expand=False)

        df3 = pd.DataFrame(pointsSplit.tolist(), columns=['rcut','alpha'], dtype=np.float64)
        df4 = pd.DataFrame(dfMean.values, columns=['err'], dtype=np.float64)
        print(df3)
        print(df4)

        minxy = df3.min()
        maxxy = df3.max()

        x = df3.iloc[:,0].to_numpy()
        y = df3.iloc[:,1].to_numpy()
        z = np.abs(df4.iloc[:,0].to_numpy())

        x2 = np.linspace(minxy[0], maxxy[0], 6500)
        y2 = np.linspace(minxy[1], minxy[1], 6500)
        print((maxxy[0] - minxy[0]))
        print((maxxy[1] - minxy[1]))
        rranges = slice(minxy[0], maxxy[0], (maxxy[0] - minxy[0])/650), slice(minxy[1], maxxy[1], (maxxy[1] - minxy[1])/650)
        print(rranges)
        X2, Y2 = np.meshgrid(x2, y2)

        X2, Y2 = np.meshgrid(x2, y2)

        F2 = interpolate.interp2d(x, y, z, kind='quintic')

        Z2 = F2(x2, y2)

        f = lambda x: np.abs(F2(*x))


        bounds = [(minxy[0], maxxy[0]),(minxy[1], maxxy[1])]
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


        ZBF = F2(bfXY[0], bfXY[1])
        ZGD = F2(gdXY[0], gdXY[1])

        print("ZBF : ", ZBF)
        print("ZGD : ", ZGD)
        ax = plt.axes(projection='3d')
        ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none');
        ax.set_title(wolfKind+"_"+potential+"_Box_"+box, fontsize=20)
        ax.set_xlabel('Alpha', fontsize=20, labelpad=20)
        ax.set_ylabel('RCut', fontsize=20, labelpad=20)
        ax.set_zlabel('Relative Error', fontsize=20, labelpad=20)

        add_point(ax, gdXY[0], gdXY[1], gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01, labelArg = "Gradient Descent")
        add_point(ax, gdXY[0], gdXY[1], 10*gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01)
        add_point(ax, gdXY[0], gdXY[1], 20*gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01)
        xbf,ybf = bf[0]
        add_point(ax, bfXY[0], bfXY[1], bf[1], fc = 'r', ec = 'r', radius=0.01, labelArg = "Brute Force")
        add_point(ax, bfXY[0], bfXY[1], 10*bf[1], fc = 'r', ec = 'r', radius=0.01)
        add_point(ax, bfXY[0], bfXY[1], 20*bf[1], fc = 'r', ec = 'r', radius=0.01)
        ax.legend(loc='best')
        plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
