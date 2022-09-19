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
from scipy.interpolate import Rbf, interp2d, RegularGridInterpolator
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

#https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
def reject_outliers(data, m=2):
    return [abs(data - np.mean(data)) < m * np.std(data)]

# Bad function, doesn't work
def reject_outliers_median(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

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

    x = np.unique(x)
    y = np.unique(y)

    #z = np.abs(df4.iloc[:,0].to_numpy())
    # I wonder if interpolation has problem with abs value
    z = df4.iloc[:,0].to_numpy()

    z = np.reshape(z, (len(x),len(y)))

    print(x)
    print(y)
    print(z)

    rranges = slice(x.min(), x.max(), (x.max() - x.min())/650), slice(y.min(), y.max(), (y.max() - y.min())/650)
    print(rranges)

    sptbf_mins = {}
    sptgd_mins = {}
    sptbf_auc = {}
    sptgd_auc = {}

    sptshgo_mins = {}
    sptda_mins = {}
    sptde_mins = {}
    sptdanx_mins = {}
    sptdenx_mins = {}

    sptshgo_auc = {}
    sptda_auc = {}
    sptde_auc = {}
    sptdanx_auc = {}
    sptdenx_auc = {}
    #F2 = interpolate.interp2d(x, y, z, kind='linear')
    #F2 = interpolate.RegularGridInterpolator((y_raw,x_raw), z_raw)
    F2 = interpolate.RegularGridInterpolator(points=(x,y), values=z, method='linear', bounds_error=False, fill_value=None)
    
   
    """
    fig = plt.figure()

    ax = fig.add_subplot(projection='3d')
    
    xx = np.linspace(x.min(), x.max(), 650)
    yy = np.linspace(y.min(), y.max(), 650)

    X, Y = np.meshgrid(xx, yy, indexing='ij')
    ax.plot_wireframe(X, Y, F2((X, Y), method='linear'), rstride=3, cstride=3,
                  alpha=0.4, color='m', label='linear interp')


    plt.legend()

    plt.show()
    """


    bounds = [(x.min(), x.max()),(y.min(), y.max())]
    f = lambda x: np.abs(F2(xi= x, method='linear'))
    bf = brute(f, rranges, full_output=True, finish=None)
    bfXY = np.array(bf[0])
    print("best rcut", bfXY[0])
    print("best alpha", bfXY[1])
  
    #x0 = (12, 0.12)
    x0 = (bfXY[0], bfXY[1])
    gd = minimize(f, x0, method='SLSQP', bounds=bounds)
    print(gd)
    gdXY = np.array(gd.x)
    print(gdXY[0])
    print(gdXY[1])
    gdJacXY = np.array(gd.jac)
    print(gdJacXY[0])
    print(gdJacXY[1])

    sptbf_mins["REF"] = bfXY
    sptbf_auc["REF"] = bf[1]
    sptgd_mins["REF"] = gd.x
    sptgd_auc["REF"] = gd.fun


    ZBF = F2(bfXY, method='linear')
    ZGD = F2(gd.x, method='linear')

    d = {'x': [gdXY[0]], 'y': [gdXY[1]], 'z':[ZGD]}


    print("Calling shgo")
    # Default method is SLSQP
    #shgoOut = shgo(f, x0, method='SLSQP', bounds=bounds)
    # Doesnt work for shgo
    #bounds = [(x.min(), x.max()),(y.min(), y.max())]
    # TypeError: shgo() got multiple values for argument 'bounds'
    # Derivative-free, so can't use gradient here to rank.
    sptshgoOut = shgo(f, bounds=bounds)
    print(sptshgoOut)
    sptshgo_mins["REF"] = sptshgoOut.x
    sptshgo_auc["REF"] = sptshgoOut.fun
    
    print("Calling dual_annealing")
    sptdual_annealingOutNoX0 = dual_annealing(f, bounds=bounds)
    sptdual_annealingOut = dual_annealing(f, bounds=bounds, x0=x0)
    print(sptdual_annealingOut)    
    print("Calling dual_annealingNoX0")
    print(sptdual_annealingOutNoX0)    
    sptda_mins["REF"] = sptdual_annealingOut.x
    sptdanx_mins["REF"] = sptdual_annealingOutNoX0.x

    sptda_auc["REF"] = sptdual_annealingOut.fun
    sptdanx_auc["REF"] = sptdual_annealingOutNoX0.fun

    print("Calling differential_evolution")
    sptdifferential_evolutionOut = differential_evolution(f, bounds=bounds, x0=x0)
    sptdifferential_evolutionOutNox0 = differential_evolution(f, bounds=bounds)
    print(sptdifferential_evolutionOut)  
    print("sptdifferential_evolutionOut.keys()", sptdifferential_evolutionOut.keys())  
    print(sptdifferential_evolutionOut.keys())   
    print("Calling differential_evolutionX0")
    print(sptdifferential_evolutionOutNox0) 
    sptde_mins["REF"] = sptdifferential_evolutionOut.x
    sptdenx_mins["REF"] = sptdifferential_evolutionOutNox0.x

    sptde_auc["REF"] = sptdifferential_evolutionOut.fun
    sptdenx_auc["REF"] = sptdifferential_evolutionOutNox0.fun



    dfGD = pd.DataFrame(data=d)


    print("ZBF : ", ZBF)
    print("ZGD : ", ZGD)

    print("sptbf_mins", sptbf_mins)
    print("sptgd_mins", sptgd_mins)
    print("sptshgo_mins", sptshgo_mins)
    print("sptda_mins", sptda_mins)
    print("sptde_mins", sptde_mins)
    print("sptdanx_mins", sptdanx_mins)
    print("sptdenx_mins", sptdenx_mins)


    print("sptbf_auc", sptbf_auc)
    print("sptgd_auc", sptgd_auc)
    print("sptshgo_auc", sptshgo_auc)
    print("sptda_auc", sptda_auc)
    print("sptde_auc", sptde_auc)
    print("sptdanx_auc", sptdanx_auc)
    print("sptdenx_auc", sptdenx_auc)

    goMethods = {}
    goMethods["sptbf"] = sptbf_mins
    goMethods["sptgd"] = sptgd_mins
    goMethods["sptshgo"] = sptshgo_mins
    goMethods["sptda"] = sptda_mins
    goMethods["sptde"] = sptde_mins
    goMethods["sptdanx"] = sptdanx_mins
    goMethods["sptdenx"] = sptdenx_mins

    goAUCs = {}

    # AUC of single points are not what I'm using here.

    goAUCs["sptbf"] = sptbf_auc
    goAUCs["sptgd"] = sptgd_auc
    goAUCs["sptshgo"] = sptshgo_auc
    goAUCs["sptda"] = sptda_auc
    goAUCs["sptde"] = sptde_auc
    goAUCs["sptdanx"] = sptdanx_auc
    goAUCs["sptdenx"] = sptdenx_auc

    smallestGrad = 100000000
    winningOptimizer = ""
    sizeOfRegionScale = 0.001
    for key, value in goMethods.items():
        sizeOfRegionX = sizeOfRegionScale*(x.max()-x.min())
        sizeOfRegionY = sizeOfRegionScale*(y.max()-y.min())
    
        xxx = np.linspace(value["REF"][0]-sizeOfRegionX, value["REF"][0]+sizeOfRegionX, 10)
        yyy = np.linspace(value["REF"][1]-sizeOfRegionY, value["REF"][1]+sizeOfRegionY, 10)


        Xpoint, Ypoint = np.meshgrid(xxx, yyy, indexing='ij')
        print("method",key, value)
        print("calculating gradient of points centered at ",value["REF"][0], " " , value["REF"][1])

        print(Xpoint.ravel())
        print(Ypoint.ravel())
        print(F2((Xpoint, Ypoint), method='linear').ravel())
        grad = np.gradient(np.stack([Xpoint.ravel(), Ypoint.ravel(), F2((Xpoint, Ypoint), method='linear').ravel()]))
        print("gradient of ", key, " : ", grad)
        gradNorm = LA.norm(grad)
        print("LA.norm(grad) ", gradNorm)
        
        if (gradNorm < smallestGrad):
            winningOptimizer = key
            smallestGrad = gradNorm

        
    if(plotSuface):
        title = model+"_"+wolfKind+"_"+potential+"_Box_"+box

        prefix = os.path.split(path)
        plotPath = os.path.join(prefix[0], title)
        xx_forplotting = np.linspace(x.min(), x.max(), 1000)
        yy_forplotting = np.linspace(y.min(), y.max(), 1000)

        X_forplotting, Y_forplotting = np.meshgrid(xx_forplotting, yy_forplotting, indexing='xy')

        iteractivefig = go.Figure()
        #iteractivefig.add_surface(autocolorscale=True, x=X, y=Y, z=F2((X, Y), method='linear'))
        iteractivefig.add_surface(autocolorscale=True, x=xx_forplotting, y=yy_forplotting, z=F2((X_forplotting, Y_forplotting), method='linear'))
        #iteractivefig.add_surface(autocolorscale=True, x=X.ravel(), y=Y.ravel(), z=F2((X, Y), method='linear').ravel())
        #iteractivefig.add_surface(x=xi_forplotting,y=yi_forplotting,z=Z2_forplotting)
        layout = go.Layout(title=title,autosize=True, margin=dict(l=65, r=65, b=65, t=65))
        iteractivefig.update_layout(layout)
        iteractivefig.update_layout(scene = dict(
                    xaxis_title='RCut',
                    yaxis_title='Alpha',
                    zaxis_title='Relative Error'),
                    width=700,
                    margin=dict(r=20, b=10, l=10, t=10))
        iteractivefig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))
        
        for key, value in goMethods.items():
            print("method",key, value)
            xvals = [item[0] for item in value.values()]
            yvals = [item[1] for item in value.values()]
            zvals = []
            for x,y in zip(xvals,yvals):
                zvals.append(F2((x, y), method='linear'))
            print("x:", xvals)
            print("y:", yvals)
            print("z:", zvals)
            iteractivefig.add_trace(
                go.Scatter3d(x=xvals,
                            y=yvals,
                            z=zvals,
                            mode='markers',
                            name=key,
                            hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                            showlegend=True)
            )
        
        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    print("GD_rcut",goMethods[winningOptimizer]["REF"][0])
    print("GD_alpha",goMethods[winningOptimizer]["REF"][1])
    print("GD_relerr",goAUCs[winningOptimizer]["REF"] )
    print("GD_grad",smallestGrad) 
    print("WINNING_OPT",winningOptimizer)
    return ( ("GD_rcut",goMethods[winningOptimizer]["REF"][0]), ("GD_alpha",goMethods[winningOptimizer]["REF"][1]), ("GD_relerr",goAUCs[winningOptimizer]["REF"][0]), ("GD_grad",smallestGrad),  ("WINNING_OPT",winningOptimizer) )
