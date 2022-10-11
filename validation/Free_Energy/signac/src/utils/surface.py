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
import scipy
from scipy.optimize import OptimizeResult, least_squares
from  scipy.interpolate import UnivariateSpline
def merge(list1, list2):
      
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list

def poly_fun(coeffs, a, x):
    predicted = 10**np.polynomial.polynomial.polyval(np.log10(a), coeffs)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, predicted) 
    return np.array([3, 4, 5]) * np.array([slope-1, r_value-1, intercept]) 

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

from pymoo.core.problem import ElementwiseProblem
from scipy.interpolate import RectBivariateSpline
from scipy import interpolate
class MyProblem(ElementwiseProblem):

    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax, LowerBoundRcut):
        super().__init__(n_var=2,
                         n_obj=3,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.tck_pd = tck_pd
        self.RCutMin = RCutMin
        self.LowerBoundRcut = LowerBoundRcut
    def _evaluate(self, x, out, *args, **kwargs):
        # Minimize Relative Error
        f1 = np.abs(self.rect_B_spline.ev(x[0], x[1]))
        # Minimize RCut
        f2 = (x[0]-self.RCutMin)
        # Minimize Gradient
        f3 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_pd))

        g1 = -(x[0]-self.LowerBoundRcut)

        out["F"] = [f1, f2, f3]
        out["G"] = [g1]



def find_minimum(path, model, wolfKind, potential, box, plotSuface=False):
    df = pd.read_csv(path,sep='\t',index_col=0)
    # remove the nan column
    df = df.iloc[: , :-1]
    dfMean = df.iloc[0, :]

    #dfMean = df.mean()
    print(dfMean)
    #quit()
    print("columns",df.columns)
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
    print(x)
    print(y)
    print(z)
    print("lenx ", len(x))
    print("leny ", len(y))

    z = np.reshape(z, (len(x),len(y)))
    #z = np.reshape(z, (len(y),len(x)))

    print(x)
    print(y)
    print(z)

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
    from scipy.interpolate import RectBivariateSpline
    from scipy import interpolate

    rect_B_spline = RectBivariateSpline(x, y, z)
    pd_RCut_varies_alpha_constant = [1,0]
    pd_RCut_constant_alpha_varies = [0,1]

    # OK - M.O. 2.0
    derivs = rect_B_spline.partial_derivative(pd_RCut_varies_alpha_constant[0],pd_RCut_varies_alpha_constant[1])
    # Bad - don't use this
    #derivs = rect_B_spline.partial_derivative(pd_RCut_constant_alpha_varies[0],pd_RCut_constant_alpha_varies[1])

    tck_pd = [derivs.tck[0], derivs.tck[1],derivs.tck[2],derivs.degrees[0],derivs.degrees[1]]
 
    print("Derivative data:")
    print(tck_pd)
    #F2 = interpolate.RegularGridInterpolator(points=(x,y), values=z, method='linear', bounds_error=True, fill_value=None)
    #f = lambda x: np.abs(F2(xi= x, method='linear'))
    exampleX = 10
    exampleY = 0.16
    val = interpolate.bisplev(exampleX, exampleY, tck_pd)
    print("deriv at ", exampleX, exampleY)
    print(val)
    exampleX2 = 10
    exampleY2 = 0.0
    val = interpolate.bisplev(exampleX2, exampleY2, tck_pd)
    print("deriv", exampleX2, exampleY2)
    print(val)


    problem = MyProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 10)

    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.operators.crossover.sbx import SBX
    from pymoo.operators.mutation.pm import PM
    from pymoo.operators.sampling.rnd import FloatRandomSampling

    algorithm = NSGA2(
        pop_size=40,
        n_offsprings=10,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True
    )

    from pymoo.termination import get_termination

    termination = get_termination("n_gen", 400)

    from pymoo.optimize import minimize

    res = minimize(problem,
                algorithm,
                termination,
                seed=1,
                save_history=True,
                verbose=True)

    X = res.X
    F = res.F

    print(X)
    print(F)

    import matplotlib.pyplot as plt
    xl, xu = problem.bounds()
    plt.figure(figsize=(7, 5))
    plt.scatter(X[:, 0], X[:, 1], s=30, facecolors='none', edgecolors='r')
    plt.xlim(xl[0], xu[0])
    plt.ylim(xl[1], xu[1])
    plt.title("Design Space")
    plt.show()

    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.title("Objective Space")
    plt.show()

    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)

    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(approx_ideal[0], approx_ideal[1], facecolors='none', edgecolors='red', marker="*", s=100, label="Ideal Point (Approx)")
    plt.scatter(approx_nadir[0], approx_nadir[1], facecolors='none', edgecolors='black', marker="p", s=100, label="Nadir Point (Approx)")
    plt.title("Objective Space")
    plt.legend()
    plt.show()

    nF = (F - approx_ideal) / (approx_nadir - approx_ideal)

    fl = nF.min(axis=0)
    fu = nF.max(axis=0)
    print(f"Scale f1: [{fl[0]}, {fu[0]}]")
    print(f"Scale f2: [{fl[1]}, {fu[1]}]")

    ### last working line

    plt.figure(figsize=(7, 5))
    plt.scatter(nF[:, 0], nF[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.title("Objective Space")
    plt.show()
    
    # if you use MO 1.0
    #weights = np.array([0.2, 0.8])
    weights = np.array([0.6, 0.1, 0.3])



    from pymoo.decomposition.asf import ASF

    decomp = ASF()

    i = decomp.do(nF, 1/weights).argmin()

    print("Best regarding ASF: Point \ni = %s\nF = %s" % (i, F[i]))
    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(F[i, 0], F[i, 1], marker="x", color="red", s=200)
    plt.title("Objective Space")
    plt.show()
    from pymoo.mcdm.pseudo_weights import PseudoWeights

    i = PseudoWeights(weights).do(nF)

    print("Best regarding Pseudo Weights: Point \ni = %s\nF = %s" % (i, F[i]))
    print(X[i])

    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(F[i, 0], F[i, 1], marker="x", color="red", s=200)
    plt.title("Objective Space")
    plt.show()


    if(plotSuface):
        title = model+"_"+wolfKind+"_"+potential+"_Box_"+box

        prefix = os.path.split(path)
        plotPath = os.path.join(prefix[0], title)
        xx_forplotting = np.linspace(x.min(), x.max(), 1000)
        yy_forplotting = np.linspace(y.min(), y.max(), 1000)

        #X_forplotting, Y_forplotting = np.meshgrid(xx_forplotting, yy_forplotting, indexing='xy')

        iteractivefig = go.Figure()
        #iteractivefig.add_surface(autocolorscale=True, x=X, y=Y, z=F2((X, Y), method='linear'))
        iteractivefig.add_surface(autocolorscale=True, x=xx_forplotting, y=yy_forplotting, z=rect_B_spline(xx_forplotting, yy_forplotting, grid=True))
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
        """
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
        """
        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)
    quit()

    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    print("GD_rcut",goMethods[winningOptimizer]["REF"][0])
    print("GD_alpha",goMethods[winningOptimizer]["REF"][1])
    print("GD_relerr",goAUCs[winningOptimizer]["REF"] )
    print("WINNING_OPT",winningOptimizer)
    return ( ("GD_rcut",goMethods[winningOptimizer]["REF"][0]), ("GD_alpha",goMethods[winningOptimizer]["REF"][1]), ("GD_relerr",F2(goMethods[winningOptimizer]["REF"])), ("WINNING_OPT",winningOptimizer) )