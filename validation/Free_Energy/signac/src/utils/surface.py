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

def pareto_frontier_multi(myArray):
    # Sort on first dimension
    myArray = myArray[myArray[:,0].argsort()]
    # Add first row to pareto_frontier
    pareto_frontier = myArray[0:1,:]
    # Test next row against the last row in pareto_frontier
    for row in myArray[1:,:]:
        if sum([row[x] >= pareto_frontier[-1][x]
                for x in range(len(row))]) == len(row):
            # If it is better on all features add the row to pareto_frontier
            pareto_frontier = np.concatenate((pareto_frontier, [row]))
    return pareto_frontier

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

def keep_efficient(pts):
    'returns Pareto efficient row subset of pts'
    # sort points by decreasing sum of coordinates
    pts = pts[pts.sum(1).argsort()[::-1]]
    # initialize a boolean mask for undominated points
    # to avoid creating copies each iteration
    undominated = np.ones(pts.shape[0], dtype=bool)
    for i in range(pts.shape[0]):
        # process each point in turn
        n = pts.shape[0]
        if i >= n:
            break
        # find all points not dominated by i
        # since points are sorted by coordinate sum
        # i cannot dominate any points in 1,...,i-1
        undominated[i+1:n] = (pts[i+1:] >= pts[i]).any(1) 
        # keep points undominated so far
        pts = pts[undominated[:n]]
    return pts

# Faster than is_pareto_efficient_simple, but less readable.
def is_pareto_efficient(costs, return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient

def cull(pts, dominates):
    dominated = []
    cleared = []
    remaining = pts
    while remaining:
        candidate = remaining[0]
        new_remaining = []
        for other in remaining[1:]:
            [new_remaining, dominated][dominates(candidate, other)].append(other)
        if not any(dominates(other, candidate) for other in new_remaining):
            cleared.append(candidate)
        else:
            dominated.append(candidate)
        remaining = new_remaining
    return cleared, dominated

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
from pymoo.util.misc import stack
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
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.LowerBoundRcut = LowerBoundRcut
        
    def _evaluate(self, x, out, *args, **kwargs):
        # Minimize Relative Error
        f1 = np.abs(self.rect_B_spline.ev(x[0], x[1]))
        # Minimize RCut
        f2 = (x[0]-self.RCutMin)
        # Minimize Gradient
        f3 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_pd))

        g1 = -(x[0]-self.LowerBoundRcut)
        #out["F"] = [f1, f2]

        out["F"] = [f1, f2, f3]
        out["G"] = [g1]

    def _calc_pareto_front(self, flatten=True, *args, **kwargs):

        num_pts = 100
        rcuts = np.linspace(self.RCutMin, self.RCutMax, num_pts)
        alphas = np.linspace(self.AlphaMin, self.AlphaMax, num_pts)
        rcuts_g, alphas_g = np.meshgrid(rcuts, alphas)

        points = np.array(list(zip(rcuts_g.ravel(), alphas_g.ravel())))


        F1_a_costs = np.array(np.abs(self.rect_B_spline.ev(rcuts_g.ravel(), alphas_g.ravel())))
        F2_a_costs = np.array(rcuts_g.ravel())-self.RCutMin
        derivLambda = lambda x : np.abs(interpolate.bisplev(x[0], x[1], self.tck_pd))

        F3_a_costs = np.array([derivLambda(i) for i in points])
        F3_a_costs = F3_a_costs
        costs = np.array(list(zip(F1_a_costs, F2_a_costs, F3_a_costs)))
        pass2Method = np.array(list(zip(costs,points)))
        boolean_array = is_pareto_efficient(costs)
        print(boolean_array)
        #pf_a = pareto_frontier_multi(costs)
        return costs[boolean_array]

class MyProblemNorm(ElementwiseProblem):

    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax, LowerBoundRcut, ParetoFront):
        super().__init__(n_var=2,
                         n_obj=3,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.tck_pd = tck_pd
        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.LowerBoundRcut = LowerBoundRcut
        self.ParetoFront = ParetoFront
        self.F1Min = 0
        self.F2Min = 0
        self.F3Min = 0
        self.F1Max = ParetoFront[:, 0].max()
        self.F2Max = ParetoFront[:, 1].max()
        self.F3Max = ParetoFront[:, 2].max()
                
    def _evaluate(self, x, out, *args, **kwargs):
        # Minimize Relative Error
        f1 = np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.F1Max
        # Minimize RCut
        f2 = (x[0]-self.RCutMin)/self.F2Max
        # Minimize Gradient
        f3 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_pd))/self.F3Max

        g1 = -(x[0]-self.LowerBoundRcut)
        #out["F"] = [f1, f2]

        out["F"] = [f1, f2, f3]
        out["G"] = [g1]

    def _calc_pareto_front(self, flatten=True, *args, **kwargs):

        num_pts = 100
        rcuts = np.linspace(self.RCutMin, self.RCutMax, num_pts)
        alphas = np.linspace(self.AlphaMin, self.AlphaMax, num_pts)
        rcuts_g, alphas_g = np.meshgrid(rcuts, alphas)

        points = np.array(list(zip(rcuts_g.ravel(), alphas_g.ravel())))


        F1_a_costs = np.array((np.abs(self.rect_B_spline.ev(rcuts_g.ravel(), alphas_g.ravel())))/self.F1Max)
        F2_a_costs = ((rcuts_g.ravel()-self.RCutMin)/self.F2Max)
        derivLambda = lambda x : np.abs(interpolate.bisplev(x[0], x[1], self.tck_pd))

        F3_a_costs = np.array([derivLambda(i) for i in points])
        F3_a_costs = F3_a_costs/self.F3Max
        costs = np.array(list(zip(F1_a_costs, F2_a_costs, F3_a_costs)))
        pass2Method = np.array(list(zip(costs,points)))
        boolean_array = is_pareto_efficient(costs)
        print(boolean_array)
        #pf_a = pareto_frontier_multi(costs)
        return costs[boolean_array]

def find_minimum(path, model, wolfKind, potential, box, plotSuface=False):
    df = pd.read_csv(path,sep='\t',index_col=0)
    # remove the nan column
    df = df.iloc[: , :-1]
    
    # You can test Vlugt's assertion that 
    
    """
    For dense liquids such as water and methanol at ambient conditions, 
    it is suﬃcient to plot Figure 1 for a single conﬁguration [93].
    
    dfMean = df.iloc[0, :]

    """

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
    #This is wrong : z = np.reshape(z, (len(y),len(x)))


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
    pd_RCut_varies_alpha_varies = [1,1]

    # Bad - don't use this
    #derivs = rect_B_spline.partial_derivative(pd_RCut_varies_alpha_constant[0],pd_RCut_varies_alpha_constant[1])
    # OK - M.O. 2.0
    #derivs = rect_B_spline.partial_derivative(pd_RCut_constant_alpha_varies[0],pd_RCut_constant_alpha_varies[1])
    # M.O. 3.0
    derivs = rect_B_spline.partial_derivative(pd_RCut_varies_alpha_varies[0],pd_RCut_varies_alpha_varies[1])


    #smallestZDeriv = np.min(derivs, axis=0)
    #largestZDeriv = np.min(derivs, axis=0)
    
    print("smallestZDeriv", derivs)
    #print("largestZDeriv", largestZDeriv)
    tck_pd = [derivs.tck[0], derivs.tck[1],derivs.tck[2],derivs.degrees[0],derivs.degrees[1]]
 
    print("Derivative data:")
    print(tck_pd)
    #F2 = interpolate.RegularGridInterpolator(points=(x,y), values=z, method='linear', bounds_error=True, fill_value=None)
    #f = lambda x: np.abs(F2(xi= x, method='linear'))
    exampleX = 10.05
    exampleY = 0.16
    val = interpolate.bisplev(exampleX, exampleY, tck_pd)
    print("deriv at ", exampleX, exampleY)
    print(val)
    exampleX2 = 10.05
    exampleY2 = 0.0
    val = interpolate.bisplev(exampleX2, exampleY2, tck_pd)
    print("deriv", exampleX2, exampleY2)
    print(val)
    
    # Create problem to get the unnormalized Pareto Front
    problemUnNorm = MyProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 10)
    pf_for_norm = problemUnNorm.pareto_front(use_cache=False, flatten=False)
    problem = MyProblemNorm(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 10, pf_for_norm)

    # Gross p value ~ 0.599426150135242
    #if (wolfKind == "GROSS"):
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.operators.crossover.sbx import SBX
    from pymoo.operators.mutation.pm import PM
    from pymoo.operators.sampling.rnd import FloatRandomSampling

    algorithm = NSGA2(
        pop_size=1000,
        n_offsprings=100,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True
    )

    from pymoo.termination import get_termination
    termination = get_termination("n_gen", 1000)

    from pymoo.optimize import minimize

    res = minimize(problem,
                algorithm,
                termination,
                seed=1,
                save_history=True,
                verbose=True)

    X = res.X
    F = res.F
    hist = res.history
    print(X)
    print(F)

    n_evals = []             # corresponding number of function evaluations\
    hist_F = []              # the objective space values in each generation
    hist_cv = []             # constraint violation in each generation
    hist_cv_avg = []         # average constraint violation in the whole population

    for algo in hist:

        # store the number of function evaluations
        n_evals.append(algo.evaluator.n_eval)

        # retrieve the optimum from the algorithm
        opt = algo.opt

        # store the least contraint violation and the average in each population
        hist_cv.append(opt.get("CV").min())
        hist_cv_avg.append(algo.pop.get("CV").mean())

        # filter out only the feasible and append and objective space values
        feas = np.where(opt.get("feasible"))[0]
        hist_F.append(opt.get("F")[feas])

    k = np.where(np.array(hist_cv) <= 0.0)[0].min()
    print(f"At least one feasible solution in Generation {k} after {n_evals[k]} evaluations.")


    # replace this line by `hist_cv` if you like to analyze the least feasible optimal solution and not the population
    vals = hist_cv_avg

    k = np.where(np.array(vals) <= 0.0)[0].min()
    print(f"Whole population feasible in Generation {k} after {n_evals[k]} evaluations.")

    """
    plt.figure(figsize=(7, 5))
    plt.plot(n_evals, vals,  color='black', lw=0.7, label="Avg. CV of Pop")
    plt.scatter(n_evals, vals,  facecolor="none", edgecolor='black', marker="p")
    plt.axvline(n_evals[k], color="red", label="All Feasible", linestyle="--")
    plt.title("Convergence")
    plt.xlabel("Function Evaluations")
    plt.legend()
    plt.show()
    """

    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)

    from pymoo.indicators.hv import Hypervolume

    metric = Hypervolume(ref_point= np.array([1.1, 1.1, 1.1]),
                        norm_ref_point=False,
                        zero_to_one=True,
                        ideal=approx_ideal,
                        nadir=approx_nadir)

    hv = [metric.do(_F) for _F in hist_F]

    plt.figure(figsize=(7, 5))
    plt.plot(n_evals, hv,  color='black', lw=0.7, label="Avg. CV of Pop")
    plt.scatter(n_evals, hv,  facecolor="none", edgecolor='black', marker="p")
    plt.title("Convergence")
    plt.xlabel("Function Evaluations")
    plt.ylabel("Hypervolume")
    #plt.show()
    
    titleHyperVolume = model+"_"+wolfKind+"_"+potential+"_Box_"+box+"HyperVolume"
    prefix = os.path.split(path)
    HyperVolumeFigPath = os.path.join(prefix[0], titleHyperVolume)
    plt.savefig(HyperVolumeFigPath)


    xl, xu = problem.bounds()
    
    """
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
    """
    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)
    """
    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(approx_ideal[0], approx_ideal[1], facecolors='none', edgecolors='red', marker="*", s=100, label="Ideal Point (Approx)")
    plt.scatter(approx_nadir[0], approx_nadir[1], facecolors='none', edgecolors='black', marker="p", s=100, label="Nadir Point (Approx)")
    plt.title("Objective Space")
    plt.legend()
    plt.show()
    """
    nF = (F - approx_ideal) / (approx_nadir - approx_ideal)

    fl = nF.min(axis=0)
    fu = nF.max(axis=0)
    print(f"Scale f1: [{fl[0]}, {fu[0]}]")
    print(f"Scale f2: [{fl[1]}, {fu[1]}]")

    """
    plt.figure(figsize=(7, 5))
    plt.scatter(nF[:, 0], nF[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.title("Objective Space")
    plt.show()
    """
    pf_a = problem.pareto_front(use_cache=False, flatten=False)

    # Creating figure
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    # Creating plot
    ax.scatter3D(pf_a[:, 0], pf_a[:, 1], pf_a[:, 2], color = "red", label="Pareto-front")
    ax.scatter3D(F[:, 0], F[:, 1], F[:, 2], color = "blue", label="Functionals")
    plt.title("Objective Space")
    plt.legend()

    titleParetoFront = model+"_"+wolfKind+"_"+potential+"_Box_"+box+"_ParetoFront"
    prefix = os.path.split(path)
    paretoFrontFigPath = os.path.join(prefix[0], titleParetoFront)
    plt.savefig(paretoFrontFigPath)
    
    sbsp=interpolate.SmoothBivariateSpline(pf_a[:, 0], pf_a[:, 1], pf_a[:, 2])
    pf_x_min = (pf_a[:, 0]).min()
    pf_x_max = (pf_a[:, 0]).max()
    pf_y_min = (pf_a[:, 1]).min()
    pf_y_max = (pf_a[:, 1]).max()
    pf_z_min = (pf_a[:, 2]).min()
    pf_z_max = (pf_a[:, 2]).max()
    
    print("pf_x_min",pf_x_min)
    print("pf_x_max",pf_x_max)
    print("pf_y_min",pf_y_min)
    print("pf_y_max",pf_y_max)
    print("pf_z_min",pf_z_min)
    print("pf_z_max",pf_z_max)
    xx_pareto = np.linspace(pf_x_min, pf_x_max, 1000)
    yy_pareto = np.linspace(pf_y_min, pf_y_max, 1000)

    X_pareto_forplotting, Y_pareto_forplotting = np.meshgrid(xx_pareto, yy_pareto, indexing="ij")

    # Pareto front is an irregular grid.
    from scipy.interpolate import griddata
    grid_z0 = griddata(list(zip(pf_a[:, 0], pf_a[:, 1])), pf_a[:, 2], (X_pareto_forplotting, Y_pareto_forplotting), method='nearest')
    grid_z1 = griddata(list(zip(pf_a[:, 0], pf_a[:, 1])), pf_a[:, 2], (X_pareto_forplotting, Y_pareto_forplotting), method='linear')

    grid_z2 = griddata(list(zip(pf_a[:, 0], pf_a[:, 1])), pf_a[:, 2], (X_pareto_forplotting, Y_pareto_forplotting), method='cubic')
    """
    plt.subplot(221)

    #plt.imshow(func(X_pareto_forplotting, Y_pareto_forplotting).T, extent=(0,1,0,1), origin='lower')

    plt.plot(pf_a[:,0], pf_a[:,1], 'k.', ms=1)

    plt.title('Original')

    plt.subplot(222)

    plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')

    plt.title('Nearest')

    plt.subplot(223)

    plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')

    plt.title('Linear')

    plt.subplot(224)

    plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')

    plt.title('Cubic')

    plt.gcf().set_size_inches(6, 6)

    plt.show()
    """

    iteractivefig = go.Figure()
    iteractivefig.add_surface(autocolorscale=True, x=X_pareto_forplotting, y=Y_pareto_forplotting, z=grid_z2)
    layout = go.Layout(title=titleParetoFront,autosize=True, margin=dict(l=65, r=65, b=65, t=65))
    iteractivefig.update_layout(layout)
    iteractivefig.update_layout(scene = dict(
                xaxis_title='F1',
                yaxis_title='F2',
                zaxis_title='F3'),
                width=700,
                margin=dict(r=20, b=10, l=10, t=10))
    iteractivefig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                highlightcolor="limegreen", project_z=True))
    iteractivefig.add_trace(
        go.Scatter3d(x=F[:, 0],
                    y=F[:, 1],
                    z=F[:, 2],
                    mode='markers',
                    name="Functionals",
                    marker=dict(
                        color='LightSkyBlue',
                        size=4,
                        line=dict(
                            color='LightSkyBlue',
                            width=0.5
                        )
                    ),
                    showlegend=True)
    )
    pio.write_html(iteractivefig, file=paretoFrontFigPath+".html", auto_open=False)

    from pymoo.indicators.igd_plus import IGDPlus

    metric = IGDPlus(pf_a, zero_to_one=True)

    igd = [metric.do(_F) for _F in hist_F]

    fig = plt.figure(figsize = (10, 7))
    plt.plot(n_evals, igd,  color='black', lw=0.7, label="Avg. CV of Pop")
    plt.scatter(n_evals, igd,  facecolor="none", edgecolor='black', marker="p")
    plt.axhline(10**-2, color="red", label="10^-2", linestyle="--")
    plt.title("Convergence")
    plt.xlabel("Function Evaluations")
    plt.ylabel("IGD+")
    plt.yscale("log")
    plt.legend()
    titleConv = model+"_"+wolfKind+"_"+potential+"_Box_"+box+"_Convergence"
    prefix = os.path.split(path)
    convFigPath = os.path.join(prefix[0], titleConv)
    plt.savefig(convFigPath)

    # if you use MO 1.0
    #weights = np.array([0.5, 0.5])
    #weights = np.array([0.333, 0.333, 0.333])
    weights = np.array([0.5, 0.25, 0.25])



    from pymoo.decomposition.asf import ASF

    decomp = ASF()

    i = decomp.do(nF, 1/weights).argmin()

    """
    print("Best regarding ASF: Point \ni = %s\nF = %s" % (i, F[i]))
    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(F[i, 0], F[i, 1], marker="x", color="red", s=200)
    plt.title("Objective Space")
    plt.show()
    """
    x_opts, y_opts = zip(X[i])
    print(x_opts)
    print(y_opts)

    from pymoo.mcdm.pseudo_weights import PseudoWeights

    i = PseudoWeights(weights).do(nF)

    print("Best regarding Pseudo Weights: Point \ni = %s\nF = %s" % (i, F[i]))
    print(X[i])
    x_popts, y_popts = zip(X[i])

    """
    plt.figure(figsize=(7, 5))
    plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
    plt.scatter(F[i, 0], F[i, 1], marker="x", color="red", s=200)
    plt.title("Objective Space")
    plt.show()
    """
    if(plotSuface):
        title = model+"_"+wolfKind+"_"+potential+"_Box_"+box

        prefix = os.path.split(path)
        plotPath = os.path.join(prefix[0], title)

        xx_forplotting = np.linspace(x.min(), x.max(), 1000)
        yy_forplotting = np.linspace(y.min(), y.max(), 1000)

        X_forplotting, Y_forplotting = np.meshgrid(xx_forplotting, yy_forplotting)
        zs = np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel()))
        Z = zs.reshape(X_forplotting.shape)

        iteractivefig = go.Figure()
        iteractivefig.add_surface(autocolorscale=True, x=X_forplotting, y=Y_forplotting, z=Z)
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
        x_opts, y_opts = zip(*X)
        print(x_opts)
        print(y_opts)
        iteractivefig.add_trace(
            go.Scatter3d(x=x_opts,
                        y=y_opts,
                        z=rect_B_spline.ev(x_opts,y_opts))
                        #,mode='markers',
                        #name=key,
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        #showlegend=True)
        )
        """
        
        #x_opts, y_opts = zip(X[i])
        #print(x_opts)
        #print(y_opts)
        iteractivefig.add_trace(
            go.Scatter3d(x=x_popts,
                        y=y_popts,
                        z=rect_B_spline.ev(x_popts,y_popts),
                        mode='markers',
                        name="M.O. 3 - psuedo",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
        )
        iteractivefig.add_trace(
            go.Scatter3d(x=x_opts,
                        y=y_opts,
                        z=rect_B_spline.ev(x_opts,y_opts),
                        mode='markers',
                        name="M.O. 3",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
        )
        iteractivefig.add_trace(
            go.Scatter3d(x=[14],
                        y=[0.12],
                        z=rect_B_spline.ev(14,0.12),
                        mode='markers',
                        name="Rahbari's choice",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
        )
        iteractivefig.add_trace(
            go.Scatter3d(x=[12.286524514052358],
                        y=[0.1893972369539617],
                        z=rect_B_spline.ev(12.286524514052358,0.1893972369539617),
                        mode='markers',
                        name="Best Known",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
        )
        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    print("GD_rcut",x_popts[0])
    print("GD_alpha",y_popts[0])
    print("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0])
    return (("GD_rcut",x_popts[0]), ("GD_alpha",y_popts[0]), ("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0]))