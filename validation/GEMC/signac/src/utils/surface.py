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
import math
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

def _find_global_function_max(self):
    print("finding largest function val")
    # Single objective lambda function (ZError).
    f_max = lambda x: -np.abs(self.rect_B_spline.ev(x[0], x[1]))
    #f_min = lambda x: np.abs(self.rect_B_spline.ev(x[0], x[1]))
    bounds = [(self.RCutMin, self.RCutMax),(self.AlphaMin, self.AlphaMax)]
    sptdifferential_evolutionOutNox0 = differential_evolution(f_max, tol=0,bounds=bounds)
    print("Calling differential_evolutionX0")
    print(sptdifferential_evolutionOutNox0) 
    return sptdifferential_evolutionOutNox0.x

from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize
class DEProblem(ElementwiseProblem):

    def __init__(self, rect_B_spline, RCutMin, RCutMax, AlphaMin, AlphaMax):
        super().__init__(n_var=2,
                         n_obj=1,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline

        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.tolerance = 0.01


    def _evaluate(self, x, out, *args, **kwargs):
        # Minimize Relative Error
        f1 = -np.abs(self.rect_B_spline.ev(x[0], x[1]))
        g1 = np.abs(self.rect_B_spline.ev(x[0], x[1])) - self.tolerance

        out["F"] = [f1]
        out["G"] = [g1]

class DEProblemDeriv(ElementwiseProblem):

    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax):
        super().__init__(n_var=2,
                         n_obj=1,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.tck = tck_pd
        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.tolerance = 0.01

    def _evaluate(self, x, out, *args, **kwargs):
        # Maximize Relative Error
        f1 = -np.abs(interpolate.bisplev(x[0], x[1], self.tck))
        g1 = np.abs(self.rect_B_spline.ev(x[0], x[1])) - self.tolerance

        out["F"] = [f1]
        out["G"] = [g1]

class MyProblem(ElementwiseProblem):


    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax, LowerBoundRcut):
        super().__init__(n_var=2,
                         n_obj=3,
                         n_ieq_constr=2,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.LowerBoundRcut = LowerBoundRcut

        self.tck_pd = tck_pd


        self.pd_RCut_constant_alpha_constant = [0,0]
        self.pd_RCut_varies_alpha_constant = [1,0]
        self.pd_RCut_varies_alpha_constant_DD = [2,0]
        self.pd_RCut_constant_alpha_varies = [0,1]
        self.pd_RCut_constant_alpha_varies_DD = [0,2]
        self.pd_RCut_varies_alpha_varies = [1,1]
        self.pd_RCut_varies_alpha_varies_DD = [2,2]
        
        self.derivs_wrt_alpha = rect_B_spline.partial_derivative(self.pd_RCut_constant_alpha_varies[0],self.pd_RCut_constant_alpha_varies[1])
        self.derivs_wrt_alpha_DD = rect_B_spline.partial_derivative(self.pd_RCut_constant_alpha_varies_DD[0],self.pd_RCut_constant_alpha_varies_DD[1])
        self.derivs_wrt_rcut = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_constant[0],self.pd_RCut_varies_alpha_constant[1])
        self.derivs_wrt_rcut_DD = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_constant_DD[0],self.pd_RCut_varies_alpha_constant_DD[1])
        self.derivs_wrt_alpha_and_rcut = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_varies[0],self.pd_RCut_varies_alpha_varies[1])
        self.derivs_wrt_alpha_and_rcut_DD = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_varies_DD[0],self.pd_RCut_varies_alpha_varies_DD[1])

        self.tck_wrt_alpha = [self.derivs_wrt_alpha.tck[0], self.derivs_wrt_alpha.tck[1],self.derivs_wrt_alpha.tck[2],self.derivs_wrt_alpha.degrees[0],self.derivs_wrt_alpha.degrees[1]]
        self.tck_wrt_alpha_DD = [self.derivs_wrt_alpha_DD.tck[0], self.derivs_wrt_alpha_DD.tck[1],self.derivs_wrt_alpha_DD.tck[2],self.derivs_wrt_alpha_DD.degrees[0],self.derivs_wrt_alpha_DD.degrees[1]]
        self.tck_wrt_rcut = [self.derivs_wrt_rcut.tck[0], self.derivs_wrt_rcut.tck[1],self.derivs_wrt_rcut.tck[2],self.derivs_wrt_rcut.degrees[0],self.derivs_wrt_rcut.degrees[1]]
        self.tck_wrt_rcut_DD = [self.derivs_wrt_rcut_DD.tck[0], self.derivs_wrt_rcut_DD.tck[1],self.derivs_wrt_rcut_DD.tck[2],self.derivs_wrt_rcut_DD.degrees[0],self.derivs_wrt_rcut_DD.degrees[1]]
        self.tck_wrt_alpha_and_rcut = [self.derivs_wrt_alpha_and_rcut.tck[0], self.derivs_wrt_alpha_and_rcut.tck[1],self.derivs_wrt_alpha_and_rcut.tck[2],self.derivs_wrt_alpha_and_rcut.degrees[0],self.derivs_wrt_alpha_and_rcut.degrees[1]]
        self.tck_wrt_alpha_and_rcut_DD = [self.derivs_wrt_alpha_and_rcut_DD.tck[0], self.derivs_wrt_alpha_and_rcut_DD.tck[1],self.derivs_wrt_alpha_and_rcut_DD.tck[2],self.derivs_wrt_alpha_and_rcut_DD.degrees[0],self.derivs_wrt_alpha_and_rcut_DD.degrees[1]]
        

        self.MyDEProblem = DEProblem(rect_B_spline, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.MyDEProblem,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.FMax = np.abs(res.F[0])

        self.DEProblemDerivWRTRcut = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_rcut, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRTRcut,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRTRcut_max = np.abs(res.F[0])

        self.DEProblemDerivWRTRcut_DD = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_rcut_DD, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRTRcut_DD,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRTRcut_DD_max = np.abs(res.F[0])   

        self.DEProblemDerivWRTAlpha = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_alpha, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRTAlpha,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRTAlpha_max = np.abs(res.F[0])
        
        self.DEProblemDerivWRTAlpha_DD = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_alpha_DD, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRTAlpha_DD,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRTAlpha_DD_max = np.abs(res.F[0])        

        self.DEProblemDerivWRT_RCut_and_Alpha = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_alpha_and_rcut, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRT_RCut_and_Alpha,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRT_RCut_and_Alpha_max = np.abs(res.F[0])
        
        self.DEProblemDerivWRT_RCut_and_Alpha_DD = DEProblemDeriv(self.rect_B_spline, self.tck_wrt_alpha_and_rcut_DD, RCutMin, RCutMax, AlphaMin, AlphaMax)
        
        algorithm = DE(
            pop_size=100,
            sampling=LHS(),
            variant="DE/rand/1/bin",
            CR=0.3,
            dither="vector",
            jitter=False
        )

        res = minimize(self.DEProblemDerivWRT_RCut_and_Alpha_DD,
                    algorithm,
                    seed=1,
                    verbose=True)

        print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
        
        self.DEProblemDerivWRT_RCut_and_Alpha_DD_max = np.abs(res.F[0])
        


                

    def neg_bspline( x ):
        global tck
        f = -interpolate.bisplev( x[0], x[1], tck, dx=0, dy=0)
        g = [-interpolate.bisplev( x[0], x[1], tck, dx=1, dy=0 ), -interpolate.bisplev( x[0], x[1], tck, dx=0, dy=1)]
        return f, g

    def _find_global_max_deriv(self, tck, label):
        print("finding largest ", label)
        # Single objective lambda function (ZError).
        f = lambda x : -np.abs(interpolate.bisplev(x[0], x[1], tck))

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


class MyDumProblem(ElementwiseProblem):

    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax, tolerance_power):
        super().__init__(n_var=2,
                         n_obj=1,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.tck_pd = tck_pd
        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.F1Min = 0
        self.F2Min = 0
        self.F3Min = 0
        self.tolerance = pow(10, -tolerance_power)
        #self.tolerance = pow(10, -2)


        self.pd_RCut_constant_alpha_constant = [0,0]
        self.pd_RCut_varies_alpha_constant = [1,0]
        self.pd_RCut_varies_alpha_constant_DD = [2,0]
        self.pd_RCut_constant_alpha_varies = [0,1]
        self.pd_RCut_constant_alpha_varies_DD = [0,2]
        self.pd_RCut_varies_alpha_varies = [1,1]
        self.pd_RCut_varies_alpha_varies_DD = [2,2]
        
        self.derivs_wrt_alpha_and_rcut = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_varies[0],self.pd_RCut_varies_alpha_varies[1])
        self.derivs_wrt_alpha_and_rcut_DD = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_varies_DD[0],self.pd_RCut_varies_alpha_varies_DD[1])

        self.tck_wrt_alpha_and_rcut = [self.derivs_wrt_alpha_and_rcut.tck[0], self.derivs_wrt_alpha_and_rcut.tck[1],self.derivs_wrt_alpha_and_rcut.tck[2],self.derivs_wrt_alpha_and_rcut.degrees[0],self.derivs_wrt_alpha_and_rcut.degrees[1]]
        self.tck_wrt_alpha_and_rcut_DD = [self.derivs_wrt_alpha_and_rcut_DD.tck[0], self.derivs_wrt_alpha_and_rcut_DD.tck[1],self.derivs_wrt_alpha_and_rcut_DD.tck[2],self.derivs_wrt_alpha_and_rcut_DD.degrees[0],self.derivs_wrt_alpha_and_rcut_DD.degrees[1]]

    def _evaluate(self, x, out, *args, **kwargs):
        f4 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut)))
        g2 = (np.abs(self.rect_B_spline.ev(x[0], x[1]))) - self.tolerance
        out["F"] = [f4]
        #out["F"] = [f1, f2]
        out["G"] = [g2]


class MyProblemNorm(ElementwiseProblem):

    def __init__(self, rect_B_spline, tck_pd, RCutMin, RCutMax, AlphaMin, AlphaMax, FMax, DEProblemDerivWRTRcut_max, DEProblemDerivWRTAlpha_max, DEProblemDerivWRT_RCut_and_Alpha_max, tolerance_power = 4):
        super().__init__(n_var=2,
                         n_obj=2,
                         n_ieq_constr=1,
                         xl=np.array([RCutMin,AlphaMin]),
                         xu=np.array([RCutMax,AlphaMax]))
        self.rect_B_spline = rect_B_spline
        self.tck_pd = tck_pd
        self.RCutMin = RCutMin
        self.RCutMax = RCutMax
        self.AlphaMin = AlphaMin
        self.AlphaMax = AlphaMax
        self.F1Min = 0
        self.F2Min = 0
        self.F3Min = 0
        self.FMax = FMax
        self.DEProblemDerivWRTRcut_max = DEProblemDerivWRTRcut_max
        self.DEProblemDerivWRTAlpha_max = DEProblemDerivWRTAlpha_max
        self.DEProblemDerivWRT_RCut_and_Alpha_max = DEProblemDerivWRT_RCut_and_Alpha_max
        self.tolerance = pow(10, -tolerance_power)
        
        self.pd_RCut_constant_alpha_constant = [0,0]
        self.pd_RCut_varies_alpha_constant = [1,0]
        self.pd_RCut_constant_alpha_varies = [0,1]
        self.pd_RCut_varies_alpha_varies = [1,1]

        self.derivs_wrt_alpha = rect_B_spline.partial_derivative(self.pd_RCut_constant_alpha_varies[0],self.pd_RCut_constant_alpha_varies[1])
        self.derivs_wrt_rcut = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_constant[0],self.pd_RCut_varies_alpha_constant[1])
        self.derivs_wrt_alpha_and_rcut = rect_B_spline.partial_derivative(self.pd_RCut_varies_alpha_varies[0],self.pd_RCut_varies_alpha_varies[1])

        self.tck_wrt_alpha = [self.derivs_wrt_alpha.tck[0], self.derivs_wrt_alpha.tck[1],self.derivs_wrt_alpha.tck[2],self.derivs_wrt_alpha.degrees[0],self.derivs_wrt_alpha.degrees[1]]
        self.tck_wrt_rcut = [self.derivs_wrt_rcut.tck[0], self.derivs_wrt_rcut.tck[1],self.derivs_wrt_rcut.tck[2],self.derivs_wrt_rcut.degrees[0],self.derivs_wrt_rcut.degrees[1]]
        self.tck_wrt_alpha_and_rcut = [self.derivs_wrt_alpha_and_rcut.tck[0], self.derivs_wrt_alpha_and_rcut.tck[1],self.derivs_wrt_alpha_and_rcut.tck[2],self.derivs_wrt_alpha_and_rcut.degrees[0],self.derivs_wrt_alpha_and_rcut.degrees[1]]


    def _evaluate(self, x, out, *args, **kwargs):
        # Minimize Relative Error
        f1 = (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax)*(np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))/self.DEProblemDerivWRT_RCut_and_Alpha_max)
        # Minimize RCut
        f2 = (x[0]-self.RCutMin)/(self.RCutMax-self.RCutMin)
        g1 = 0
        """
        # Minimize dWRT RCut
        f3 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_rcut))/self.DEProblemDerivWRTRcut_max
        # Minimize dWRT Alpha
        f4 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha))/self.DEProblemDerivWRTAlpha_max
        # Minimize Gradient
        f5 = np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))/self.DEProblemDerivWRT_RCut_and_Alpha_max
        # This way I don't wind up on the side of a hill
        # Gradient <= RCutLoss
        g1 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_rcut))/self.DEProblemDerivWRTRcut_max) - (x[0]-self.RCutMin)/(self.RCutMax-self.RCutMin)
        g2 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha))/self.DEProblemDerivWRTAlpha_max) - (x[0]-self.RCutMin)/(self.RCutMax-self.RCutMin)
        g3 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))/self.DEProblemDerivWRT_RCut_and_Alpha_max) - (x[0]-self.RCutMin)/(self.RCutMax-self.RCutMin)
        # Gradient <= RelError
        #g4 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_rcut))/self.DEProblemDerivWRTRcut_max) - (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax)
        # Gradient <= RelError
        #g5 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha))/self.DEProblemDerivWRTAlpha_max) - (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax)
        # Gradient <= RelError
        #g6 = (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))/self.DEProblemDerivWRT_RCut_and_Alpha_max) - (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax)
        # RelError <= D_wrt_rcut
        g4 = (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax) - (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_rcut))/self.DEProblemDerivWRTRcut_max)
        # RelError <= D_wrt_alpha
        g5 = (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax) - (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha))/self.DEProblemDerivWRTAlpha_max)
        # RelError <= Gradient
        g6 = (np.abs(self.rect_B_spline.ev(x[0], x[1]))/self.FMax) - (np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))/self.DEProblemDerivWRT_RCut_and_Alpha_max)
        g7 = np.abs(self.rect_B_spline.ev(x[0], x[1]))-self.tolerance
        """
        #out["F"] = [f1, f2, f3, f4, f5]
        out["F"] = [f1, f2]
        out["G"] = [g1]
        #out["G"] = [g1, g2, g3, g4, g5, g6, g7]
        #out["G"] = [g1, g3, g5]


    def _calc_pareto_front(self, flatten=True, *args, **kwargs):

        num_pts = 100
        rcuts = np.linspace(self.RCutMin, self.RCutMax, num_pts)
        alphas = np.linspace(self.AlphaMin, self.AlphaMax, num_pts)
        rcuts_g, alphas_g = np.meshgrid(rcuts, alphas)

        points = np.array(list(zip(rcuts_g.ravel(), alphas_g.ravel())))


        F1_a_costs = np.array((np.abs(self.rect_B_spline.ev(rcuts_g.ravel(), alphas_g.ravel())))/self.FMax)
        F2_a_costs = ((rcuts_g.ravel()-self.RCutMin)/(self.RCutMax-self.RCutMin))
        deriv1Lambda = lambda x : np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_rcut))
        deriv2Lambda = lambda x : np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha))
        deriv3Lambda = lambda x : np.abs(interpolate.bisplev(x[0], x[1], self.tck_wrt_alpha_and_rcut))

        F3_a_costs = np.array([deriv1Lambda(i) for i in points])
        F3_a_costs = F3_a_costs/self.DEProblemDerivWRTRcut_max

        F4_a_costs = np.array([deriv2Lambda(i) for i in points])
        F4_a_costs = F4_a_costs/self.DEProblemDerivWRTAlpha_max

        F5_a_costs = np.array([deriv2Lambda(i) for i in points])
        F5_a_costs = F5_a_costs/self.DEProblemDerivWRT_RCut_and_Alpha_max

        costs = np.array(list(zip(F1_a_costs, F2_a_costs, F3_a_costs, F4_a_costs, F5_a_costs)))
        pass2Method = np.array(list(zip(costs,points)))
        boolean_array = is_pareto_efficient(costs)
        print(boolean_array)
        #pf_a = pareto_frontier_multi(costs)
        return costs[boolean_array]
        
def neg_bspline( x ):
    global tck
    f = -interpolate.bisplev( x[0], x[1], tck, dx=0, dy=0)
    g = [-interpolate.bisplev( x[0], x[1], tck, dx=1, dy=0 ), -interpolate.bisplev( x[0], x[1], tck, dx=0, dy=1)]
    return f, g

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

    # Spline and pymoo prefer convex surfaces.
    # All are convex at zero intercept when nonabs except for VDSP and VWICDSP
    if(wolfKind == "VLUGT" and potential == "DSP"):
        z = np.abs(df4.iloc[:,0].to_numpy())
    elif(wolfKind == "VLUGTWINTRACUTOFF" and potential == "DSP"):
        z = np.abs(df4.iloc[:,0].to_numpy())
    else:
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
    
    RCutMin = x.min()
    val = interpolate.bisplev(exampleX2, exampleY2, tck_pd)
    print("deriv", exampleX2, exampleY2)
    print(val)
    
    from pymoo.termination import get_termination
    # Create problem to get the unnormalized Pareto Front
    # Since I only use gradient, it's a single objective, no need to scale.
    #problemUnNorm = MyProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 10)
    #pf_for_norm = problemUnNorm.pareto_front(use_cache=False, flatten=False)

    tolPower = 0
    for tolPow in range(8, -1, -1):
        print("Tolerance = ", pow(10, -tolPow))
        tolPower = tolPow
        #problem = MyProblemNorm(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), problemUnNorm.FMax, problemUnNorm.DEProblemDerivWRTRcut_max, problemUnNorm.DEProblemDerivWRTAlpha_max, problemUnNorm.DEProblemDerivWRT_RCut_and_Alpha_max, tolPower)
        
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.operators.crossover.sbx import SBX
        from pymoo.operators.mutation.pm import PM
        from pymoo.operators.sampling.rnd import FloatRandomSampling

        algorithm = NSGA2(
            pop_size=400,
            n_offsprings=100,
            sampling=FloatRandomSampling(),
            crossover=SBX(prob=0.9, eta=15),
            mutation=PM(eta=20),
            eliminate_duplicates=True
        )
        termination = get_termination("n_gen", 400)


        from pymoo.optimize import minimize
        prob = MyDumProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), tolPower)

        res = minimize(prob,
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

        k = 0
        try:
            k = np.where(np.array(hist_cv) <= 0.0)[0].min()
        except:
            print(f"No feasible solution in Generation {k} after {n_evals[k]} evaluations.")
            print(f"Increase tolerance.")
            continue
        print(f"At least one feasible solution in Generation {k} after {n_evals[k]} evaluations.")


        # replace this line by `hist_cv` if you like to analyze the least feasible optimal solution and not the population
        vals = hist_cv_avg

        try:
            k = np.where(np.array(vals) <= 0.0)[0].min()
            print(f"Whole population feasible in Generation {k} after {n_evals[k]} evaluations.")
            break
        except:
            print(f"Whole population not feasible after {n_evals[k]} evaluations.")
            print(f"Increase tolerance.")
            continue


    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)

    from pymoo.indicators.hv import Hypervolume

    """

    metric = Hypervolume(ref_point= np.array([0.5, 0.5, 0.5, 0.5, 0.5]),
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
    """

    #xl, xu = problem.bounds()

    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)

    nF = (F - approx_ideal) / (approx_nadir - approx_ideal)

    fl = nF.min(axis=0)
    fu = nF.max(axis=0)
    #print(f"Scale f1: [{fl[0]}, {fu[0]}]")
    #print(f"Scale f2: [{fl[1]}, {fu[1]}]")


    #pf_a = problem.pareto_front(use_cache=False, flatten=False)
    """
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

    """
    """
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
    """
    # if you use MO 1.0
    weights = np.array([1.0])
    #weights = np.array([0.5,0.5])
    #weights = np.array([0.333, 0.333, 0.333])
    #weights = np.array([0.2, 0.2, 0.2, 0.2, 0.2])



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
    #x_opts, y_opts = zip(X)
    #x_opts, y_opts = zip(X[i])
    if (np.isscalar(X[i])):
        x_opts, y_opts = zip(X)
    else:
        x_opts, y_opts = zip(X[i])    

    #print(x_opts)
    #print(y_opts)
    from pymoo.mcdm.pseudo_weights import PseudoWeights

    i = PseudoWeights(weights).do(nF)

    print("Best regarding Pseudo Weights: Point \ni = %s\nF = %s" % (i, F[i]))
    if (np.isscalar(X[i])):
        x_popts, y_popts = zip(X)
    else:
        x_popts, y_popts = zip(X[i])

    #x_popts, y_popts = zip(X)
    print(x_popts)
    print(y_popts)

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
                    margin=dict(r=20, b=10, l=10, t=10),
                    font=dict(size=18))
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
        """
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
        """
        pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    print("GD_rcut",x_popts[0])
    print("GD_alpha",y_popts[0])
    print("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0])
    return (("GD_rcut",x_popts[0]), ("GD_alpha",y_popts[0]), ("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0]))


def plot_all_surfaces(pr, job, file, model, wolfKind, potential, box, plotSuface=False):
    title = model+"_"+wolfKind+"_"+potential+"_Box_"+box+"_all_surfaces"
    plotPath = job.fn(title)
    iteractivefig = go.Figure()

    ewald_sp = job.statepoint()
    ewald_sp['replica_number_int']=0
    ewald_sp['electrostatic_method']="Wolf"
    ewald_sp['wolf_potential']="Calibrator"
    ewald_sp['wolf_model']="Calibrator"
    jobs = list(pr.find_jobs({"sp.replica_number_int": 0, "sp.electrostatic_method": "Wolf", "sp.wolf_potential": "Calibrator", "sp.wolf_model": "Calibrator"}))
    for ewald_job in jobs:
        df = pd.read_csv(ewald_job.fn(file),sep='\t',index_col=0)
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

        # Spline and pymoo prefer convex surfaces.
        # All are convex at zero intercept when nonabs except for VDSP and VWICDSP
        if(wolfKind == "VLUGT" and potential == "DSP"):
            z = np.abs(df4.iloc[:,0].to_numpy())
        elif(wolfKind == "VLUGTWINTRACUTOFF" and potential == "DSP"):
            z = np.abs(df4.iloc[:,0].to_numpy())
        else:
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
        
        RCutMin = x.min()
        val = interpolate.bisplev(exampleX2, exampleY2, tck_pd)
        print("deriv", exampleX2, exampleY2)
        print(val)
        
        from pymoo.termination import get_termination
        # Create problem to get the unnormalized Pareto Front
        # Since I only use gradient, it's a single objective, no need to scale.
        #problemUnNorm = MyProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 10)
        #pf_for_norm = problemUnNorm.pareto_front(use_cache=False, flatten=False)
        """
        tolPower = 0
        for tolPow in range(8, -1, -1):
            print("Tolerance = ", pow(10, -tolPow))
            tolPower = tolPow
            #problem = MyProblemNorm(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), problemUnNorm.FMax, problemUnNorm.DEProblemDerivWRTRcut_max, problemUnNorm.DEProblemDerivWRTAlpha_max, problemUnNorm.DEProblemDerivWRT_RCut_and_Alpha_max, tolPower)
            
            from pymoo.algorithms.moo.nsga2 import NSGA2
            from pymoo.operators.crossover.sbx import SBX
            from pymoo.operators.mutation.pm import PM
            from pymoo.operators.sampling.rnd import FloatRandomSampling

            algorithm = NSGA2(
                pop_size=400,
                n_offsprings=100,
                sampling=FloatRandomSampling(),
                crossover=SBX(prob=0.9, eta=15),
                mutation=PM(eta=20),
                eliminate_duplicates=True
            )
            termination = get_termination("n_gen", 400)


            from pymoo.optimize import minimize
            #prob = MyDumProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), problemUnNorm.FMax, problemUnNorm.DEProblemDerivWRTRcut_max, problemUnNorm.DEProblemDerivWRTRcut_DD_max, problemUnNorm.DEProblemDerivWRTAlpha_max, problemUnNorm.DEProblemDerivWRTAlpha_DD_max, problemUnNorm.DEProblemDerivWRT_RCut_and_Alpha_max, problemUnNorm.DEProblemDerivWRT_RCut_and_Alpha_DD_max, tolPower)
            prob = MyDumProblem(rect_B_spline, tck_pd, x.min(), x.max(), y.min(), y.max(), 0, 0, 0, 0, 0, 0, 0, tolPower)

            res = minimize(prob,
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

            k = 0
            try:
                k = np.where(np.array(hist_cv) <= 0.0)[0].min()
            except:
                print(f"No feasible solution in Generation {k} after {n_evals[k]} evaluations.")
                print(f"Increase tolerance.")
                continue
            print(f"At least one feasible solution in Generation {k} after {n_evals[k]} evaluations.")


            # replace this line by `hist_cv` if you like to analyze the least feasible optimal solution and not the population
            vals = hist_cv_avg

            try:
                k = np.where(np.array(vals) <= 0.0)[0].min()
                print(f"Whole population feasible in Generation {k} after {n_evals[k]} evaluations.")
                break
            except:
                print(f"Whole population not feasible after {n_evals[k]} evaluations.")
                print(f"Increase tolerance.")
                continue


        approx_ideal = F.min(axis=0)
        approx_nadir = F.max(axis=0)

        from pymoo.indicators.hv import Hypervolume

        xl, xu = prob.bounds()

        approx_ideal = F.min(axis=0)
        approx_nadir = F.max(axis=0)

        nF = (F - approx_ideal) / (approx_nadir - approx_ideal)

        fl = nF.min(axis=0)
        fu = nF.max(axis=0)
        #print(f"Scale f1: [{fl[0]}, {fu[0]}]")
        #print(f"Scale f2: [{fl[1]}, {fu[1]}]")


        #pf_a = problem.pareto_front(use_cache=False, flatten=False)

        # if you use MO 1.0
        weights = np.array([0.5,0.5])
        #weights = np.array([0.333, 0.333, 0.333])
        #weights = np.array([0.2, 0.2, 0.2, 0.2, 0.2])



        from pymoo.decomposition.asf import ASF

        decomp = ASF()

        i = decomp.do(nF, 1/weights).argmin()

        x_opts, y_opts = zip(X)
        #x_opts, y_opts = zip(X[i])
        

        print(x_opts)
        print(y_opts)
        from pymoo.mcdm.pseudo_weights import PseudoWeights

        i = PseudoWeights(weights).do(nF)

        print("Best regarding Pseudo Weights: Point \ni = %s\nF = %s" % (i, F[i]))
        print(X[i])
        #x_popts, y_popts = zip(X[i])
        x_popts, y_popts = zip(X)
        print(x_popts)
        print(y_popts)
        """
        if(plotSuface):


            xx_forplotting = np.linspace(x.min(), x.max(), 100)
            yy_forplotting = np.linspace(y.min(), y.max(), 100)

            X_forplotting, Y_forplotting = np.meshgrid(xx_forplotting, yy_forplotting)
            #zs = np.exp(np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel())))
            zs = np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel()))
            Z = zs.reshape(X_forplotting.shape)

            from scipy import optimize
            nx, ny = 100, 100
            guess = np.zeros((nx, ny), float)
            fR = lambda x: rect_B_spline.ev(x[0], x[1])
            print("about to call root")
            sol = optimize.root(fR, guess, method='krylov')
            print(sol.x)

            import matplotlib.pyplot as plt
            plt.pcolormesh(xx_forplotting, yy_forplotting, sol.x, shading='gouraud')
            plt.colorbar()
            plt.savefig(plotPath)
            iteractivefig.add_surface(autocolorscale=True, x=X_forplotting, y=Y_forplotting, z=Z)


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
            """
    """
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
    pio.write_html(iteractivefig, file=plotPath+".html", auto_open=False)

    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    #print("GD_rcut",x_popts[0])
    #print("GD_alpha",y_popts[0])
    #print("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0])
    #return (("GD_rcut",x_popts[0]), ("GD_alpha",y_popts[0]), ("GD_relerr",rect_B_spline.ev(x_popts,y_popts)[0]))
