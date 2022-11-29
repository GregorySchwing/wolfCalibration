#!/usr/bin/python

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import glob
from matplotlib import cm
from scipy.optimize import minimize, brute, shgo, dual_annealing, differential_evolution, root
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
    #if(wolfKind == "VLUGT" and potential == "DSP"):
    #    z = np.abs(df4.iloc[:,0].to_numpy())
    #elif(wolfKind == "VLUGTWINTRACUTOFF" and potential == "DSP"):
    #    z = np.abs(df4.iloc[:,0].to_numpy())
    #else:
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
    """
    xx_quick_forplotting = np.linspace(x.min(), x.max(), 1000)
    yy_quick_forplotting = np.linspace(y.min(), y.max(), 1000)

    X_quick_forplotting, Y_quick_forplotting = np.meshgrid(xx_quick_forplotting, yy_quick_forplotting)
    #zs = np.exp(np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel())))
    zsquick_ = np.array(rect_B_spline.ev(X_quick_forplotting.ravel(), Y_quick_forplotting.ravel()))
    Zquick = zsquick_.reshape(X_quick_forplotting.shape)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X_quick_forplotting, Y_quick_forplotting, Zquick)
    title = model+"_"+wolfKind+"_"+potential+"_Box_"+box
    ax.set_title(title)

    ax.set_xlabel('Rc')
    ax.set_ylabel('Alpha')
    ax.set_zlabel('Error')

    plt.show()
    """
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
    derivs_wrt_alpha_and_rcut = rect_B_spline.partial_derivative(pd_RCut_varies_alpha_varies[0],pd_RCut_varies_alpha_varies[1])
    tck_wrt_alpha_and_rcut = [derivs_wrt_alpha_and_rcut.tck[0], derivs_wrt_alpha_and_rcut.tck[1],derivs_wrt_alpha_and_rcut.tck[2],derivs_wrt_alpha_and_rcut.degrees[0],derivs_wrt_alpha_and_rcut.degrees[1]]

    foundSoln = False
    tolPower = 0
    from scipy.optimize import NonlinearConstraint, Bounds

    rosen = lambda a : np.abs(rect_B_spline.ev(a[0], a[1]))
    rosenD = lambda a : np.abs(interpolate.bisplev(a[0], a[1], tck_pd))
    nlc = NonlinearConstraint(rosenD, -np.inf, 1.0)
    bounds = [(x.min(), x.max()), (y.min(), y.max())]
    result = differential_evolution(rosen, bounds, updating='deferred',constraints=(nlc),seed=1)
    print(result.x, result.fun)
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
        

        iteractivefig.add_trace(
            go.Scatter3d(x=[result.x[0]],
                        y=[result.x[1]],
                        z=rect_B_spline.ev(result.x[0],result.x[1]),
                        mode='markers',
                        name="M.O. 3 - psuedo",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
        )
        iteractivefig.add_trace(
            go.Scatter3d(x=[result.x[0]],
                        y=[result.x[1]],
                        z=rect_B_spline.ev(result.x[0],result.x[1]),
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

    print("GD_rcut",result.x[0])
    print("GD_alpha",result.x[1])
    print("GD_relerr",rect_B_spline.ev(result.x[0],result.x[1]))
    return (("GD_rcut",result.x[0]), ("GD_alpha",result.x[1]), ("GD_relerr",rect_B_spline.ev(result.x[0],result.x[1])))


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

        result.x[0], result.x[1] = zip(X)
        #result.x[0], result.x[1] = zip(X[i])
        

        print(result.x[0])
        print(result.x[1])
        from pymoo.mcdm.pseudo_weights import PseudoWeights

        i = PseudoWeights(weights).do(nF)

        print("Best regarding Pseudo Weights: Point \ni = %s\nF = %s" % (i, F[i]))
        print(X[i])
        #result.x, y_popts = zip(X[i])
        result.x, y_popts = zip(X)
        print(result.x)
        print(y_popts)
        """
        if(plotSuface):


            xx_forplotting = np.linspace(x.min(), x.max(), 100)
            yy_forplotting = np.linspace(y.min(), y.max(), 100)

            X_forplotting, Y_forplotting = np.meshgrid(xx_forplotting, yy_forplotting)
            #zs = np.exp(np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel())))
            zs = np.array(rect_B_spline.ev(X_forplotting.ravel(), Y_forplotting.ravel()))
            Z = zs.reshape(X_forplotting.shape)
            """
            ZeroSlice = np.where(((Z < 0.0001) or (Z > 0.0001)), Z)

            plt.scatter(ZeroSlice)
            plt.savefig("test")
            """
            iteractivefig.add_surface(autocolorscale=True, x=X_forplotting, y=Y_forplotting, z=Z)


            """
            #result.x[0], result.x[1] = zip(X[i])
            #print(result.x[0])
            #print(result.x[1])
            iteractivefig.add_trace(
            go.Scatter3d(x=result.x,
                        y=y_popts,
                        z=rect_B_spline.ev(result.x,y_popts),
                        mode='markers',
                        name="M.O. 3 - psuedo",
                        #hovertext=["REF"] if len(xvals) == 1 else [str(x) for x in scales],
                        showlegend=True)
            )
            iteractivefig.add_trace(
                go.Scatter3d(x=[result.x[0]],
                            y=[result.x[1]],
                            z=rect_B_spline.ev(result.x[0],result.x[1]),
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
    #plt.colorbar()
    #plt.savefig(plotPath)
    # Using any of the single point BF/GD methods is obviously a bad idea.
    #    return (("BF_rcut",bfXY[0]), ("BF_alpha",bfXY[1]), ("BF_relerr",ZBF), ("GD_rcut",gdXY[0]), ("GD_alpha",gdXY[1]), ("GD_relerr",ZGD), ("GD_jac_rcut",gdJacXY[0]), ("GD_jac_alpha",gdJacXY[1]))
    # The question is which of the above optimizations to use.  For now, I am going with "REF" AUC as the metric.

    #print("GD_rcut",result.x[0])
    #print("GD_alpha",result.x[1])
    #print("GD_relerr",rect_B_spline.ev(result.x,y_popts)[0])
    #return (("GD_rcut",result.x[0]), ("GD_alpha",result.x[1]), ("GD_relerr",rect_B_spline.ev(result.x,y_popts)[0]))
