import os
import subprocess
import shutil
import re
import numpy as np
import pandas as pd
from pymbar import timeseries
from scipy.optimize import minimize, fmin_cobyla, basinhopping, Bounds, shgo
from scipy.interpolate import RectBivariateSpline
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.interpolate as si
from skopt import dump, load, expected_minimum

class Calibrator:
    def __init__(self, gomc, wolf_model, wolf_potential, target_y, initial_x, template_directory, template_control_file_name_str, \
                conffile, forcefield, min, max, num_iters=5):
        self.gomc = gomc
        self.wolf_model = wolf_model
        self.wolf_potential = wolf_potential
        self.target_y = target_y
        self.initial_x = initial_x
        self.template_directory = template_directory
        self.template_control_file_name_str = template_control_file_name_str
        self.conffile = conffile
        self.forcefield = forcefield
        self.min = min
        self.max = max
        self.traj = pd.DataFrame()
        self.iteration = 0
        self.num_iters = num_iters
        self.x = 0
        self.fun = 0
        self.obj_calls = {}

    # >= x0-0.05
    def constr1(self,x):
        return (x) - (self.initial_x - 0.05)
    # <= x0+0.05
    def constr2(self,x):
        return (self.initial_x + 0.05) - (x)
    # >= 0.0
    def constr3(self,x):
        return x
    # <= 0.5
    def constr4(self,x):
        return 0.50 - x
    def copy_template(self, x):
        print('creating template files')

        shutil.copyfile("{}{}.conf".format(self.template_directory,self.template_control_file_name_str), "{}{}.conf".format(self.conffile, self.iteration))
        shutil.copyfile("{}{}".format(self.template_directory,self.forcefield), self.forcefield)
        with open("{}{}.conf".format(self.conffile, self.iteration), "a") as myfile:
            defPotLine = "InitStep\t{zero}\n".format(zero=0)
            myfile.write(defPotLine)
            defPotLine = "Wolf\t{freq}\n".format(freq=True)
            myfile.write(defPotLine)
            defPotLine = "WolfKind\t{freq}\n".format(freq=self.wolf_model)
            myfile.write(defPotLine)
            defPotLine = "WolfPotential\t{freq}\n".format(freq=self.wolf_potential)
            myfile.write(defPotLine)  
            defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=0, val=np.round(x[0], 4))
            myfile.write(defAlphaLine)

    def run_simulation(self):
        run_command = "{} +p{} {}.conf > out_{}.dat".format(
            str(self.gomc),
            str(4),
            "{}{}".format(self.conffile, self.iteration),
            "{}{}".format(self.conffile, self.iteration),
        )
        print('gomc gomc_calibration_run_ensemble run_command = ' + str(run_command))

        exec_run_command = subprocess.Popen(
            run_command, shell=True, stderr=subprocess.STDOUT
        )
        os.waitpid(exec_run_command.pid, 0)  # os.WSTOPPED) # 0)
        # test if the simulation actualy finished before checkin and adding 1 to the equilb counter
        # test if the simulation actualy finished before checkin and adding 1 to the equilb counter

    def extract_target(self, x):
        EnRegex = re.compile("ENER_0")
        print('extract_target')

        steps = []
        energies = []
        densities = []
        with open("out_{}{}.dat".format(self.conffile, self.iteration), 'r', encoding='utf8') as f:
            for line in f:
                if EnRegex.match(line):
                    try:
                        steps.append(float(line.split()[1]))
                        #energies.append(float(line.split()[2]))
                        # Use Total_Elec to avoid underreporting error.
                        energies.append(float(line.split()[7]))

                    except:
                        print(line)
                        print("An exception occurred") 
        nskip = 1
        steps_np = np.array(steps)
        energies_np = np.array(energies)
        t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil = energies_np[t0:]
        A_t_equil_steps = steps_np[t0:]
        #energies_np[:t0] = np.nan
        df=pd.DataFrame({'steps':steps_np, np.round(x, 4) :energies_np})
        if (self.traj.empty):
            self.traj = df
        else:
            self.traj = pd.merge(self.traj, df, on='steps', how='outer')
        return A_t_equil.mean()
        
    def extract_reference_target(self):
        EnRegex = re.compile("ENER_0")
        
        steps = []
        energies = []
        densities = []
        with open("{}/out_{}{}.dat".format(self.template_directory,self.conffile, 0), 'r', encoding='utf8') as f:
            for line in f:
                if EnRegex.match(line):
                    try:
                        steps.append(float(line.split()[1]))
                        #energies.append(float(line.split()[2]))
                        # Use Total_Elec to avoid underreporting error.
                        energies.append(float(line.split()[7]))

                    except:
                        print(line)
                        print("An exception occurred") 
        nskip = 1
        steps_np = np.array(steps)
        energies_np = np.array(energies)
        t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil = energies_np[t0:]
        A_t_equil_steps = steps_np[t0:]
        df=pd.DataFrame({'steps':steps_np, 'Ewald':energies_np})
        if (self.traj.empty):
            self.traj = df
        else:
            self.traj = pd.merge(self.traj, df, on='steps', how='outer')
        return A_t_equil.mean()
        
    def loss(self,y_hat):
        return 100.0*(np.abs((self.target_y-y_hat)/self.target_y))

    # objective function
    def objective(self, x):
        #if (x[0] in self.obj_calls):
        #    return self.obj_calls[x[0]]
        
        Calibrator.copy_template(self,x)
        Calibrator.run_simulation(self)
        y_hat = Calibrator.extract_target(self, x[0])
        #self.obj_calls[x[0]] = Calibrator.loss(self,y_hat)
        print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}'.format(self.iteration, x[0], self.target_y, y_hat, Calibrator.loss(self,y_hat)))
        with open("calibration.log", "a") as myfile:
            line = "{0: 3.6f}   {1: 3.6f}\n".format(x[0], Calibrator.loss(self,y_hat))
            myfile.write(line)

        self.iteration = self.iteration + 1
        print("returning")
        return Calibrator.loss(self,y_hat)
    
    def calibrate(self):
        f = lambda x: Calibrator.objective(self,x)
        c1 = lambda x: Calibrator.constr1(self,x)
        c2 = lambda x: Calibrator.constr2(self,x)
        c3 = lambda x: Calibrator.constr3(self,x)
        c4 = lambda x: Calibrator.constr4(self,x)

        cons = [{'type':'ineq', 'fun':c1},        
        {'type':'ineq', 'fun':c2},
        {'type':'ineq', 'fun':c3},        
        {'type':'ineq', 'fun':c4}]


        x0 = np.array([self.initial_x], dtype=np.double)
        res = minimize(f, x0, method='COBYLA', constraints=cons, options={'rhobeg': 0.0025, 'disp': True, 'tol': 0.00125,'catol': 0.000,'maxiter': self.num_iters})
        print(res)
        self.x = np.round(res.x, 4)
        print(self.x)
        alpha = pd.DataFrame()
        alpha["alpha"]=self.x
        alpha.to_csv("best_alpha.csv", header=True)
        print(alpha)
        Calibrator.extract_reference_target(self)
        self.traj.to_csv('alpha_v_mc_steps.csv', header=True, index='steps', sep=' ')
        plot = self.traj.plot(figsize=(10,5), grid=True, x='steps')
        fig = plot.get_figure()
        fig.savefig('plot_alphas.png', bbox_inches='tight')

    
    def calibrate_shgo(self):

        f = lambda x: Calibrator.objective(self,x)

        bounds = Bounds(0.0, 0.5)

        x0 = np.array([self.initial_x], dtype=np.double)
        res = shgo(f, bounds)
        print(res)
        self.x = np.round(res.x, 4)
        print(self.x)
        alpha = pd.DataFrame()
        alpha["alpha"]=self.x
        alpha.to_csv("best_alpha.csv", header=True)
        print(alpha)
        Calibrator.extract_reference_target(self)
        self.traj.to_csv('alpha_v_mc_steps.csv', header=True, index='steps', sep=' ')
        plot = self.traj.plot(figsize=(10,5), grid=True, x='steps')
        fig = plot.get_figure()
        fig.savefig('plot_alphas.png', bbox_inches='tight')


    
    def calibrate_global(self):

        title = "Wolf_Calibration_"+self.wolf_model+"_"+self.wolf_potential+"_BOX_0_"+self.conffile+"0.dat"
        print("Reading ", title)
        df = pd.read_csv(self.template_directory+title,sep='\t',index_col=0, header=None)
        print(df)
        # remove the nan column
        df = df.iloc[: , :-1]
        alphas = np.arange (0.0, 0.5, 0.005)
        df.columns = alphas
        # You can test Vlugt's assertion that 
        """
        For dense liquids such as water and methanol at ambient conditions, 
        it is suﬃcient to plot Figure 1 for a single conﬁguration [93].
        
        dfMean = df.iloc[0, :]
        """

        dfMean = df.mean()
        print(dfMean)
        plot = dfMean.plot(figsize=(10,5), grid=True)
        fig = plot.get_figure()
        fig.savefig('plot_ewald_surface_{}.png'.format(title), bbox_inches='tight')

        sampleEwaldSurface = si.interp1d(dfMean.index, dfMean)
        print(sampleEwaldSurface(0.25))


        f = lambda x: Calibrator.objective(self,x)

        bounds = Bounds(0.0, 0.5)

        x0 = np.array([self.initial_x], dtype=np.double)
        #res = minimize(f, x0, method='COBYLA', constraints=cons, options={'rhobeg': 0.0025, 'disp': True, 'tol': 0.00125,'catol': 0.000,'maxiter': self.num_iters})
        minimizer_kwargs = {"method":"L-BFGS-B", "jac":None, "bounds":bounds, "tol":0.01}
        res = basinhopping(f, x0, minimizer_kwargs=minimizer_kwargs,niter=200,niter_success=40,stepsize=0.05,seed=1)

        print(res)
        self.x = np.round(res.x, 4)
        print(self.x)
        alpha = pd.DataFrame()
        alpha["alpha"]=self.x
        alpha.to_csv("best_alpha.csv", header=True)
        print(alpha)
        Calibrator.extract_reference_target(self)
        self.traj.to_csv('alpha_v_mc_steps.csv', header=True, index='steps', sep=' ')
        plot = self.traj.plot(figsize=(10,5), grid=True, x='steps')
        fig = plot.get_figure()
        fig.savefig('plot_alphas.png', bbox_inches='tight')


    
    def calibrate_skopt(self):
        prefix = "Wolf_Calibration_"+self.wolf_model+"_"+self.wolf_potential+"_BOX_0_"+self.conffile+"0"
        title = "{}.dat".format(prefix)
        print("Reading ", title)
        df = pd.read_csv(self.template_directory+title,sep='\t',index_col=0, header=None)
        print(df)
        # remove the nan column
        df = df.iloc[: , :-1]
        alphas = np.arange (0.0, 0.5, 0.005)
        df.columns = alphas
        # You can test Vlugt's assertion that 
        """
        For dense liquids such as water and methanol at ambient conditions, 
        it is suﬃcient to plot Figure 1 for a single conﬁguration [93].
        
        dfMean = df.iloc[0, :]
        """

        dfMean = df.mean()
        dfMean = dfMean.abs()
        print(dfMean)
        f = lambda x: Calibrator.objective(self,x)

        bounds = Bounds(0.0, 0.5)

        x0 = np.array([self.initial_x], dtype=np.double)

        from skopt import Optimizer
        import matplotlib.pyplot as plt
        from skopt.plots import plot_convergence, plot_evaluations, plot_objective, plot_regret
        opt = Optimizer([(0, 0.5)], "GP", acq_func="EI",
                    acq_optimizer="sampling",
                    initial_point_generator="lhs",
                    n_initial_points=5,
                    random_state=123)
        
        for cal_run in range(self.num_iters):
            if (cal_run == 0):
                x = [self.initial_x]
            else:
                x = opt.ask()
            res = opt.tell(x, f(x))
            dump(res, "./{}_checkpoint.pkl".format(prefix), compress=3)

        print("sampling for expected minimum")
        em = expected_minimum(res, n_random_starts=10000, random_state=123)
        print(em)
        em_l = list(em)
        self.x = np.round(em_l[0][0], 4)
        print(self.x)

        dump(res, '{}.pkl'.format(prefix))

        _ = plot_objective(res, minimum='expected_minimum', n_minimum_search=10000, dimensions=["Alpha"], size=4)
        _.plot(dfMean.index, dfMean.values)
        _.legend(loc="upper left", labels=["Wolf","Wolf Minimum","Ewald"])
        ylim = _.get_ylim()
        print(ylim)
        ylim_l = list(ylim)
        if (ylim_l[0] >= 0.0):
            ylim_l[0] = 0.0
        _.set_ylim(tuple(ylim_l))
        _.set_title("Objective")
        _.set_ylabel("Relative Error")
        plt.savefig('plot_objective.png', bbox_inches='tight')
        plt.show()
        _ = plot_convergence(res)
        _.set_title("Convergence")
        plt.savefig('plot_convergence.png', bbox_inches='tight')
        #plt.show()
        _ = plot_evaluations(res)
        _.set_title("Evaluations")
        plt.savefig('plot_evaluations.png', bbox_inches='tight')
        #plt.show()
        _ = plot_regret(res)
        _.set_title("Regret")
        plt.savefig('plot_regret.png', bbox_inches='tight')
        #plt.show()

        alpha = pd.DataFrame()
        alpha["alpha"]=self.x
        alpha.to_csv("best_alpha.csv", header=True)
        print(alpha)
        Calibrator.extract_reference_target(self)
        self.traj.to_csv('alpha_v_mc_steps.csv', header=True, index='steps', sep=' ')
        plot = self.traj.plot(figsize=(10,5), grid=True, x='steps')
        fig = plot.get_figure()
        fig.savefig('plot_alphas.png', bbox_inches='tight')

"""
    def calibrate(self):
        f = lambda x: Calibrator.objective(self,x)
        c1 = lambda x: Calibrator.constr1(self,x)
        c2 = lambda x: Calibrator.constr2(self,x)
        c3 = lambda x: Calibrator.constr3(self,x)
        c4 = lambda x: Calibrator.constr4(self,x)

        x0 = np.array([self.initial_x], dtype=np.double)
        res = \
            fmin_cobyla(f, 
                    x0, 
                    cons=[],
                    rhobeg=0.01,
                    rhoend=0.0025, 
                    disp=3, 
                    catol=0.0,
                    maxfun=self.num_iters) 
        print(res)
        self.x = res

        Calibrator.extract_reference_target(self)
        self.traj.to_csv('alpha_v_mc_steps.csv', header=True, index='steps', sep=' ')
        plot = self.traj.plot(figsize=(10,5), grid=True, x='steps')
        fig = plot.get_figure()
        fig.savefig('plot_alphas.png', bbox_inches='tight')

"""