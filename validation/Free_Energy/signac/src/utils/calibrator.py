    import os
    import subprocess
    import shutil
    import re
    import numpy as np
    import pandas as pd
    from pymbar import timeseries
    from scipy.optimize import fmin_l_bfgs_b, minimize, fmin_cobyla
    from scipy.optimize import least_squares
    import matplotlib.pyplot as plt

    class Calibrator:
        def __init__(self, gomc, wolf_model, wolf_potential, target_y, initial_x, template_directory, template_control_file_name_str, \
                    conffile, forcefield, min, max, num_iters=10):
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

        def constr(self, x):
            return (float((self.initial_x - 0.05) < x)) + float(( x > (self.initial_x + 0.05)))

        def copy_template(self, x):
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
                defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=0, val=x[0])
                myfile.write(defAlphaLine)

        def run_simulation(self):
            run_command = "{} +p{} {}.conf > out_{}.dat".format(
                str(self.gomc),
                str(1),
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
            df=pd.DataFrame({'steps':steps_np, np.round(x, 3) :energies_np})
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
            if (x[0] in self.obj_calls):
                return self.obj_calls[x[0]]
            Calibrator.copy_template(self,x)
            Calibrator.run_simulation(self)
            y_hat = Calibrator.extract_target(self, x[0])
            self.obj_calls[x[0]] = Calibrator.loss(self,y_hat)
            print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}'.format(self.iteration, x[0], self.target_y, y_hat, Calibrator.loss(self,y_hat)))
            with open("calibration.log", "a") as myfile:
                line = "{0: 3.6f}   {1: 3.6f}\n".format(x[0], Calibrator.loss(self,y_hat))
                myfile.write(line)

            self.iteration = self.iteration + 1
            return Calibrator.loss(self,y_hat)
        
        def calibrate(self):
            f = lambda x: Calibrator.objective(self,x)
            c = lambda x: Calibrator.constr(self,x)

            x0 = np.array([self.initial_x], dtype=np.double)
            res = \
                fmin_cobyla(f, 
                        x0, 
                        cons=[c],
                        rhobeg=0.01,
                        rhoend=0.005, 
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
