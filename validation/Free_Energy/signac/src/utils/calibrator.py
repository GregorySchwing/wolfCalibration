import os
import subprocess
import shutil
import re
import numpy as np
from pymbar import timeseries
from scipy.optimize import fmin_bfgs
class Calibrator:
    def __init__(self, gomc, wolf_model, wolf_potential, target_y, initial_x, template_directory, template_control_file_name_str, \
                 conffile, forcefield, num_iters=20):
        self.gomc = gomc
        self.wolf_model = wolf_model
        self.wolf_potential = wolf_potential
        self.target_y = target_y
        self.initial_x = initial_x
        self.template_directory = template_directory
        self.template_control_file_name_str = template_control_file_name_str
        self.conffile = conffile
        self.forcefield = forcefield
        self.iteration = 0
        self.num_iters = num_iters

    
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
            defAlphaLine = "WolfAlpha\t{box}\t{val}\n".format(box=0, val=x)
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

    def extract_target(self):
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
        nskip = 20
        steps_np = np.array(steps)
        energies_np = np.array(energies)
        t0, g, Neff_max = timeseries.detectEquilibration(energies_np, nskip=nskip) # compute indices of uncorrelated timeseries
        A_t_equil = energies_np[t0:]
        A_t_equil_steps = steps_np[t0:]

        return A_t_equil.mean()
        


    # objective function
    def objective(self, x):
        Calibrator.copy_template(self,x)
        Calibrator.run_simulation(self)
        y_hat = Calibrator.extract_target(self)
        print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}'.format(self.iteration, x, self.target_y, y_hat, np.abs(self.target_y-y_hat)))
        self.iteration = self.iteration + 1
        return np.abs(self.target_y-y_hat)
    
    def calibrate(self):
        f = lambda x: Calibrator.objective(self,x)
        x0 = np.array([self.initial_x], dtype=np.double)
        [xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflg] = \
            fmin_bfgs(f, 
                    x0, 
                    #callback=callbackF, 
                    maxiter=3, 
                    full_output=True, 
                    retall=False)
