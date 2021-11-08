import os
import numpy as np
import shutil

# range for floats with np.arange()
count = 0
methods = [ "Gross" , "Vlugt", "Cassandra" ]
potentials = [ "DSP" , "DSF"]
for Rcut in np.arange (10, 16, 1):
    directory_name = "rCutCoul_{}/Ewald".format(Rcut)
    confFileName = "water_RCutCoul_{}_Ewald.conf".format(Rcut)
    try:
        os.makedirs(directory_name)
    except OSError as error:
            print(error)
    try:
        shutil.copyfile("water_RCutCoul_AAA_Ewald.conf", directory_name+"/"+confFileName)
    except OSError as error:
        print(error)
    sedString1 = "sed -i \'s/AAA/{}/\' {}".format(Rcut, directory_name+"/"+confFileName)
    os.system(sedString1)
    os.chdir(directory_name)
    command = "~/GOMC/bin/GOMC_CPU_NPT {} > Ewald_log.txt".format(confFileName)
    os.system(command)
    os.chdir("../..")


    print(Rcut, end=', ')
    for alpha in np.arange (0.005, 0.305, 0.005):
        directory_name = "rCutCoul_{}/Alpha_{:.3f}".format(Rcut, alpha)
        try:
            os.makedirs(directory_name)
        except OSError as error:
                print(error)
        for method in methods:
            for potential in potentials:
                confFileName = "RCutCoul_{}_Alpha_{:.3f}_Method_{}_Potential_{}.conf".format(Rcut, alpha, method,potential)

                try:
                    shutil.copyfile("water_RCutCoul_AAA_Alpha_ZZZ_Method_XXX_Potential_YYY.conf", directory_name+"/"+confFileName)
                except OSError as error:
                    print(error)

                sedString1 = "sed -i \'s/AAA/{}/\' {}".format(Rcut, directory_name+"/"+confFileName)
                sedString2 = "sed -i \'s/ZZZ/{:.3f}/\' {}".format(alpha, directory_name+"/"+confFileName)
                sedString3 = "sed -i \'s/XXX/{}/\' {}".format(method, directory_name+"/"+confFileName)
                sedString4 = "sed -i \'s/YYY/{}/\' {}".format(potential, directory_name+"/"+confFileName)
                os.system(sedString1)
                os.system(sedString2)
                os.system(sedString3)
                os.system(sedString4)
                os.chdir(directory_name)
                command = "~/GOMC/bin/GOMC_CPU_NPT {} > {}_{}_log.txt".format(confFileName, method, potential)
                os.system(command)
                os.chdir("../..")
