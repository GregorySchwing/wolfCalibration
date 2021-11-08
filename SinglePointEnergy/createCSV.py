import os
import numpy as np
import shutil
import subprocess
import csv

# range for floats with np.arange()
count = 0
methods = [ "Gross" , "Cassandra" , "Vlugt" ]
potentials = [ "DSP" , "DSF"]
totalElecIndex = 0
for Rcut in np.arange (10, 16, 1):
    filename = "{}_error.csv".format(Rcut)
    with open(filename, 'w', encoding='UTF8') as f:
        directory_name = "rCutCoul_{}/Ewald".format(Rcut)
        logFileName = "Ewald_log.txt"
        logFilePath = directory_name+"/"+logFileName
        proc = subprocess.Popen(["grep", "-P", "ETITLE:", logFilePath], stdout=subprocess.PIPE)
        (out, err) = proc.communicate()
        g = out.split()  # tokenize the string   
        g.pop(0)
        totalElecIndex = g.index("TOTAL_ELECT")
        proc = subprocess.Popen(["grep", "-P", "ENER_0:", logFilePath], stdout=subprocess.PIPE)
        (out, err) = proc.communicate()
        g = out.split()  # tokenize the string           
        g.pop(0)
        # write the data
        energyValuesEwald = np.asarray(g, dtype=np.double)
        ewaldTotElec = energyValuesEwald[totalElecIndex]
        ewaldAbs = np.abs(ewaldTotElec)
        writer = csv.writer(f, delimiter=' ')
        # write the header
        header = ["Alpha"]
        for method in methods:
            for potential in potentials:
                colName = "{}_{}".format(method, potential)
                header.append(colName)
        #writer.writerow(header)
        print(Rcut, end=', ')
        for alpha in np.arange (0.005, 0.305, 0.005):
            errValues = []
            errValues.append(round(alpha,3))
            directory_name = "rCutCoul_{}/Alpha_{:.3f}".format(Rcut, alpha)
            for method in methods:
                for potential in potentials:
                    logFileName = "{}_{}_log.txt".format(method, potential)
                    logFilePath = directory_name+"/"+logFileName
                    proc = subprocess.Popen(["grep", "-P", "ENER_0:", logFilePath], stdout=subprocess.PIPE)
                    (out, err) = proc.communicate()
                    g = out.split()  # tokenize the string           
                    g.pop(0)
                    energyValues = np.asarray(g, dtype=np.double)
                    totElect = energyValues[totalElecIndex]
                    ERR = (np.abs(totElect) - ewaldAbs)/ewaldAbs
                    errValues.append(ERR)
            writer.writerow(errValues)
            
