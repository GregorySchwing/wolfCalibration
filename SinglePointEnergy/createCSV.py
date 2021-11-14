import os
import numpy as np
import shutil
import subprocess
import csv
import pandas as pd
import pickle
# range for floats with np.arange()
count = 0
methods = [ "Gross" , "Vlugt", "Cassandra" ]
potentials = [ "DSP" , "DSF"]
totalElecIndex = 0
RCutDF = list()
for Rcut in np.arange (10, 16, 1):
    filename = "{}_error.csv".format(Rcut)
    with open(filename, 'w', encoding='UTF8') as f:
        directory_name = "rCutCoul_{}/Ewald".format(Rcut)
        logFileName = "Ewald_log.txt"
        logFilePath = directory_name+"/"+logFileName
        proc = subprocess.Popen(["grep", "-P", "\s+ETITLE:", logFilePath], stdout=subprocess.PIPE)
        (out, err) = proc.communicate()
        g = out.split()  # tokenize the string   
        g.pop(0)
        print(g)
        totalElecIndex = g.index(b'TOTAL_ELECT')
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
        #header = []
        for method in methods:
            for potential in potentials:
                colName = "{}_{}".format(method, potential)
                header.append(colName)
        writer.writerow(header)
        dataFrameList = list()
        print(Rcut, end=', ')
        for alpha in np.arange (0.005, 0.305, 0.005):
            errValues = []
            indices = []            
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
            dataFrameList.append(errValues)
            writer.writerow(errValues)

        df = pd.DataFrame(dataFrameList, columns=header)
        df.set_index('Alpha')
        RCutDF.append(df)
        print(df)

RCut = np.arange (10, 16, 1)
zipbObj = zip(RCut, RCutDF)
dfAll = pd.concat(dict(zipbObj),axis=1, names=["RCut", "Method"])
print(dfAll)
dfAll.to_pickle('test_pickle.p')
