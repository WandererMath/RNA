import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from itertools import combinations 
import os
import sys

os.chdir(sys.argv[1])
PATH="norm.csv"
TITLE="Normalized Gene Counts Scatter Plot"
PNG_DIR="PNG"

os.makedirs(PNG_DIR, exist_ok=True)

samples=[]
with open(PATH, "r") as f:
    line1=iter(next(f).split(","))
    next(line1)
    for s in line1:
        samples.append(s)

print(samples)
samples_index=[i for i in range(1,len(samples)+1)]
print(samples_index)
tasks=combinations(samples_index,2)

print(tasks)
def main(task):
    i1, i2=task
    path_png=os.path.join(PNG_DIR,samples[i1-1]+"-"+samples[i2-1]+".png")


    X_ori=[]
    Y_ori=[]
    with open(PATH, "r") as f:
        next(f)
        for line in f:
            try:
                data=line.split(",")
                X_ori.append(float(data[i1]))
                Y_ori.append(float(data[i2]))
                if float(data[i1])==0 or float(data[i2])==0:
                    continue
                    
                #X.append(log(float(data[i1]),10))
                #Y.append(log(float(data[i2]),10))
            except:
                pass

    plt.scatter(X_ori, Y_ori)
    correlation_coefficient, p_value = pearsonr(X_ori, Y_ori)
    plt.title(TITLE+"\nR="+str(correlation_coefficient))
    plt.xlabel(samples[i1-1])
    plt.ylabel(samples[i2-1])
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(path_png)
    print(correlation_coefficient)

    plt.clf()

for task in tasks:
    main(task)