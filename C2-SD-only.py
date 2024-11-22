from subprocess import run, PIPE
import os, sys
from pickle import dump, load




import matplotlib.pyplot as plt
from numpy import array, vectorize
from scipy import stats
from statistics import mean, mode

# Anti-SD seq on rRNA
ASD='ACCUCCUUA'


EXTEND=50

PATH_INPUT=sys.argv[1]
PATH_OUTPUT=sys.argv[2]
FILE_FNA=sys.argv[3]
FILE_GTF=sys.argv[4]



F1=os.path.join(PATH_INPUT, "BmCp.txt")
F2=os.path.join(PATH_INPUT, "BpCm.txt")
F3=os.path.join(PATH_INPUT, "all_true_genes.txt")

            



os.chdir(PATH_OUTPUT)
with open("data.pickle",'rb') as f:
    stat1=load(f)
    stat2=load(f)
    stat3=load(f)


E1, spacing1=list(array(stat1).transpose())
E2, spacing2=list(array(stat2).transpose())
E3, spacing3=list(array(stat3).transpose())

print(mean(E1), mean(E2), mean (E3))

print(stats.ttest_ind(E1, E2))
print(stats.ttest_ind(spacing1, spacing2))



################


space=[spacing1, spacing2, spacing3]
modes=[mode(s) for s in space]
N=(len(s) for s in space)

peak_c=[]
for s, M in zip(space,modes):
    x=0
    for t in s:
        if abs(t-M)<=3:
            x+=1
    peak_c.append(x)

peak_freq=[c/n for c, n in zip(peak_c, N)]
print(peak_freq)
