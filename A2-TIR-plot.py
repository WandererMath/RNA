import sys, os
from pickle import dump, load

import matplotlib.pyplot as plt

from common import *

EXTEND=30

PATH_INPUT=sys.argv[1]
PATH_OUTPUT=sys.argv[2]
FILE_FNA=sys.argv[3]
FILE_GTF=sys.argv[4]


letters=['A', 'T', 'C', 'G']

os.chdir(PATH_OUTPUT)

with open("TIR_Data.pickle","rb") as f:
    stat1=load(f)
    stat2=load(f)
    statb=load(f)


stat1r=[a/b for a, b in zip(stat1, statb)]
stat2r=[a/b for a, b in zip(stat2, statb)]
stat1r2r=[a/b for a, b in zip(stat1r, stat2r)]
X=list(range(-EXTEND, EXTEND))
for i in range(len(letters)):
    offset=i*2*EXTEND
    Y=stat1r[offset:2*EXTEND+offset]
    plt.plot(X, Y, label='C Up, B Down')

    Y=stat2r[offset:2*EXTEND+offset]
    plt.plot(X, Y, label='C Down, B Up')
    plt.legend()
    plt.ylabel("Frequency of {} / Frequency of {} of All Genes".format(letters[i], letters[i]))
    plt.xlabel("Position of the Sequence ( 0 is the gene start position )")
    plt.title("Nucleotide Relative Frequency - {}".format(letters[i]))
    plt.savefig(letters[i]+'.pdf')
    plt.cla()

    Y=stat1r2r[offset:2*EXTEND+offset]
    plt.plot(X, Y)
    plt.title("{} Relative Frequency Ratio (Stimulated Set/Reduced Set)".format(letters[i]))
    plt.xlabel("Position of the Sequence ( 0 is the gene start position )")
    plt.savefig(letters[i]+'_Ratio.pdf')
    plt.cla()
