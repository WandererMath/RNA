import sys, os
from pickle import dump, load
from scipy import stats
import matplotlib.pyplot as plt


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
    N1=load(f)
    N2=load(f)
    N3=load(f)


stat1r=[a/b for a, b in zip(stat1, statb)]
stat2r=[a/b for a, b in zip(stat2, statb)]
stat1r2r=[a/b for a, b in zip(stat1r, stat2r)]

X=list(range(-EXTEND, EXTEND))

def to_t(Y, N):
    return [int(N1*a)*[1]+int(N1*(1-a))*[0] for a in Y]

for j in range(len(letters)):
    offset=j*2*EXTEND
    Y1=stat1[offset:2*EXTEND+offset]

    Y2=stat2[offset:2*EXTEND+offset]

    Ref=statb[offset:2*EXTEND+offset]

    Y1_t=to_t(Y1, N1)
    Y2_t=to_t(Y2, N2)
    Ref_t=to_t(Ref, N3)

    T1=[]
    T2=[]
    #print(Y1_t, Ref_t)
    for i in range(len(Y1_t)):

        # Paired samples t-test
        _, p1 = stats.ttest_ind(Ref_t[i], Y1_t[i], equal_var=False)
        _, p2 = stats.ttest_ind(Ref_t[i], Y2_t[i], equal_var=False)
        T1.append(p1)
        T2.append(p2)
    plt.plot(X, T1, label='Stimulated')
    plt.plot(X, T2, label="Reduced")
    plt.ylabel("p-value")
    plt.xlabel("Position of the Sequence Relative to the Gene Start Position")
    plt.legend()
    plt.title("p-values - Nucleotide {}".format( letters[j]))
    plt.savefig(letters[j]+"-t-test.pdf")
    plt.cla()

    #print(letters[i],"p-value", p_value)