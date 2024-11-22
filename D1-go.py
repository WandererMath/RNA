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
FILE_GTF=sys.argv[3]



F1=os.path.join(PATH_INPUT, "BmCp.txt")
F2=os.path.join(PATH_INPUT, "BpCm.txt")
F3=os.path.join(PATH_INPUT, "all_true_genes.txt")


GENES={}
with open(FILE_GTF, 'r') as f:
    for line in f:
        flag1=False
        flag2=False
        flag3=False
        line=line.split()
        for s in line:
            if s=="gene_id":
                flag1=True
                continue
            if flag1:
                id=s.replace('\"', '').replace(";", "")
                if not id in GENES:
                    GENES[id]=[]
                flag1=False
                flag2=True
                continue
            
            # Applend Ontology
            if flag2:
                if s=="Ontology_term":
                    flag3=True
                    continue
                if flag3:
                    flag3=False
                    ontology=s.replace("\"", "").replace(";", "")
                    GENES[id].append(ontology)
                    continue

def stat(file):
    result={}
    with open(file, 'r') as f:
        for line in f:
            gene=line.replace('\n', '')
            for ontology in GENES[gene]:
                if ontology not in result:
                    result[ontology]=0
                else:
                    result[ontology]+=1
    return result





stimulated=stat(F1)
reduced=stat(F2)

stimulated_sorted_list=sorted(list(stimulated),reverse=True, key=lambda x: stimulated[x])
reduced_sorted_list=sorted(list(reduced),reverse=True, key=lambda x: reduced[x])

os.chdir(PATH_OUTPUT)
with open("GENES.pickle", 'wb') as f:
    dump(GENES, f)

with open("stimulated.txt", 'w') as f:
    for ontology in stimulated_sorted_list:
        f.write(ontology+','+str(stimulated[ontology])+'\n')

with open("reduced.txt", 'w') as f:
    for ontology in reduced_sorted_list:
        f.write(ontology+','+str(reduced[ontology])+'\n')



      
#print(GENES)