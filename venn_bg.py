import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import sys
import os

# Arg1 input_dir, Arg2 output_dir
original_dir = os.getcwd()
os.chdir(sys.argv[1])
OUT_DIR=sys.argv[2]

PATH_AB="result-1.tsv"
PATH_AC="result-2.tsv"
#threshold=float(sys.argv[1])
threshold=1
p_threshold=0.05

ID=0
PADJ=3
LOG2FC=-1




def analysis2(path):
    up=[]
    down=[]
    c=[]
    with open(path, 'r') as f:
        next(f)
        for line in f:
            line=line.split()
            id=line[ID]
            padj=float(line[PADJ])
            l2fc=float(line[LOG2FC])
            if padj<=p_threshold and l2fc>=threshold:
                up.append(id)
            if padj<=p_threshold and l2fc<=-threshold:
                down.append(id) 
            if l2fc>=0:
                c.append(id)
    return up, down,c

B_p, B_m,c=analysis2(PATH_AB)



os.chdir(original_dir)
os.makedirs(OUT_DIR, exist_ok=True)
os.chdir(OUT_DIR)


def save(genes, path):
    with open (path, "w") as f:
        for gene in genes:
            f.write(gene+"\n")

genes=[c]
paths=['bg4Bm2.txt']
for gene, path in zip(genes, paths):
    save(gene, path)