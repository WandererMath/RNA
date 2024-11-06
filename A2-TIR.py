import sys, os
from pickle import dump

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import matplotlib.pyplot as plt

from common import *

EXTEND=30

PATH_INPUT=sys.argv[1]
PATH_OUTPUT=sys.argv[2]
FILE_FNA=sys.argv[3]
FILE_GTF=sys.argv[4]




# Load your GTF file and create a database (it will be created if it doesn't exist)
db = gffutils.create_db(FILE_GTF, dbfn='GTF.db', force=True, keep_order=True, merge_strategy='merge')

# Specify your gene ID

def get_start(gene_id):
    # Retrieve the gene entry by its ID
    try:
        gene = db[gene_id]
        start_position= gene.start
        return start_position
        #print(f"The start position of gene {gene_id} is {start_position}")
    except KeyError:
        print(f"Gene ID {gene_id} not found in the GTF file.")

def get_strand_direction(gene):
    return db[gene].strand


def get_seq(start, end, strand="+"):
    with open(FILE_FNA, 'r') as fna_file:
        for record in SeqIO.parse(fna_file, 'fasta'):
            # Get the sequence from start to end position
            sequence = record.seq[start - 1:end]  # Adjust to 0-based indexing
            #print(f"Extracted sequence from position {start} to {end}:")
            #print(sequence)
            if strand=="-":
                s=Seq(sequence)
                s=s.reverse_complement()
                sequence=str(s)
            return sequence
        

F1=os.path.join(PATH_INPUT, "BmCp.txt")
F2=os.path.join(PATH_INPUT, "BpCm.txt")
F3=os.path.join(PATH_INPUT, "all_true_genes.txt")
# 1A*2EXTEND T C G
stat1=[0 for i in range(4*2*EXTEND)]
stat2=[0 for i in range(4*2*EXTEND)]
statb=[0 for i in range(4*2*EXTEND)]
letters=['A', 'T', 'C', 'G']


def count(f, stat):
    N=0
    for line in f:
        N+=1
        # remove \n
        gene=line.split()[0]
        start=get_start(gene)
        d=get_strand_direction(gene)
        seq=get_seq(start-EXTEND, start+EXTEND-1, strand=d)
        if not len(seq)==2*EXTEND:
            raise Exception
        for i in range(len(letters)):
            for n in range(len(seq)):
                if letters[i]==seq[n]:
                    stat[i*2*EXTEND+n]+=1
    return N

with open(F1,'r') as f1,open(F2,'r') as f2, open(F3,'r') as fb:
    print("\n\n")
    N1=count(f1, stat1)
    N2=count(f2, stat2)
    N3=count(fb, statb)

stat1=[a/N1 for a in stat1]
stat2=[a/N2 for a in stat2]
statb=[a/N3 for a in statb]

stat1r=[a/b for a, b in zip(stat1, statb)]
stat2r=[a/b for a, b in zip(stat2, statb)]

stat1r2r=[a/b for a, b in zip(stat1r, stat2r)]

os.chdir(PATH_OUTPUT)

with open("TIR_Data.pickle","wb") as f:
    dump(stat1, f)
    dump(stat2, f)
    dump(statb, f)
    dump(N1, f)
    dump(N2, f)
    dump(N3, f)


X=list(range(-EXTEND, EXTEND))
for i in range(len(letters)):
    offset=i*2*EXTEND
    Y=stat1r[offset:2*EXTEND+offset]
    plt.plot(X, Y, label='C Up, B Down')

    Y=stat2r[offset:2*EXTEND+offset]
    plt.plot(X, Y, label='C Down, B Up')
    plt.legend()
    plt.ylabel("#{} / #{} of All Genes".format(letters[i], letters[i]))
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

