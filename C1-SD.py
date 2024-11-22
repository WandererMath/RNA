from subprocess import run, PIPE
import os, sys
from pickle import dump


from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import matplotlib.pyplot as plt
from numpy import array

# Anti-SD seq on rRNA
ASD='ACCUCCUUA'


EXTEND=50

PATH_INPUT=sys.argv[1]
PATH_OUTPUT=sys.argv[2]
FILE_FNA=sys.argv[3]
FILE_GTF=sys.argv[4]

# Load your GTF file and create a database (it will be created if it doesn't exist)
db = gffutils.create_db(FILE_GTF, dbfn='GTF.db', force=True, keep_order=True, merge_strategy='merge')



F1=os.path.join(PATH_INPUT, "BmCp.txt")
F2=os.path.join(PATH_INPUT, "BpCm.txt")
F3=os.path.join(PATH_INPUT, "all_true_genes.txt")

def cmd(s1, s2):
    return "echo \"{}\n{}\" | RNAduplex".format(s1, s2)


# Returns positions of seq2
def duplex(seq1, seq2):
    result = run(cmd(seq1, seq2),shell=True, stdout=PIPE, stderr=PIPE, text=True)
    out=result.stdout
    E=float(out.split('(')[-1].split(')')[0])
    A2=out.split(':')[-1].split('(')[0]
    s2, e2=A2.split(',')
    s2=int(s2)
    e2=int(e2)
    return float(E), s2, e2







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
        
def spacings(file):
    # E, spacing pairs
    stat=[]
    with open(file,'r') as f:
        for line in f:
            # remove \n
            gene=line.split()[0]
            #print(gene)
            start=get_start(gene)
            d=get_strand_direction(gene)
            seq=get_seq(start-EXTEND, start-1+EXTEND, strand=d)

            E, s, _=duplex(ASD, seq)
            stat.append([E, s-EXTEND-1])
    return stat
            


stat1=spacings(F1)
stat2=spacings(F2)
stat3=spacings(F3)

os.chdir(PATH_OUTPUT)
with open("data.pickle",'wb') as f:
    dump(stat1, f)
    dump(stat2, f)
    dump(stat3, f)

E1, spacing1=list(array(stat1).transpose())
E2, spacing2=list(array(stat2).transpose())
E3, spacing3=list(array(stat3).transpose())

plt.hist(spacing1, bins=EXTEND)
plt.savefig("Stimulated.pdf")
plt.cla()

plt.hist(spacing2, bins=EXTEND)
plt.savefig("Reduced.pdf")
plt.cla()

plt.hist(spacing3, bins=EXTEND)
plt.savefig("All.pdf")
plt.cla()

plt.hist(E1)
plt.savefig("Stimulated_G.pdf")
plt.cla()

plt.hist(E2)
plt.savefig("Reduced_G.pdf")
plt.cla()

plt.hist(E3)
plt.savefig("All_G.pdf")
plt.cla()