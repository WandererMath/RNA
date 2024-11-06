import sys

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils

from common import *

EXTEND=25

FILE_GENES=sys.argv[1]
FILE_OUTPUT=sys.argv[2]
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
        




with open(FILE_GENES,'r') as f, open(FILE_OUTPUT,'w') as fw:
    for line in f:
        # remove \n
        gene=line.split()[0]
        start=get_start(gene)
        d=get_strand_direction(gene)
        seq=get_seq(start-EXTEND, start+EXTEND, strand=d)
        #seq=get_seq(start-EXTEND, start+EXTEND)

        fw.write(f">{gene}\n{seq}\n")
        #fw.write(f">{gene} strand={d}\n{seq}\n")
