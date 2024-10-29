import sys, os

AvsC=int(sys.argv[1])
PATH_M=sys.argv[2]
PATH_R=sys.argv[3]
PATH_OUTPUT=sys.argv[4]

if AvsC:
    FILE="genecount_matrix-2.csv"
    FILE_META="meta-2.csv"
    FILE_TABLE='matrix-2.csv'
else:
    FILE="genecount_matrix-1.csv"
    FILE_META="meta-1.csv"
    FILE_TABLE='matrix-1.csv'

FILE_M=os.path.join(PATH_M, FILE)
FILE_R=os.path.join(PATH_R, FILE)

OUT_META=os.path.join(PATH_OUTPUT, FILE_META)
OUT_TABLE=os.path.join(PATH_OUTPUT, FILE_TABLE)

with open(FILE_M, "r") as fM,\
open(FILE_R,"r") as fR,\
open(OUT_META, "w") as fMETA,\
open(OUT_TABLE, "w") as fTABLE:
    #[(name_brief, type, condition),...]
    head_info_M=[]
    head_info_R=[]

    # Generate Meta
    fMETA.write('Samples,Data_Type,Conditions\n')
    
    head_M=next(fM).split(",")
    data_type="RNA-seq"
    for name in head_M[1:]:
        if "A0" in name:
            s_index=name.find("A0")
            head_info_M.append(["M"+name[s_index: s_index+3], data_type,"Control"])
        else:
            if "B0" in name:
                s_index=name.find("B0")
            else:
                s_index=name.find("C0")
            head_info_M.append(["M"+name[s_index: s_index+3], data_type,"DrugTreated"])
    print(head_info_M)
    for sample in head_info_M:
        fMETA.write(",".join(sample)+'\n')


    head_M=next(fR).split(",")
    data_type="Ribo-seq"
    for name in head_M[1:]:
        if "A0" in name:
            s_index=name.find("A0")
            head_info_R.append(["R"+name[s_index: s_index+3], data_type,"Control"])
        else:
            if "B0" in name:
                s_index=name.find("B0")
            else:
                s_index=name.find("C0")
            head_info_R.append(["R"+name[s_index: s_index+3], data_type,"DrugTreated"])
    print(head_info_R)
    
    
    for sample in head_info_R:
        fMETA.write(",".join(sample)+'\n')
    
    # Generate Count Matrix
    fTABLE.write("Entry\t"+"\t".join(head_info_M[i][0] for i in range(len(head_info_M)))+'\t'+'\t'.join(head_info_R[i][0] for i in range(len(head_info_R)))+"\n")
    for lineM, lineR in zip(fM, fR):
        #remove \n
        lineM=lineM[:-1]
        lineR=lineR[:-1]
        Ms=lineM.split(",")
        fTABLE.write(Ms[0])
        for c in Ms[1:]:
            fTABLE.write("\t"+c)
        
        Rs=lineR.split(',')
        for c in Rs[1:]:
            fTABLE.write("\t"+c)
        fTABLE.write("\n")
