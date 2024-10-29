import os
import sys
from copy import deepcopy


# Enter 1 or 0
AvsC=int(sys.argv[1])
#folder_path = "../gene_counts_feature"
folder_path = sys.argv[2]
output_path=sys.argv[3]

os.makedirs(output_path, exist_ok=True)

def get_txt_files(folder_path):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".txt"):
                txt_files.append(os.path.join(root, file))
    return txt_files


file_names = get_txt_files(folder_path)


for file in deepcopy(file_names):
    if AvsC:
        if not ("A0" in file or "C0" in file):
            file_names.remove(file)
    else:
        if not ("A0" in file or "B0" in file):
            file_names.remove(file)

print(file_names)
print(len(file_names))

if AvsC:
    fw=open(os.path.join(output_path,"genecount_matrix-2.csv"),"w")
    #fw.write("GeneID,A01,A02,A03,C01,C02,C03\n")
else:
    fw=open(os.path.join(output_path,"genecount_matrix-1.csv"),"w")
    #fw.write("GeneID,A01,A02,A03,B01,B02,B03\n")

# Write table header
fw.write("GeneID")
for file in file_names:
    fw.write(","+file)
fw.write("\n")

F=[]
for file in file_names:
    F.append(open(file, "r"))

# Skip first two lines
for f in F:
    next(f)
    next(f)


while True:
    try:
        counts=[]
        flag=True
        for f in F:
            line=next(f)
            # For the first file, append gene_id
            if flag:
                flag=False
                counts.append(line.split()[0])
            # Append Gene count
            counts.append(line.split()[-1])
        fw.write(",".join(counts)+"\n")

    except StopIteration:
        break

fw.close()
for f in F:
    f.close()


#generate metadata
print("Generating metadata")
if AvsC:
    PATH_METADATA=os.path.join(output_path,"metadata-2.csv")
else:
    PATH_METADATA=os.path.join(output_path,"metadata-1.csv")

with open(PATH_METADATA, "w") as f:
    f.write("id, dex\n")
    for file in file_names:
        f.write(file)
        if "A0" in file:
            f.write(",control\n")
        else:
            f.write(",treated\n")