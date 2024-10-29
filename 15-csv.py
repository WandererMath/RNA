import os
import sys



folder_path = sys.argv[1]
output_path=sys.argv[2]


def get_txt_files(folder_path):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".txt"):
                txt_files.append(os.path.join(root, file))
    return txt_files


file_names = get_txt_files(folder_path)




print(file_names)
print(len(file_names))


fw=open(os.path.join(output_path,"genecounts.csv"),"w")
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

PATH_METADATA=os.path.join(output_path,"metadata.csv")

with open(PATH_METADATA, "w") as f:
    f.write("id, dex\n")
    for file in file_names:
        f.write(file)
        if "A0" in file:
            f.write(",control\n")
        else:
            f.write(",treated\n")