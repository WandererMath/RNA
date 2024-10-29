import os
import sqlite3
DB="../ref/GENE_ID_NAME.db"

def get_txt_files(folder_path):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".txt"):
                txt_files.append(os.path.join(root, file))
    return txt_files

OUT_DIR="./output"
os.makedirs(OUT_DIR, exist_ok=True)
files=get_txt_files(".")
files_out=[os.path.join(OUT_DIR, file.split("/")[-1]) for file in files]
print(files_out)

with sqlite3.connect(DB) as CONN:
    CURSOR=CONN.cursor()
    for file, out_file in zip(files, files_out):
        print(out_file)
        with open(file, "r") as fr, open(out_file,"w") as fw: 
            for line in fr:
                line=line.split("\n")[0]
                CURSOR.execute("SELECT NAME FROM id_name where gene_id = '{}'".format(line))
                name=CURSOR.fetchall()[0][0]
                #name already contains \n
                fw.write(line+"\t"+ name)
                print(line+"\t"+ name)

    

            
    CURSOR.close()