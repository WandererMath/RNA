import os
def get_files(folder_path, end):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(end):
                txt_files.append(os.path.join(root, file))
    return txt_files