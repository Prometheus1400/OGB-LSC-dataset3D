import os
import numpy as np
from numpy import empty, loadtxt
from database_api import Writer
from tqdm import tqdm, trange
from utils import get_coordinates


db_writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/smile_coord.db")
error_log = open("error_log.txt", "w")
write_log = open("write_log", "w")


# each tar file contains information for ~25,000 molecules
path_to_tars = "/data3/kaleb.dickerson2001/Datasets/PubChem3D"
file_list = []
for file in os.listdir(path_to_tars):
    if file.endswith(".tar.gz"):
        file_list.append(file)
file_list.sort()

# array to store smile,coordinate pairs before writing to database
SIZE_ARR = 10000
arr = np.empty(SIZE_ARR, dtype=tuple)
for file in tqdm(file_list):
    start, end = int(file[9:18]), int(file[19:28])
    start_pad = str(start).zfill(9)
    end_pad = str(end).zfill(9)
    
    # extract the files
    os.system(f"cd {path_to_tars} && tar xf {file}")
    path_to_extracted = "".join([path_to_tars, f"/Compound_{start_pad}_{end_pad}"])
    # set array size to 0
    size = 0

    for i in range(start, end+1):
        i_pad = str(i).zfill(9)
        # path to specific molecule directory
        mol_dir_path = "".join([path_to_extracted, f"/{i_pad}"])

        # checks if molecule number directory exists
        exists = os.path.isdir(mol_dir_path)
        if exists:
            # base path for files inside the molecule directory
            file_base_path = "".join([mol_dir_path, f"/{i_pad}"])
            CID = i
            coords = get_coordinates(f"{file_base_path}.b3lyp_6-31g(d).inp", error_log=error_log)

            # write pair to array
            arr[size] = (CID, coords)
            size += 1
            # writes to database in batches of 10000
            if size == SIZE_ARR:
                db_writer.add_entries(arr[:size])
                size = 0

    db_writer.add_entries(arr[:size])
    # remove extracted file, would quickly consume too much storage otherwise
    os.system(f"cd {path_to_tars} && rm -r Compound_{start_pad}_{end_pad}")
    write_log.write(f"Finished processing: {file}\n")
error_log.close()