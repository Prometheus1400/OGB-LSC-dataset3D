import os
import numpy as np
from numpy import empty, loadtxt
from database_api import Writer
from tqdm import tqdm, trange
from my_thread import ThreadWithReturn


db_writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/smile_coord.db")
error_log = open("error_log.txt", "w")

def get_smiles(file_path) -> str:
    """
    input: .info file
    return: smiles string
    """
    # simply returns smiles
    try:
        smile = loadtxt(file_path, dtype=str)[3]
        return smile
    except:
        error_log.write(f"error in get_smiles at file: {file_path}\n")
        return "NONE"

def get_coordinates(file_path) -> np.array:
    """
    input: .xyz file
    extracts information as array of tuples
    return: array of tuples
    """
    try:
        data = loadtxt(file_path,delimiter="\n",dtype=str)[1:]
        length = len(data)
        to_return = empty(length, dtype=tuple)
        for i in range(length):
            line = data[i].split()
            to_return[i] = (line[1], line[2], line[3])
        return to_return
    except:
        error_log.write(f"error in get_coordinates at file: {file_path}\n")
        return -1


# each tar file contains information for ~25,000 molecules
path_to_tars = "/data3/kaleb.dickerson2001/Datasets/PubChem3D"
file_list = []
for file in os.listdir(path_to_tars):
    if file.endswith(".tar.xz"):
        file_list.append(file)
file_list.sort()


for file in tqdm(file_list):
    start, end = int(file[9:18]), int(file[19:28])
    start_pad = str(start).zfill(9)
    end_pad = str(end).zfill(9)
    
    # extract the files
    os.system(f"cd {path_to_tars} && tar xf {file}")
    path_to_extracted = "".join([path_to_tars, f"/Compound_{start_pad}_{end_pad}"])
    # array to store smile,coordinate pairs before writing to database
    arr = np.empty(5000, dtype=tuple)
    size = 0

    for i in range(start, end+1):
        i_pad = str(i).zfill(9)
        # path to specific molecule directory
        mol_dir_path = "".join([path_to_extracted, f"/{i_pad}"])
        # base path for files inside the molecule directory
        file_base_path = "".join([mol_dir_path, f"/{i_pad}"])
        # checks if molecule number directory exists
        exists = os.path.isdir(mol_dir_path)

        if exists:
            # two different threads to read the information in from files
            t1 = ThreadWithReturn(target=get_smiles, args=(f"{file_base_path}.20160829.info",))
            t2 = ThreadWithReturn(target=get_coordinates, args=(f"{file_base_path}.initial.xyz",))
            t1.start(),t2.start()
            smile = t1.join()
            coords = t2.join()

            # write pair to array
            arr[size] = (smile, str(coords))
            size += 1
            # writes to database in batches of 5000
            if size == 5000:
                db_writer.add_entries(arr[:size])
                size = 0

    db_writer.add_entries(arr[:size])
    # remove extracted file, would quickly consume too much storage otherwise
    os.system(f"cd {path_to_tars} && rm -r Compound_{start_pad}_{end_pad}")
error_log.close()