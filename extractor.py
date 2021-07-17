import os
import pandas as pd
from pandas import read_fwf
import numpy as np
from numpy import empty
from database_api import Writer
from tqdm import tqdm
from my_thread import ThreadWithReturn

#KDD cup pytorch geometric compatible dataset
# data_df = pd.read_csv('/data3/kaleb.dickerson2001/Datasets/KDD-pyg-dataset/pcqm4m_kddcup2021/raw/data.csv.gz')
# extractor_log = open("extractor_log.txt", "w")
# removed_log = open("removed_log.txt", "w")
db_writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/smile_coord.db")


def isolate_smiles(file_path) -> str:
    """
    input: .info file
    clears out extraneous information
    return: smiles string
    """
    # simply returns smiles
    smile = pd.read_fwf(file_path).columns[3]
    return smile

def get_coordinates(file_path) -> list:
    """
    input: .xyz file
    extracts information as list of tuples
    return: list of tuples
    """
    file = read_fwf(file_path)
    arr_x = file.iloc[:,[1]].astype('float32').to_numpy().flatten()
    arr_y = file.iloc[:,[2]].astype('float32').to_numpy().flatten()
    arr_z = file.iloc[:,[3]].astype('float32').to_numpy().flatten()
    length = len(arr_x)

    to_return = empty(length, dtype=tuple)
    for i in range(length):
        to_return[i] = (arr_x[i],arr_y[i],arr_z[i])
    return to_return



# def smile_in_KDD_dataset(smile, index) -> bool:
#     """
#     specifies if molecule is present in KDD dataset
#     """
#     i = data_df[data_df.smiles == smile].first_valid_index()
#     if i == None:
#         removed_log.write(f"deleted unpresent {index} | {smile}\n")
#         return False
#     else:
#         return True

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
    
    os.system(f"cd {path_to_tars} && tar xf {file}")
    path_to_extracted = "".join([path_to_tars, f"/Compound_{start_pad}_{end_pad}"])
    # delete all files but smiles and coordinates
    # os.system(f"cd {path_to_extracted} && find . -name *.PM6.* -type f -delete")

    for i in range(start, end+1):
        i_pad = str(i).zfill(9)
        # path to specific molecule directory
        mol_dir_path = "".join([path_to_extracted, f"/{i_pad}"])
        # base path for files inside the molecule directory
        file_base_path = "".join([mol_dir_path, f"/{i_pad}"])
        # checks if molecule number directory exists
        exists = os.path.isdir(mol_dir_path)

        if exists:
            # returns smile string and removes other info from the file
            # smile = isolate_smiles(f"{file_base_path}.20160829.info")
            # coords = get_coordinates(f"{file_base_path}.initial.xyz")
            t1 = ThreadWithReturn(target=isolate_smiles, args=(f"{file_base_path}.20160829.info",))
            t2 = ThreadWithReturn(target=get_coordinates, args=(f"{file_base_path}.initial.xyz",))
            t1.start(),t2.start()
            smile = t1.join()
            coords = t2.join()

            db_writer.add_entry(smile, coords)
            # if the molecule isnt presesnt in KDD dataset delete it
            # if not smile_in_KDD_dataset(smile, i_pad):
            #     os.system(f"cd {path_to_extracted} && rm -r {i_pad}")
        # else:
            # if molecule file doesnt exist log it
            # extractor_log.write(f"No Molecule {i_pad}\n")
    os.system(f"cd {path_to_tars} && rm -r Compound_{start_pad}_{end_pad}")
# extractor_log.close()
# removed_log.close()
