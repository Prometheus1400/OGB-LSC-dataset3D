from typing import List
from rdkit import Chem
import pubchempy
from numpy import loadtxt
from threading import Thread
import numpy as np
import os
import pubchempy as pcp
from tqdm import tqdm
from cclib.io import ccread

class ThreadWithReturn(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None
    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                                **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return

def get_coordinates(file_path, error_log=None) -> str:
    """
    input: .inp file containing optimized coords
    return: string of coords seperated by space
    """
    try:
        data = loadtxt(file_path,delimiter="\n",dtype=str)[9:][:-1]
        length = len(data)
        # to_return = empty(length, dtype=tuple)
        to_return = ""
        first = False
        for i in range(length):
            line = data[i].split()
            s = " ".join(line[2:])
            to_return += " " + s
        return to_return[1:] #remove leading space
    except:
        if error_log != None:
            error_log.write(f"error in get_coordinates at file: {file_path}\n")
        return -1

def split_into(iterable, num_splits):
    chunk = np.array_split(np.array(iterable),num_splits)
    for i in range(len(chunk)):
        chunk[i] = list(chunk[i])
    
    return chunk

def smile2CID(smile):
    compounds = pubchempy.get_compounds(smile, namespace='smiles')
    cid = compounds[0].cid
    return cid

def getHOMOLUMO(log_file):
    """
    input: path to log file for B3LYP optimized geometry
    output: HOMO/LUMO gap as float
    """
    data = ccread(log_file)
    homo = data.homos[0]
    energies = data.moenergies[0]
    homolumogap = energies[homo+1] - energies[homo]
    return homolumogap




# def combine_databases(target_db:str, path:str) -> None:
#     """
#     combines all databases from target directory with
#     the target database
#     """
#     files = []
#     for file in os.listdir(path):
#         if file.endswith(".db") and file != target_db:
#             files.append(file)
#     writer = Writer(target_db)
#     for db in files:
#         writer.combine_with(db)