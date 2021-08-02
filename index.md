#### Table of contents
1. [Overview of PubChemQC](#1)  
    1.1 [Getting Mol Object with RDkit](#1.1)  
    1.2 [Obtaining HOMO/LUMO Gap with cclib](#1.2)  
2. [Closing Comments](#2)

## Overview of PubChemQC  <a id="1"></a>
The OGB-LSC dataset is a based of the PubChemQC dataset (which is itself a subset of PubChem). OGB-LSC contains 3,803,453 molecules and PubChemQC contains 3,981,230. The PubChemQC dataset was generated using only CID (Chemical ID), InChI, and isomeric SMILES. All the calculations were performed from this information. Dataset located at **/mnt/dive/shared/kaleb/Datasets/PubChemQC**  
The main contributions of the PubChemQC dataset are as follows:  
* Molecular structures optimized by density functional theory (DFT)
* Calculated the excited states for over 2 million molecules using time-dependent DFT (TDDFT)  
Note that PubChemQC does not provide SMILES or InChI representations directly, this information is in The PubChem Project

### Getting Mol Object with RDkit <a id="1.1"></a>
PubChemQC has a .mol file for each molecule. This file can be conveniently read using RDkit to generate a molecule object with the 3D information like so:
```py
from rdkit import Chem

mol = Chem.MolFromMolFile(mol_file) #mol_file is path to .mol file
```
The atoms can then be iterated, and the coordinates obtained as follows:
```py
conf = mol.GetConformer()
sub = mol.GetSubstructMatch(mol)
for s in sub:
    atom = mol.GetAtoms()[s]
    atom_coords = list(conf.GetAtomPosition(s))))
```
However, for many of the .mol files, RDkit will generate an error in the following form "Explicit valence for atom \_\_\_\_ is greater than permitted". I am waiting on a response from the RDkit developers why this might be the case, and how to proceed with the calculations.  
I do still think RDkit is the way to go, because each atom is guaranteed to have to correct coordinates, compared to mapping the coordinates manually which is error prone as it is easy to assign coordinates to the wrong atoms of same atomic number.

### Obtaining HOMO/LUMO gap <a id="1.2"></a>
Each molecule directory contains a log.xz file. This can be used directly to obtain the HOMO-LUMO gap like so:
```py
from cclib import ccread #ccread requires open_babel to be installed
import lzma

with lzma.open(log_file, mode='rt') as file:
    data = ccread(file)
homo = data.homos[0]
energies = data.moenergies[0]
homolumogap = energies[homo+1] - energies[homo]
```
The log files are quite large, and parsing them with ccread is fairly slow.

## Closing Comments  <a id="2"></a>
Trying to map the coordinates from PubChemQC to the molecules in OGB-LSC is tricky, I think generating our own dataset may be a better option.  
On the other hand, reading and parsing millions of files from the shared drive is extremely slow, and there is still the issue with RDkit.

One last note is that I contacted the OGB team and they responded saying that they are "likely to release DFT-calculated 3D structure for training molecules in the next month or so".