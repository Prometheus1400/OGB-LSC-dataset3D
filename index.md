All data is hosted on the shared drive.

#### Table of contents
1. [Overview of PubChemQC](#1)  
    1.1 [Obtaining HOMO/LUMO Gap](#1.1)
2. [OGB-LSC Dataset With 3D Information](#2)
3. [Local Database](#3)

## Overview of PubChemQC  <a id="1"></a>
The OGB-LSC dataset is a based of the PubChemQC dataset (which is itself a subset of PubChem). OGB-LSC contains 3,803,453 molecules and PubChemQC contains 3,981,230. The PubChemQC dataset was generated using only CID (Chemical ID), InChI, and isomeric SMILES. All the calculations were performed from this information.
The main contributions of the PubChemQC dataset are as follows:  
* Molecular structures optimized by density functional theory (DFT)
* Calculated the excited states for over 2 million molecules using time-dependent DFT (TDDFT)
Note that PubChemQC does not provide SMILES or InChI representations directly, this information is in The PubChem Project

#### Obtaining HOMO/LUMO gap <a id="1.1"></a>
## OGB-LSC Dataset With 3D Information <a id="2"></a>

I have been working on extracting the DFT calculated 3D information from the PubChemQC dataset, and augmenting the OGB-LSC dataset with it. It is still a work and progress and it is currently only 80% the size of the original. Consider this dataset in 'beta' and I do not guarantee it's correctness yet.   
The dataset is pytorch geometric compatible and located at **/mnt/dive/shared/kaleb/Datasets/OGB-LSC-3D**. This should already be processed and be compatible with the OGB package. Note that the number of features per node is now 12 (from 9). The x, y, and z coordinates are at the end of the feature vector. There is a pickled split dictionary you can use instead of the built in one.  
Instead of:
```py
split_idx = dataset.get_idx_split()
```
Use:
```py
split_idx = pickle.load(open("/mnt/dive/shared/kaleb/Datasets/OGB-LSC-3D/pcqm4m_kddcup2021/split_idx3D.p", "rb" ))
train_loader = DataLoader(dataset[split_idx["train"]], batch_size=batch_size, shuffle=True)
valid_loader = DataLoader(dataset[split_idx["valid"]], batch_size=batch_size, shuffle=False)
test_loader = DataLoader(dataset[split_idx["test"]], batch_size=batch_size, shuffle=False)
```

## Local Database <a id="3"></a>

In addition to this there is a large database located at **/mnt/dive/shared/kaleb/Datasets/PubChem3D** which contains smiles strings and corresponding coordinates. There is also a simple API in the same directory you can use to query this database with the smiles string and recieve the coodinates as an array of tuples. -1 indicates that smiles string was not found in database.

Here is an example how to use the API. 
```py
import sys
sys.path.insert(1,'/data3/kaleb.dickerson2001/Datasets/PubChem3D/db_api.py')
from db_api import Reader
reader = Reader('/data3/kaleb.dickerson2001/Datasets/PubChem3D/smile_coord.db')
coordinates = reader.get_coords("smiles string here")
```

