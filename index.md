All data is hosted on the shared drive.

# Overview of PubChemQC
The OGB-LSC dataset is a based of the PubChemQC dataset (which is itself a subset of PubChem). OGB-LSC contains 3,803,453 molecules and PubChemQC contains 3,981,230.
The main contributions of the PubChemQC dataset are as follows:
* Molecular structures optimized by density functional theory (DFT)
* Calculated the excited states for over 2 million molecules using time-dependent DFT (TDDFT)
Note that PubChemQC does not provide SMILES or InChI representations directly, this information is in The PubChem Project

## OGB-LSC Dataset With 3D Information

There is a pytorch geometric compatible version of the OGB-LSC dataset augmented with 3D information located at **/mnt/dive/shared/kaleb/Datasets/OGB-LSC-3D**. This should already be processed and be compatible with your OGB package. Note that the number of features per node is now 12 (from 9). The x, y, and z coordinates are at the end of the feature vector.

## Smiles and Coordinates Database

In addition to this there is a large database located at **/mnt/dive/shared/kaleb/Datasets/PubChem3D** which contains smiles strings and corresponding coordinates. There is also a simple API in the same directory you can use to query this database with the smiles string and recieve the coodinates as an array of tuples. -1 indicates that smiles string was not found in database.

Here is an example how to use the API. 
```py
import sys
sys.path.insert(1,'/data3/kaleb.dickerson2001/Datasets/PubChem3D/db_api.py')
from db_api import Reader
reader = Reader('/data3/kaleb.dickerson2001/Datasets/PubChem3D/smile_coord.db')
coordinates = reader.get_coords("smiles string here")
```

