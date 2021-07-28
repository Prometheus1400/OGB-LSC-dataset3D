import sys
sys.path.insert(1,'/home/ugrads/k/kaleb.dickerson2001/pcqm4m_dataset_with_coordinates/ogb')
from ogb.lsc import PygPCQM4MDataset
from ogb.utils import smiles2graph
import torch
from torch_geometric.data import DataLoader
import numpy as np
from tqdm import tqdm
import torch
import pickle

# for this to work, processes folder must be empty
# uses modified ogb code
dataset = PygPCQM4MDataset(
    root="/data3/kaleb.dickerson2001/Datasets/OGB-LSC-3D",
    smiles2graph=smiles2graph,
)

split_idx = {'train': list(np.arange(0, 2458702)), 'valid': list(np.arange(2458702, 2741027)), 'test': list(np.arange(2741027, 3023352))}
pickle.dump(split_idx, open( "/data3/kaleb.dickerson2001/Datasets/OGB-LSC-3D/pcqm4m_kddcup2021/split_idx3D.p", "wb" ))

batch_size = 256
train_loader = DataLoader(dataset[split_idx["train"]], batch_size=batch_size, shuffle=True)
valid_loader = DataLoader(dataset[split_idx["valid"]], batch_size=batch_size, shuffle=False)
test_loader = DataLoader(dataset[split_idx["test"]], batch_size=batch_size, shuffle=False)

# count = 0

# print(dic)
# print(count)
