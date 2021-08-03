import os
import numpy as np
from database_api import Writer
from tqdm import tqdm, trange
from utils import get_coordinates, split_into, ThreadWithReturn
from multiprocessing import Pool
from random import shuffle
from rdkit import Chem
import lzma
from cclib.io import ccread
from torch_geometric.data import InMemoryDataset
from torch_geometric.data import Data
import torch
import os.path as osp
import shutil

allowable_features = {
    'possible_atomic_num_list' : list(range(1, 119)) + ['misc'],
    'possible_chirality_list' : [
        'CHI_UNSPECIFIED',
        'CHI_TETRAHEDRAL_CW',
        'CHI_TETRAHEDRAL_CCW',
        'CHI_OTHER'
    ],
    'possible_degree_list' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'misc'],
    'possible_formal_charge_list' : [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 'misc'],
    'possible_numH_list' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 'misc'],
    'possible_number_radical_e_list': [0, 1, 2, 3, 4, 'misc'],
    'possible_hybridization_list' : [
        'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'misc'
        ],
    'possible_is_aromatic_list': [False, True],
    'possible_is_in_ring_list': [False, True],
    'possible_bond_type_list' : [
        'SINGLE',
        'DOUBLE',
        'TRIPLE',
        'AROMATIC',
        'misc'
    ],
    'possible_bond_stereo_list': [
        'STEREONONE',
        'STEREOZ',
        'STEREOE',
        'STEREOCIS',
        'STEREOTRANS',
        'STEREOANY',
    ], 
    'possible_is_conjugated_list': [False, True],
}

def safe_index(l, e):
    """
    Return index of element e in list l. If e is not present, return the last index
    """
    try:
        return l.index(e)
    except:
        return len(l) - 1

def atom_to_feature_vector(atom, coord):
    """
    Converts rdkit atom object to feature list of indices
    :param mol: rdkit atom object
    :return: list
    """
    atom_feature = [
            safe_index(allowable_features['possible_atomic_num_list'], atom.GetAtomicNum()),
            allowable_features['possible_chirality_list'].index(str(atom.GetChiralTag())),
            safe_index(allowable_features['possible_degree_list'], atom.GetTotalDegree()),
            safe_index(allowable_features['possible_formal_charge_list'], atom.GetFormalCharge()),
            safe_index(allowable_features['possible_numH_list'], atom.GetTotalNumHs()),
            safe_index(allowable_features['possible_number_radical_e_list'], atom.GetNumRadicalElectrons()),
            safe_index(allowable_features['possible_hybridization_list'], str(atom.GetHybridization())),
            allowable_features['possible_is_aromatic_list'].index(atom.GetIsAromatic()),
            allowable_features['possible_is_in_ring_list'].index(atom.IsInRing()),
            coord[0],   #X coord
            coord[1], #Y coord
            coord[2]   #Z coord
            ]
    return atom_feature

def bond_to_feature_vector(bond):
    """
    Converts rdkit bond object to feature list of indices
    :param mol: rdkit bond object
    :return: list
    """
    bond_feature = [
                safe_index(allowable_features['possible_bond_type_list'], str(bond.GetBondType())),
                allowable_features['possible_bond_stereo_list'].index(str(bond.GetStereo())),
                allowable_features['possible_is_conjugated_list'].index(bond.GetIsConjugated()),
            ]
    return bond_feature

def getHOMOLUMO(log_file):
    """
    input: molecule compressed log file
    output: homo-lumo gap as float
    """
    try:
        with lzma.open(log_file, mode='rt') as file:
            data = ccread(file)
        homo = data.homos[0]
        energies = data.moenergies[0]
        homolumogap = energies[homo+1] - energies[homo]
        return homolumogap
    except Exception:
        return None

def mol2graph(mol):
    """
    Converts rdki Mol object to graph Data object
    :input: Mol object (str)
    :return: graph object
    """
    conf = mol.GetConformer()
    sub = mol.GetSubstructMatch(mol)
    atom_features = []
    for s in sub:
        atom_features.append(atom_to_feature_vector(mol.GetAtoms()[s], list(conf.GetAtomPosition(s))))

    # sys.exit()
    #changed to float64 from int64
    x = np.array(atom_features, dtype = np.float64)

    # bonds
    num_bond_features = 3  # bond type, bond stereo, is_conjugated
    if len(mol.GetBonds()) > 0: # mol has bonds
        edges_list = []
        edge_features_list = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            edge_feature = bond_to_feature_vector(bond)

            # add edges in both directions
            edges_list.append((i, j))
            edge_features_list.append(edge_feature)
            edges_list.append((j, i))
            edge_features_list.append(edge_feature)

        # data.edge_index: Graph connectivity in COO format with shape [2, num_edges]
        edge_index = np.array(edges_list, dtype = np.int64).T

        # data.edge_attr: Edge feature matrix with shape [num_edges, num_edge_features]
        edge_attr = np.array(edge_features_list, dtype = np.int64)

    else:   # mol has no bonds
        edge_index = np.empty((2, 0), dtype = np.int64)
        edge_attr = np.empty((0, num_bond_features), dtype = np.int64)

    graph = dict()
    graph['edge_index'] = edge_index
    graph['edge_feat'] = edge_attr
    graph['node_feat'] = x
    graph['num_nodes'] = len(x)

    return graph

def make_dataset(dir_list):
    if len(dir_list) == 0:
        return

    data_list = []
    for directory in tqdm(dir_list):
        start, end = int(directory[9:18]), int(directory[19:28])
        start_pad = str(start).zfill(9)
        end_pad = str(end).zfill(9)
        path_to_mols = "".join([path2tars, f"/Compound_{start_pad}_{end_pad}"])

        count = 0
        for i in range(start, end+1):
            i_pad = str(i).zfill(9)

            # path to specific molecule directory
            mol_dir_path = "".join([path_to_mols, f"/{i_pad}"])

            # checks if molecule number directory exists
            exists = os.path.isdir(mol_dir_path)
            if exists:
                mol_file = "".join([mol_dir_path, "/", i_pad, ".mol"])
                log_file = "".join([mol_dir_path, "/", i_pad, ".b3lyp_6-31g(d).log.xz"])
                if os.path.exists(mol_file) and os.path.exists(log_file):
                    mol = Chem.MolFromMolFile(mol_file)
                    if mol != None:
                        graph = mol2graph(mol)
                        homolumogap = getHOMOLUMO(log_file)
                        if homolumogap != None:
                            assert(len(graph['edge_feat']) == graph['edge_index'].shape[1])
                            assert(len(graph['node_feat']) == graph['num_nodes'])

                            data = Data()
                            data.__num_nodes__ = int(graph['num_nodes'])
                            data.edge_index = torch.from_numpy(graph['edge_index']).to(torch.int64)
                            data.edge_attr = torch.from_numpy(graph['edge_feat']).to(torch.int64)
                            data.x = torch.from_numpy(graph['node_feat']).to(torch.float64)
                            data.y = torch.Tensor([homolumogap])
                            data_list.append(data)
    return data_list

def make_threads(list_obj):
    if __name__ == "__main__":
        NUM_THREADS = 8
        _dir_list_fragmented = split_into(list_obj, NUM_THREADS)
        T = []
        for i in range(NUM_THREADS):
            T.append(ThreadWithReturn(target=make_dataset, args=(_dir_list_fragmented[i],)))
            T[i].start()
        
        arr = np.empty(NUM_THREADS, dtype=object)
        for i,t in enumerate(T):
            arr[i] = t.join()
        _data_list = []
        for lists in arr:
            for data_obj in lists:
                _data_list.append(data_obj)
        return _data_list

class MyDataset(InMemoryDataset):
    def __init__(self, root, path2tars, make_dataset, transform=None, pre_transform=None):
        self.folder = root + "/OGB-LSC3D"
        self.path2tars = path2tars
        self.make_dataset = make_dataset
        super(MyDataset, self).__init__(self.folder, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return "dont_download.txt"
    @property
    def processed_file_names(self):
        return "geometric_data_processed.pt"
    def download(self):
        print("trying to download")

    def process(self):
        print("trying to process")
        dirnames = next(os.walk(self.path2tars), (None, [], None))[1]
        dirnames.sort()
        shuffle(dirnames)

        if __name__ == "__main__":
            NUM_PROCESSES = 32
            dir_list_fragmented = split_into(dirnames, NUM_PROCESSES)
            # T = []
            # for i in range(NUM_PROCESSES):
            #     T.append(ThreadWithReturn(target=make_dataset, args=(dir_list_fragmented[i],)))
            #     T[i].start()
            
            # arr = np.empty(NUM_PROCESSES, dtype=object)
            # for i,t in enumerate(T):
            #     arr[i] = t.join()
            # data_lists = arr
            with Pool() as pool:
                    data_lists = pool.map(make_threads, dir_list_fragmented)

            data_list = []
            for lists in data_lists:
                for data_obj in lists:
                    data_list.append(data_obj)
            data, slices = self.collate(data_list)
            print('Saving...')
            torch.save((data, slices), self.processed_paths[0])
            print('processing complete')

path2tars = "/mnt/dive/shared/kaleb/Datasets/PubChemQC"
dataset = MyDataset(root="/home/ugrads/k/kaleb.dickerson2001/OGB-LSC-dataset3D",
                    path2tars=path2tars,
                    make_dataset=make_dataset)
                    