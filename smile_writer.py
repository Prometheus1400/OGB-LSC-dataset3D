from database_api import Reader, Writer
import numpy as np
from utils import split_into, ThreadWithReturn, CIDs2smiles
import time

reader = Reader("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
total = reader.get_all_CIDs()

# for cid in total:
#     res = CIDs2smiles([cid])

#     if len(res) != 0:
#         print(res[0][1], end=", ")

# smiles_with_CIDS = CIDs2smiles([71169753, 71174498, 54558593, 54561827, 11954886, 11957815, 11958203, 11961360, 11969494, 11969616, 11971926, 11897997, 40473929, 44517949, 44518964, 11948253, 11344370, 11344865, 25234587, 25246089, 25248327, 11858085, 11862664, 40479360, 40479368, 40480217, 40480612, 40488388, 10654344, 10654365, 25220786, 25221190, 25224339, 71177452, 71178201, 71179330, 71182576, 71187016, 71187826, 71188822])
# writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
# writer.add_smiles(smiles_with_CIDS)


del(reader)
if __name__ == "__main__":
    NUM_THREADS = 5
    total_split = split_into(total, NUM_THREADS)

    T = []
    for i in range(NUM_THREADS):
        T.append(ThreadWithReturn(target=CIDs2smiles, args=(total_split[i],)))
        T[i].start()
    
    arr = np.empty(NUM_THREADS, dtype=object)
    for i,t in enumerate(T):
        arr[i] = t.join()

    smiles_with_CIDS = np.concatenate(arr)
    
    time.sleep(30)
    writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
    writer.add_smiles(smiles_with_CIDS, type='canonical')