import sqlite3
import utils
import numpy as np

class Writer:
    def __init__(self, path):
        self.path = path
        self.conn = sqlite3.connect(path)
        self.c = self.conn.cursor()

        self.make_table()

    def make_table(self):
        try:
            self.c.execute("""CREATE TABLE entries (
                            CID integer,
                            coords text
                            )""")
            self.conn.commit()
            print("Created new table 'entries'")
        except Exception:
            print("Connected to existing database with 'entries' table")

    def add_entry(self, smile, coords):
        self.c.execute(f"INSERT INTO entries VALUES ('{smile}','{coords}')")
        self.conn.commit()
    
    def add_entries(self, combined):
        self.c.executemany(f"INSERT INTO entries VALUES (?, ?)", combined)
        self.conn.commit()
    
    def combine_with(self, other_db):
        self.c.execute(f"attach '{other_db}' as toMerge")
        self.c.execute("BEGIN")
        self.c.execute("insert into entries select * from toMerge.entries")
        self.c.execute("COMMIT")
        self.c.execute("detach toMerge")
    
    def add_smiles(self,combined, type='isomeric'):
        #combined is (new smiles, old smiles)
        self.c.executemany(f"UPDATE entries SET {type}_SMILES = ? WHERE {type}_SMILES = ?", combined)
        self.conn.commit()
    
    def clear(self):
        self.c.execute("DELETE FROM entries")
        self.conn.commit()

    def __del__(self):
        self.conn.close()

class Reader:
    def __init__(self, path):
        self.path = path
        try:
            self.conn = sqlite3.connect(path)
            self.c = self.conn.cursor()
            print("Connected to existing database with 'entries' table")
        except Exception:
            print("database doesnt exist at this path")
    
    def get_coords(self, CID):
        """
        input:str smiles representation for molecule
        return:np.array[tuple] each tuple contains x,y,z coordinate for each atom
        """
        self.c.execute(f"SELECT coords FROM entries WHERE CID='{CID}'")
        try:
            coords = self.c.fetchone()[0]
            coords = coords.split()

            num_splits = len(coords) // 3
            coords = utils.split_into(coords, num_splits)
            return coords
        except Exception:
            print(f"CID : {CID} NOT PRESENT OR COORDS NOT AVAILABLE")
            return -1

    def get_all_CIDs(self):
        self.c.execute("SELECT CID FROM entries WHERE canonical_SMILES IS NULL")
        total = self.c.fetchall()

        to_return = np.empty(len(total), dtype=int)
        for i,tup in enumerate(total):
            to_return[i] = tup[0]

        return to_return
    
    def get_all_SMILES(self):
        self.c.execute("SELECT isomeric_SMILES FROM entries WHERE isomeric_SMILES NOT NULL")
        total = self.c.fetchall()

        to_return = []
        for i,tup in enumerate(total):
            to_return.append(tup[0])

        return to_return
    
    def __del__(self):
        self.conn.close()