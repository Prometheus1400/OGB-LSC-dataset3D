import sqlite3
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
                            smile text,
                            coords blob
                            )""")
            self.conn.commit()
        except Exception:
            print("Connected to existing database with 'entries' table")

    def add_entry(self, smile, coords):
        self.c.execute(f"INSERT INTO entries VALUES ('{smile}','{coords}')")
        self.conn.commit()
    
    def add_entries(self, combined):
        self.c.executemany(f"INSERT INTO entries VALUES (?, ?)", combined)
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
    
    def get_coords(self, smile):
        """
        input:str smiles representation for molecule
        return:np.array[tuple] each tuple contains x,y,z coordinate for each atom
        """
        self.c.execute(f"SELECT coords FROM entries WHERE smile='{smile}'")

        try:
            coords = self.c.fetchone()[0]
            return coords
        except Exception:
            return -1

    def add_index(self, index_name, column):
        # this was most likely already performed
        print("adding indexes...")
        self.c.execute(f"CREATE INDEX {index_name} ON entries({column})")
        self.conn.commit()
        print("done")
    
    def __del__(self):
        self.conn.close()