import sqlite3

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

        # try:
        #     coords = self.c.fetchone()[0]
        #     return coords
        # except Exception:
        #     return -1
        return self.c.fetchall()
    
    def __del__(self):
        self.conn.close()