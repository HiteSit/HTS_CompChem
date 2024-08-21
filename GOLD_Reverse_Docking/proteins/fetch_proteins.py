import glob
import os

from pymol import cmd
def fetch_proteins(proteins_path, slice):

    os.chdir(proteins_path)

    txt = glob.glob("*.txt")

    with open(txt[0]) as f:
        pdb = f.read()

    pdb = [x.strip() for x in pdb.split(",")]
    pdb_less = pdb[::slice]

    count_removed = 0
    for p in pdb_less:
        try:
            cmd.fetch(p, type="pdb", state=1)
            cmd.select("To_Remove", "all & ! (polymer.protein | organic)")
        except:
            print(f"Impossible to fetch {p}")

        if cmd.select("Probas", "organic") == 0:
            print(f"{p} has no organic")
            cmd.remove("To_Remove")
            cmd.save(filename=p.lower() + "_cavity.pdb", format="pdb")
            cmd.delete("all")

            count_removed += 1
            continue
        else:
            cmd.remove("To_Remove")
            cmd.save(filename=p + ".pdb", format="pdb")
            cmd.delete("all")

    os.chdir("../")

    print(f"Out of: {len(pdb_less)} only {len(pdb_less) - count_removed} had cristal ligand")