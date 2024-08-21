import argparse
import glob
import multiprocessing
import os
import re
import time

import numpy as np
import pandas as pd

from multiprocessing import Pool

from ccdc import io, conformer, docking, molecule, protein, screening
from ccdc.io import MoleculeReader, EntryReader, MoleculeWriter, EntryWriter
from ccdc.molecule import Molecule
from ccdc.docking import Docker
from ccdc.conformer import ConformerGenerator, ConformerSettings
from ccdc.screening import Screener
from ccdc.protein import Protein

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

class Analyse_Docking():
    def __init__(self):
        pass
    
    def grab_conf_files(self):
        # Define the regex pattern to search for
        pattern = r".*.conf"

        # Loop through all directories and subdirectories using os.walk and find files matching the regex pattern
        conf_files = []
        for root, dirs, files in os.walk("."):
            for file in files:
                if re.match(pattern, file):
                    conf_files.append(os.path.join(root, file))
                    
        return conf_files

    def write_csv(self, conf, results_dict):

        try:
            conf_base = conf.split("/")[-1].split(".")[0]

            settings = Docker.Settings.from_file(conf)
            results = Docker.Results(settings)
            ligands_scored = results.ligands

            poses = {}
            for num, score in enumerate(ligands_scored):
                myscore = score.scoring_term("fitness")
                poses[num+1] = myscore

                hbond_count = 0
                try:
                    hbonds = score.hbonds()
                    for bond in hbonds:
                        if bond.strength > 0.5 and bond.intermolecular is False:
                            hbond_count += 1
                            print(f"{conf} worked nice")
                except:
                    print(f"{conf} no Hbond detected")

                poses[num+1]["Hbond"] = hbond_count

            results_dict[conf_base] = poses
        except:
            print(f"{conf} run not completed")
    
    def run(self):
        num_cpu = multiprocessing.cpu_count()
        
        conf_files = self.grab_conf_files()

        with multiprocessing.Manager() as manager:
            results_dict = manager.dict()
            with multiprocessing.Pool(num_cpu) as pool:
                pool.starmap(self.write_csv, [(conf, results_dict) for conf in conf_files])

            results_dict = dict(results_dict)
        
        df = pd.DataFrame.from_dict({(i,j): results_dict[i][j] for i in results_dict.keys() for j in results_dict[i].keys()}, orient='index')

        return df