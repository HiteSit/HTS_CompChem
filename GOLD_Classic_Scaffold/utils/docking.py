import argparse
import glob
import multiprocessing
import os
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

class Prepare_Files():
    def __init__(self, protein, ligands_to_dock, crystal_ligand):
        self.protein = os.path.abspath(protein)
        self.protein_basename = os.path.basename(self.protein).split('.')[0]
        
        self.crystal_ligand = os.path.abspath(crystal_ligand)
        
        self.ligands_to_dock = os.path.abspath(ligands_to_dock)
        self.ligands_prep_folder = "ligands_prep"
        self.output_folder = "output"
        self.output_folder_merged = "output_merged"
    
    def work_dir(self, workdir):
        self.workdir = workdir
        os.chdir(self.workdir)
        
    def startup_folders(self):
        os.makedirs(self.ligands_prep_folder, exist_ok=True)
        os.makedirs(self.output_folder, exist_ok=True)
        
    def prepare_ligands(self):
        preparator = Docker.LigandPreparation()
        preparator.settings.add_hydrogens = True
        
        ligands_to_prep = EntryReader(self.ligands_to_dock)
        
        ligands_basename = []
        ligands_prep_path = []

        for lig in ligands_to_prep:
            temp = lig.molecule.identifier
            ligands_basename.append(temp)
            
            ligand_path = os.path.join(self.ligands_prep_folder, temp + "_prep.mol2")
            ligand_path = os.path.abspath(ligand_path)
            ligands_prep_path.append(ligand_path)
            
            with MoleculeWriter(ligand_path) as W:
                prepared_lig = preparator.prepare(lig)
                W.write(prepared_lig.molecule)
                
        self.ligands_basename = ligands_basename
        self.ligands_prep_path = ligands_prep_path
        
    def prepare_protein(self):
        prot_entry = Protein.from_file(self.protein)

        prot_entry.remove_all_waters()
        prot_entry.remove_unknown_atoms()
        prot_entry.add_hydrogens()
        lig_prot = prot_entry.ligands
        prot_prepared_filename = self.protein_basename + "_clean.pdb"

        for l in lig_prot:
            prot_entry.remove_ligand(l.identifier)
        with EntryWriter(prot_prepared_filename) as writer:
            writer.write(prot_entry)

        self.prot_abs = os.path.abspath(prot_prepared_filename)

class Classical_Docking(Prepare_Files):
    def __init__(self, protein, ligands_to_dock, crystal_ligand):
        super().__init__(protein, ligands_to_dock, crystal_ligand)
    
    def prepare_docker(self, ligand_prepared):
        docker = Docker()
        settings = docker.settings

        settings.fitness_function = 'goldscore'
        settings.rescore_function = "chemscore"
        settings.autoscale = 100.
        settings.early_termination = False
        settings.write_options = ["NO_LINK_FILES"]

        settings.diverse_solutions = (True, 10, 1.5)
        settings.save_lone_pairs = False

        native_ligand = MoleculeReader(self.crystal_ligand)[0]

        settings.add_protein_file(self.prot_abs)
        settings.add_ligand_file(ligand_prepared, 10)
        settings.reference_ligand_file = self.crystal_ligand
        settings.binding_site = settings.BindingSiteFromLigand(self.prot_abs, native_ligand, 10.0)

        return docker
    
    def write_conf(self, params):
        ligand_prepared = params[0]
        basename = params[1]

        my_docker = self.prepare_docker(ligand_prepared)
        settings = my_docker.settings

        lig_basename = basename.split(".")[0]

        settings.output_directory = os.path.join(lig_basename)
        settings.output_file = lig_basename + "_merged.mol2"

        my_docker.dock(f"output/{lig_basename}.conf")
        
    def run_docking(self, n_cpu):
        params = [(ligand_prepared, basename) for ligand_prepared, basename in zip(self.ligands_prep_path, self.ligands_basename)]
        with Pool(n_cpu) as p:
            p.map(self.write_conf, params)
            
    def extract_files(self):
        import shutil
        os.makedirs(self.output_folder_merged, exist_ok=True)
        
        to_move = []
        for root, dirs, files in os.walk(self.output_folder):
            for file in files:
                if file.endswith("_merged.mol2") and not file.endswith("pts_merged.mol2"):
                    to_move.append(os.path.join(root, file))
                    
        for file in to_move:
            shutil.copy2(file, self.output_folder_merged)

class Scaffold_Docking(Prepare_Files):
    def __init__(self, protein, ligands_to_dock, crystal_ligand, scaffold):
        super().__init__(protein, ligands_to_dock, crystal_ligand)
        self.scaffold = os.path.abspath(scaffold)
    
    def prepare_docker(self, ligand_prepared):
        docker = Docker()
        settings = docker.settings

        settings.fitness_function = 'goldscore'
        settings.rescore_function = "chemscore"
        settings.autoscale = 100.
        settings.early_termination = False
        settings.write_options = ["NO_LINK_FILES"]

        settings.diverse_solutions = (True, 10, 1.5)
        settings.save_lone_pairs = False
        
        scaffold = MoleculeReader(self.scaffold)[0]
        native_ligand = MoleculeReader(self.crystal_ligand)[0]

        settings.add_protein_file(self.prot_abs)
        settings.add_ligand_file(ligand_prepared, 10)
        settings.reference_ligand_file = self.crystal_ligand
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))
        settings.binding_site = settings.BindingSiteFromLigand(self.prot_abs, native_ligand, 10.0)

        return docker
    
    def write_conf(self, params):
        ligand_prepared = params[0]
        basename = params[1]

        my_docker = self.prepare_docker(ligand_prepared)
        settings = my_docker.settings

        lig_basename = basename.split(".")[0]

        settings.output_directory = os.path.join(lig_basename)
        settings.output_file = lig_basename + "_merged.mol2"

        my_docker.dock(f"output/{lig_basename}.conf")
        
    def run_docking(self, n_cpu):
        params = [(ligand_prepared, basename) for ligand_prepared, basename in zip(self.ligands_prep_path, self.ligands_basename)]
        with Pool(n_cpu) as p:
            p.map(self.write_conf, params)
            
    def extract_files(self):
        import shutil
        os.makedirs(self.output_folder_merged, exist_ok=True)
        
        to_move = []
        for root, dirs, files in os.walk(self.output_folder):
            for file in files:
                if file.endswith("_merged.mol2") and not file.endswith("pts_merged.mol2"):
                    to_move.append(os.path.join(root, file))
                    
        for file in to_move:
            shutil.copy2(file, self.output_folder_merged)