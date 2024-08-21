from pymol import cmd
import os, re, shutil, glob
import subprocess
import argparse

def execute_bash(command):
    # Run the command using subprocess
    try:
        subprocess.run(command, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        
class Tleap_Preparation:
    def __init__(self, prt_file, lig_file):
        self.prt_abs = os.path.abspath(prt_file)
        self.prt_name = os.path.basename(self.prt_abs)
        self.prt_basename = self.prt_name.split(".")[0]
        
        self.lig_file = os.path.abspath(lig_file)
        self.lig_name = os.path.basename(self.lig_file)
        self.lig_basename = self.lig_name.split(".")[0]
    
    def work_dir(self, work_dir):
        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)
    
    def pdb2pqr(self):
        command = ["pdb2pqr", self.prt_abs, self.prt_basename + "_H.pdb", "--with-ph", "7.2", "--ff", "AMBER", "--ffout", "AMBER", "--titration-state-method", "propka"]
        execute_bash(command)
    
    def obabel(self):
        self.ligand_obabel_sdf = self.lig_basename + "_obabel.sdf"
        self.ligand_obabel_pdb = self.lig_basename + "_obabel.pdb"
        
        command = ["obabel", self.lig_file, "-O", self.ligand_obabel_sdf, "-p"]
        execute_bash(command)
        command = ["obabel", self.lig_file, "-O", self.ligand_obabel_pdb, "-p"]
        execute_bash(command)
    
    def rdkit_sanitize(self):
        from rdkit import Chem
        
        self.ligand_rdkit_sdf = self.lig_basename + "_rdkit.sdf"
        
        sdf_suppl = Chem.SDMolSupplier(self.ligand_obabel_sdf)
        sdf_writer = Chem.SDWriter(self.ligand_rdkit_sdf)
        
        for mol in sdf_suppl:
            try:
                Chem.SanitizeMol(mol)
                mol = Chem.AddHs(mol, addCoords=True)
                
                formal_charge = Chem.GetFormalCharge(mol)
                
                sdf_writer.write(mol)
            except:
                pass
        
        self.formal_charge = str(formal_charge)
        print(self.formal_charge)
    
    def create_complex(self):
        cmd.reinitialize()
        cmd.load(self.prt_basename + "_H.pdb")
        cmd.load(self.lig_basename + "_obabel.pdb")
        cmd.create("complex", f"{self.prt_basename}_H or {self.lig_basename}_obabel")
        cmd.save("Init.pdb", "complex")
        
    def antechamber_ligand(self):
        self.ligand_amber_mol2 = self.lig_basename + "_amber.mol2"
        self.ligand_amber_frcmod = self.lig_basename + "_amber.frcmod"
        
        command_antechamber = ["antechamber", "-i", self.ligand_rdkit_sdf, "-fi", "sdf", "-o", self.ligand_amber_mol2, "-fo", "mol2", "-c", "bcc", "-at", "gaff", "-nc", self.formal_charge]
        execute_bash(command_antechamber)
        
        command_parmchk2 = ["parmchk2", "-i", self.ligand_amber_mol2, "-f", "mol2", "-o", self.ligand_amber_frcmod, "-a", "Y"]
        execute_bash(command_parmchk2)
    
    def run_tleap_box(self, padding):
        self.padding = padding
        
        tleap_in = f'''
        source leaprc.protein.ff14SB
        source leaprc.water.tip3p
        source leaprc.gaff
        loadamberparams frcmod.ionsjc_tip3p

        loadamberparams {self.ligand_amber_frcmod}

        lig = loadmol2 {self.ligand_amber_mol2}
        prot = loadpdb {self.prt_basename}_H.pdb
        sys = combine {{lig prot}}
        set default PBRadii mbondi2

        solvateBox sys TIP3PBOX {self.padding}
        charge sys
        savepdb sys Init_Wat.pdb

        quit

        '''

        with open("tleap.in", "w") as f:
            f.write(tleap_in)

        command_tleap = ["tleap", "-f", "tleap.in"]
        execute_bash(command_tleap)
    
    def run_tleap_ions(self, conc):
        def calc_num_ions(conc):
            Co = conc
            
            with open("leap.log") as f:
                data = f.read()
            
            Nw = float(re.findall(r'Added (\d+) residues', data)[0])
            Q = float(re.findall(r'Total perturbed charge:\s*([+-]?\d+(?:\.\d+)?)', data)[0])
            
            No = (Nw * Co) / 56
            Nplus = round(No - (Q / 2))
            Nminus = round(No + (Q / 2))
            
            return Nplus, Nminus
        
        Nplus, Nminus = calc_num_ions(conc)
        
        tleap_in = f'''
        source leaprc.protein.ff14SB
        source leaprc.water.tip3p
        source leaprc.gaff
        loadamberparams frcmod.ionsjc_tip3p

        loadamberparams {self.ligand_amber_frcmod}

        lig = loadmol2 {self.ligand_amber_mol2}
        prot = loadpdb {self.prt_basename}_H.pdb
        sys = combine {{lig prot}}
        set default PBRadii mbondi2

        solvateBox sys TIP3PBOX {self.padding}
        charge sys
        savepdb sys Init_Wat.pdb
        
        addIonsRand sys Na+ {Nplus}
        addIonsRand sys Cl- {Nminus}
        savepdb sys Init_Ions.pdb
        
        saveamberparm sys system.prmtop system.inpcrd
        quit

        '''

        with open("tleap.in", "w") as f:
            f.write(tleap_in)

        command_tleap = ["tleap", "-f", "tleap.in"]
        execute_bash(command_tleap)