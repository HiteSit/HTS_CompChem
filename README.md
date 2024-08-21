# Contents
- [Dock_Score_AdGPU](Dock_Score_AdGPU)
- [GOLD_Classic_Scaffold](GOLD_Classic_Scaffold)
- [GOLD_Reverse_Docking](GOLD_Reverse_Docking)

- [OpenMM_OpenFF_VanillaMD_MMPBSA](OpenMM_OpenFF_VanillaMD_MMPBSA)
- [OpenMM_Tleap_VanillaMD](OpenMM_Tleap_VanillaMD)

- [General_Utils](General_Utils)

## OpenMM Tleap
This repository includes a tleap script to generate the topology and parameter files for OpenMM simulations.

The script initiates with a call to tleap to prepare the topology file. This process is fully automated, ensuring the topology, the simulation box size, and the ion concentration are smoothly established. The forcefield applied to the macromolecule can be customized within the tleap section of the script, as well as that for the ligand.

After topology generation, the script facilitates a classical molecular dynamics simulation using a Langevin integrator. The equilibration phase consists of both NVT (constant number of particles, volume, and temperature) and NPT (constant number of particles, pressure, and temperature) ensembles. During the NVT phase, restraints are automatically applied to the solute (excluding the ions) to stabilize the water molecules. For the NPT phase, restraints can be specified through a dictionary, allowing for selective constraints on either macromolecular atoms (e.g., the backbone) or entire residues, such as the ligand.

Finally, a standard molecular dynamics simulation is conducted at 300K.

## OpenMM OpenFF VanillaMD MMPBSA

Give a look to the [example_notebook.ipynb](OpenMM_OpenFF_VanillaMD_MMPBSA%2Fexamples%2Fexample_notebook.ipynb)

## GeneralUtils
This Python library provides a set of tools for the preparation, protonation, and standardization of small molecules. It combines the functionality of several cheminformatics libraries, including RDKit, OpenEye, and DataMol, to facilitate molecular preprocessing for computational chemistry and molecular modeling tasks.

### Dependencies
- `mamba install datamol`
- `mamba install rdkit`
- `mamba install -c openeye openeye-toolkits`
