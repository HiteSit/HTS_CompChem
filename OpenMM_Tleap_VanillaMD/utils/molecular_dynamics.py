import openmm
from openmm.app import *
from openmm import *
from openmm.unit import *

import argparse
import numpy as np
import sys
import time
from functools import wraps
from sys import stdout

# import MDAnalysis as mda
# from MDAnalysis.tests.datafiles import PDB, XTC, DCD

class Molecular_Dynamics:
    def __init__(self, prmtop_path, inpcrd_path, n_gpu:str):
        self.prmtop = AmberPrmtopFile(prmtop_path)
        self.inpcrd = AmberInpcrdFile(inpcrd_path)

        self.n_gpu = n_gpu
        self.platform = Platform.getPlatformByName('CUDA')
        self.proprieties = {'Precision': 'mixed', 'CudaDeviceIndex': self.n_gpu}

    def create_system(self):
        '''
        Generate the system
        '''
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    
    def choose_integrator_params(self, temp:int, delta:float):
        '''
        Choose the integrator
        '''
        self.integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picosecond, delta*picoseconds)
    
    def setup_simulation(self):
        '''
        Setup the simulation
        '''
        self.simulation = Simulation(
            self.prmtop.topology, self.system, self.integrator,
            platform=self.platform,
            platformProperties=self.proprieties
        )

        # Setup the initial position based and the PBC based on the inpcrd file
        self.simulation.context.setPositions(self.inpcrd.positions)
        if self.inpcrd.boxVectors is not None:
            self.simulation.context.setPeriodicBoxVectors(*self.inpcrd.boxVectors)
        
    def setup_positional_restraints(self, k:int):
        '''
        Setup the behave of the classical positional restraints baed on the formula:
        >> k*periodicdistance(x, y, z, x0, y0, z0)^2
        '''
        restraints = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        self.system.addForce(restraints)
        restraints.addGlobalParameter('k', k*kilojoules_per_mole/nanometer)

        restraints.addPerParticleParameter('x0')
        restraints.addPerParticleParameter('y0')
        restraints.addPerParticleParameter('z0')

        return restraints
    
    def remove_all_restraints(self):
        '''
        Remove the external force that generate the positional restraint
        '''
        for num, force in enumerate(self.system.getForces()):
            if isinstance(force, openmm.CustomExternalForce):
                self.system.removeForce(num)
    
    class Minimization:            
        def minimize(self):
            '''
            Minimize the system
            '''
            min_state = self.simulation.context.getState(getEnergy=True, getPositions=True)
            energy_before_min = min_state.getPotentialEnergy()
            print(f"Initial energy: {energy_before_min}")

            print("Beginning Minimization")
            self.simulation.minimizeEnergy()
            print("End Minimization")

            min_state = self.simulation.context.getState(getEnergy=True)
            energy_after_min = min_state.getPotentialEnergy()
            print(f"After Minimization: {energy_after_min}")
        
    class Nvt:
        def restraints_water(self):
            restraints = self.setup_positional_restraints(1000)
            for residue in self.prmtop.topology.residues():
                if residue.name not in ["HOH"]:
                    for atom in residue.atoms():
                        restraints.addParticle(atom.index, self.inpcrd.positions[atom.index])
        
        def setup_reporter(self, basename, nvt_steps, steps_saving, steps_log, append):
            nvt_state = self.simulation.context.getState(getEnergy=True)
            time_elapsed = nvt_state.getTime().value_in_unit(picoseconds)
            steps_elapsed = int(time_elapsed / 0.002)

            total_steps = nvt_steps + steps_elapsed
        
            self.simulation.reporters.append(DCDReporter(f"{basename}.dcd", steps_saving, append=append))

            self.simulation.reporters.append(StateDataReporter(stdout, steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))
            self.simulation.reporters.append(StateDataReporter(f"{basename}.log", steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))       

        def run(self, steps, temp):
            self.integrator.setTemperature(temp)
            self.simulation.step(steps)

    class Npt:
        def add_barostat(self):
            self.system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

        def restraint_backbone(self, restr_force, atoms_to_res:set):
            # Restraints the atoms
            restraint = self.setup_positional_restraints(restr_force)
            for atom in self.prmtop.topology.atoms():
                if atom.name in atoms_to_res:
                    restraint.addParticle(atom.index, self.inpcrd.positions[atom.index])
            
            # Restraints the residues
            for residue in self.prmtop.topology.residues():
                if residue.name in ["MOL"]:
                    for atom in residue.atoms():
                        restraint.addParticle(atom.index, self.inpcrd.positions[atom.index])
        
        def setup_reporter(self, basename, npt_steps, steps_saving, steps_log, append):
            npt_state = self.simulation.context.getState(getEnergy=True)
            time_elapsed = npt_state.getTime().value_in_unit(picoseconds)
            steps_elapsed = int(time_elapsed / 0.002)

            total_steps = npt_steps + steps_elapsed
        
            self.simulation.reporters.append(DCDReporter(f"{basename}.dcd", steps_saving, append=append))

            self.simulation.reporters.append(StateDataReporter(stdout, steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))
            self.simulation.reporters.append(StateDataReporter(f"{basename}.log", steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))

        def run(self, steps):
            self.simulation.step(steps)
    
    class Plain_Md():
        def setup_reporter(self, basename, md_steps, steps_saving, steps_log, append):
            md_state = self.simulation.context.getState(getEnergy=True)
            time_elapsed = md_state.getTime().value_in_unit(picoseconds)
            steps_elapsed = int(time_elapsed / 0.002)

            total_steps = md_steps + steps_elapsed
        
            self.simulation.reporters.append(DCDReporter(f"{basename}.dcd", steps_saving, append=append))

            self.simulation.reporters.append(StateDataReporter(stdout, steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))
            self.simulation.reporters.append(StateDataReporter(f"{basename}.log", steps_log, step=True,
                    potentialEnergy=True, temperature=True, volume=True, remainingTime=True, progress=True, speed=True, totalSteps=total_steps))  
            
            self.simulation.reporters.append(CheckpointReporter(basename + ".chk", steps_saving))

        def run(self, basename, steps):

            np.random.seed(33)
            self.simulation.context.setVelocitiesToTemperature(300*kelvin, np.random.randint(1, 1e6))
            self.simulation.step(steps)

            self.simulation.saveCheckpoint(f"{basename}.chk")
