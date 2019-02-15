from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import mdtraj as md

# Setup and run a Gromacs simulation in OpenMM
gro = GromacsGroFile('input/input.gro')
top = GromacsTopFile('input/input.top',includeDir='/home/gmeek/software/anaconda3/share/gromacs/top')
system = top.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.5*nanometers)
force = CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output/output.pdb', 1))
simulation.reporters.append(StateDataReporter('output/energies.dat', 1, step=True, potentialEnergy=True, temperature=True))
simulation.step(100000)

# Analyze the trajectory output using MDTraj

traj = md.load('output/output.pdb')

# Calculate the Na-Cl stretching distance

distances = md.compute_distances(traj,np.ndarray((1,2)))

# Calculate the potential energy



# Write the distances to file

np.savetxt('output/distances.dat',distances)


