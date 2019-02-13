from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

gro = GromacsGroFile('input/input.gro')
top = GromacsTopFile('input/input.top',includeDir='/home/gmeek/software/anaconda3/share/gromacs/top')
system = top.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.5*nanometers)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output/output.pdb', 1))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
