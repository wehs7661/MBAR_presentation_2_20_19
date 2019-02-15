from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import mdtraj as md
from openmmtools.testsystems import SodiumChlorideCrystal
from pymbar import MBAR, timeseries

# User input
kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)
bin_range = range(2,100)
output_file_name = 'output/output.dat'
total_steps = 100000
print_frequency = 1 # Number of steps to skip when printing output

perform_simulation = False

if perform_simulation == True:
# Define a system
 system = SodiumChlorideCrystal()
# Define an integrator
 integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
# Define a simulation context
 simulation = Simulation(system.topology, system.system, integrator)
# Assign particle positions for this context
 simulation.context.setPositions(system.positions)
# Set the simulation type to energy minimization
 simulation.minimizeEnergy()
# Write simulation data every # steps
 simulation.reporters.append(StateDataReporter(output_file_name, print_frequency, step=True, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True))
# Run the simulation for 100000 steps
 simulation.step(total_steps)
# Use the Na-Cl potential energies as input to generate an ensemble of decorrelated samples
output_obj = open(output_file_name,'r')
U_temp = np.array([l.split(',')[3] for l in output_obj.readlines()])
output_obj.close()
U_all = np.array(np.delete(U_temp,0,0),dtype=float)
# Use pymbar's timeseries() function to analyze U_all, identify the equilibration time, and evaluate the statistical inefficiency for all points in the timeseries
[t0, g, Neff_max] = timeseries.detectEquilibration(U_all)
print(str(t0)+' '+str(g)+' '+str(Neff_max))
U_equil = U_all[t0:]
indices = timeseries.subsampleCorrelatedData(U_equil, g=g) # Determine indices of uncorrelated samples
print("Sub-sampled simulation data with "+str(len(indices))+" uncorrelated samples")
uncorrelated_data = energies[indices]
np.savetxt('output/uncorrelated_U.dat',uncorrelated_data)
# For each number of bins (thermodynamic states) in 'bin_range' we construct a matrix of weights, from which we will compute the MBAR-predicted potential energy differences and variance of the expectation value of the free energy
for nbins in bin_range:
# We'll compute the MBAR-predicted free energy differences and variance for a range of total thermodynamic states, nstates (nbins)

exit()
