# This Python script contains example NaCl simulations in OpenMM, 
# followed by analysis with MBAR, and WHAM.

# The protocol is broken up into the following sections:

# 1) Provide user input
# 2) Prepare and run OpenMM simulation(s)
# 3) Sub-sample simulation data to identify uncorrelated samples
# 4) Optimize 'weights' for each thermodynamic state using MBAR
# 5) Calculate free energy surface with MBAR using Langevin dynamics simulation data
# 6) Calculate PMF with MBAR using Langevin dynamics simulation data
# 7) Calculate PMF with MBAR after umbrella sampling of Langevin dynamics simulation data
# 8) Calculate PMF with WHAM after umbrella sampling of Langevin dynamics simulation data


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import mdtraj as md
from openmmtools.testsystems import SodiumChlorideCrystal
from pymbar import MBAR, timeseries

#############
#
# 1) Provide user input
#
#############

simulation_time_step = 0.002 # Units = picoseconds
simulation_temperature = 300 # Units = Kelvin
kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)
state_range = range(2,20)
output_file_name = 'NaCl/output/output.dat'
num_simulations = 1
simulation_steps = 10000
print_frequency = 1 # Number of steps to skip when printing output
nskip = 10 # Number of steps to skip when reading timeseries data to find the equilibration time
data_exists = True

#############
#
# 2) Prepare and run OpenMM simulation
#
#############

# An OpenMM simulation requires three input objects: a 'system', an 'integrator', and a 'context'

total_simulation_time = simulation_time_step * simulation_steps # Units = picoseconds
if data_exists == False:
 system = SodiumChlorideCrystal() # Define a system
 integrator = LangevinIntegrator(simulation_temperature, total_simulation_time, simulation_time_step) # Define an integrator
 simulation = Simulation(system.topology, system.system, integrator) # Define a simulation 'context'
 simulation.reporters.append(StateDataReporter(output_file_name, print_frequency, \
 step=True, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True))
# If we are running more than one simulation, we randomize the initial coordinates for 'Cl':
#if num_simulations == 1:
 simulation.context.setPositions(system.positions) # Assign particle positions for this context
#else:
# simulation.context.setPositions([[0.,0.,0.],[]]
 simulation.minimizeEnergy() # Set the simulation type to energy minimization
 print("Performing openMM simulation for "+str(simulation_steps)+" steps.")
 simulation.step(simulation_steps) # Run the simulation
if data_exists == True:
 print("Reading existing simulation output from: "+output_file_name)
#############
#
# 3) Sub-sample simulation data to identify uncorrelated samples
#
#############

output_obj = open(output_file_name,'r')
# E_total_all stores the total energies from NaCl simulation output
E_total_all_temp = np.array([l.split(',')[3] for l in output_obj.readlines()])
output_obj.close()
output_obj = open(output_file_name,'r')
# T_all stores the Temperatures from NaCl simulation output
T_all_temp = np.array([l.split(',')[4] for l in output_obj.readlines()])
output_obj.close()
E_total_all = np.array(np.delete(E_total_all_temp,0,0),dtype=float)
T_all = np.array(np.delete(T_all_temp,0,0),dtype=float)
# Use pymbar's timeseries() function to analyze U_all, identify the equilibration time (t0), and evaluate the statistical inefficiency (g) for all points in the timeseries
# THIS STEP IS TIME-CONSUMING.  TO SPEED UP INCREASE 'nskip'
[t0, g, Neff_max] = timeseries.detectEquilibration(E_total_all,nskip=nskip)
print("The total simulation time was "+str(total_simulation_time)+" picoseconds.")
print("The equilibration time was "+str(t0*simulation_time_step)+" picoseconds.")
E_total_equil = E_total_all[t0:]
T_equil = T_all[t0:]
uncorrelated_energies = timeseries.subsampleCorrelatedData(E_total_equil, g=g) # Determine indices of uncorrelated samples
print("Sub-sampled simulation data with "+str(len(uncorrelated_energies))+" uncorrelated samples")
U_uncorrelated = E_total_equil[uncorrelated_energies] # Uncorrelated total energies 
np.savetxt('NaCl/output/uncorrelated_total_energies.dat',uncorrelated_energies)
T_uncorrelated = T_equil[uncorrelated_energies] # Uncorrelated temperatures
# Calculate the reduced potential energies for the uncorrelated samples
U_reduced = np.array([U_uncorrelated[index]/(T_uncorrelated[index]*kB) for index in range(0,len(uncorrelated_energies))])

#############
#
# 4) Optimize 'weights' for each thermodynamic state using MBAR
#
#############

# Now we are ready to use U_reduced to optimize a set of weights for each thermodynamic 'state', 
# We split the full range of temperatures sampled during the simulation into bins (one 'bin' for each 'state')
# Bin widths are calculated as: ( T_max - T_min ) / nbins
T_max = max(T_uncorrelated)
T_min = min(T_uncorrelated)
T_ranges_for_each_num_states = np.array([])
distributions_for_each_num_states = np.array([])
free_energies_for_each_num_states = np.array([])
weights_for_each_num_states = np.array([])
for num_states in state_range:
# First we 'bin' (count the total number of samples) in each thermodynamic state
 T_step_size = (T_max - T_min) / num_states
 print("The maximum and minimum sampled temperatures were: "+str(T_max)+" and "+str(T_min)+", respectively.")
 print("Binning the samples into "+str(num_states)+" thermodynamic states using ")
 print("temperature windows of "+str(T_step_size)+" K.")
 state_ranges = np.array([[T_min+(i*T_step_size),T_min+((i+1)*T_step_size)] for i in range(0,num_states)])
 state_index = 0
 for state in range(0,num_states):
  T_state_center = sum(state_ranges[state_index])/2.0
  state_counts = np.array([0 for i in num_states])
  
 u_kn = np.array([[U_uncorrelated[index]/T_state_center for index in U_uncorrelated] for state in num_states])
 for T in T_uncorrelated:
  state_index = 0
  for state in state_ranges:
   if T >= state[0] and sample_T < state[1]:
    state_counts[state_index] = state_counts[state_index] + 1
   state_index = state_index + 1
 print("The distribution with "+str(num_states)+" states is:"+str(state_counts))
 distributions_for_each_num_states.append(state_counts)
 T_ranges_for_each_num_states.append(state_ranges)
# Run MBAR to calculate weights for each state
 mbar = pymbar.MBAR(U_reduced_inter_state, state_counts, verbose=False, relative_tolerance=1e-12)
 weights_for_each_num_states.append(mbar.getWeights)

#############
#
# 5) Calculate free energy surface with MBAR using Langevin dynamics simulation data
#
#
#############

# Get the dimensionless free energy differences for all states
 free_energies_for_each_num_states.append(mbar.getFreeEnergyDifferences)
 
# 
