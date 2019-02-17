# This Python script contains an example NaCl simulation in OpenMM, 
# followed by analysis with MDTraj, MBAR, and WHAM.

# The protocol is broken up into the following sections:

# 1) Provide user input
# 2) Prepare and run OpenMM simulation
# 3) Analyze OpenMM simulation data with MDTraj
# 4) 

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
total_simulation_time = 1 # Units = picoseconds
simulation_temperature = 300 # Units = Kelvin
kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)
T_window_range = range(2,20)
output_file_name = 'output/output.dat'
num_simulations
simulation_steps = 10000
print_frequency = 1 # Number of steps to skip when printing output
nskip = 10 # Number of steps to skip when reading timeseries data to find the equilibration time
data_exists = False

#############
#
# 2) Prepare and run OpenMM simulation
#
#############

# An OpenMM simulation requires three input objects: a 'system', an 'integrator', and a 'context'

system = SodiumChlorideCrystal() # Define a system
integrator = LangevinIntegrator(simulation_temperature, total_simulation_time, simulation_time_step) # Define an integrator
simulation = Simulation(system.topology, system.system, integrator) # Define a simulation 'context'
print(system.positions)
simulation.context.setPositions(system.positions) # Assign particle positions for this context
simulation.minimizeEnergy() # Set the simulation type to energy minimization

if data_exists == False:
# Declare data that we want to print during the simulation
 simulation.reporters.append(StateDataReporter(output_file_name, print_frequency, \
  step=True, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True))
# Run the simulation
 simulation.step(total_steps)

#############
#
# 3) Analyze OpenMM simulation data with MDTraj
#
#############

# Use simulation output to generate an ensemble of decorrelated samples
output_obj = open(output_file_name,'r')
# E_total_all stores the total energies from NaCl simulation output
E_total_all_temp = np.array([l.split(',')[3] for l in output_obj.readlines()])
# T_all stores the Temperatures from NaCl simulation output
T_all_temp = np.array([l.split(',')[4] for l in output_obj.readlines()])
output_obj.close()
E_total_all = np.array(np.delete(E_total_all_temp,0,0),dtype=float)
T_all = np.array(np.delete(T_all_temp,0,0),dtype=float)
# Use pymbar's timeseries() function to analyze U_all, identify the equilibration time (t0), and evaluate the statistical inefficiency (g) for all points in the timeseries
# THIS STEP IS TIME-CONSUMING.  TO SPEED UP INCREASE 'nskip'
[t0, g, Neff_max] = timeseries.detectEquilibration(E_total_all,nskip=nskip)
print("The equilibration time was "+str(t0*simulation_time_step))
E_total_equil = E_total_all[t0:]
T_equil = T_all[t0:]
uncorrelated_energies = timeseries.subsampleCorrelatedData(E_total_equil, g=g) # Determine indices of uncorrelated samples
print("Sub-sampled simulation data with "+str(len(indices))+" uncorrelated samples")
U_uncorrelated = U_equil[indices] # Uncorrelated total energies 
np.savetxt('output/uncorrelated_U.dat',uncorrelated_data)
T_uncorrelated = T_equil[indices] # Uncorrelated temperatures
# Calculate the reduced potential energies for the uncorrelated samples
U_reduced = np.array([U_uncorrelated[index]/(T_uncorrelated[index]*kB)])
E_
# Now we are ready to use U_reduced to optimize a set of weights, whose values indicate the probability of occupying a specific thermodynamic state
# We split the full range of temperatures sampled during the simulation into 'bins', 
# whose size is calculated as: ( T_max - T_min ) / nbins
T_max = max(T_uncorrelated)
T_min = min(T_uncorrelated)
window_ranges_for_each_num_T_windows = np.array([])
distributions_for_each_num_T_windows = np.array([])
for num_T_windows in T_window_range:
# First we bin (count the total number) of samples in each thermodynamic state
 T_step_size = (T_max - T_min) / num_T_windows
 print("Binning the samples using "+str(num_T_windows)+" temperature windows of "+str(T_step_size)+" K.")
 window_ranges = np.array([[T_min+(i*T_step_size),T_min+((i+1)*T_step_size)] for i in range(0,num_T_windows)])
 T_window_center = sum(window_ranges[window_index])/2.0
 window_counts = np.array([0 for i in num_T_windows])
 U_reduced_inter_window = np.array([[U_uncorrelated[index]/T_window_center for index in U_uncorrelated] for window in num_T_windows]) # Equivalent to u_kn in mbar.py
 for sample_T in T_uncorrelated:
  window_index = 0
  for window in window_ranges:
   if sample_T >= window[0] and sample_T < window[1]:
    window_counts[window_index] = window_counts[window_index] + 1
   window_index = window_index + 1
 print("The window distribution with "+str(num_T_windows)+" windows is:"+str(window_counts))
# Run MBAR to calculate weights for each T_window
 mbar = pymbar.MBAR(U_reduced_inter_window, window_counts, verbose=False, relative_tolerance=1e-12)
 distributions_for_each_num_T_windows.append(window_counts)
 window_ranges_for_each_num_T_windows.append(window_ranges)
np.savetxt('output/distributions_by_num_windows.dat',distributions_for_each_num_T_windows)
np.savetxt('output/window_ranges_by_num_windows.dat',window_ranges_for_each_num_T_windows) 
 
