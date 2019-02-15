# Adaption applicable for the case of NaCl system taking ion-pair distance as the reaction coordinate

from __future__ import print_function
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
# Constants.
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = 43 # number of umbrellas
N_max = 2501 # maximum number of snapshots/simulation (time: 10 ns)
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
d_min = 0.25 # min for PMF (nm)
d_max = 0.91 # max for PMF (nm)
nbins = 200 # number of bins for 1D PMF

# Allocate storage for simulation data
N_k = numpy.zeros([K], dtype = int) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K]) # K_k[k] is the spring constant (in kJ/mol/nm**2) for umbrella simulation k
d0_k = numpy.zeros([K]) # d0_k[k] is the spring center location (in nm) for umbrella simulation k
d_kn = numpy.zeros([K,N_max]) # d_kn[k,n] is the ion-pair distance (in nm) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max]) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = numpy.zeros([K])

# Read in umbrella spring constants and centers.
infile = open('data/centers.dat', 'r')
lines = infile.readlines()
infile.close()
for k in range(K):
    # Parse line k.
    line = lines[k]
    tokens = line.split()
    d0_k[k] = float(tokens[0]) # spring center locatiomn (in nm), 1st column in centers.dat
    K_k[k] = float(tokens[1])  # spring constant (in kJ/mol/nm**2), 2nd column in centers.dat
    if len(tokens) > 2:
        T_k[k] = float(tokens[2])  # temperature the kth simulation was run at.

beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
DifferentTemperatures = True
if (min(T_k) == max(T_k)):
    DifferentTemperatures = False            # if all the temperatures are the same, then we don't have to read in energies.
# Read the simulation data
for k in range(K):
    # Read ion-pair distance data.
    filename = 'data/pullx-umbrella%d.xvg' % k
    print("Reading %s..." % filename)
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    # Parse data.
    n = 0
    m = 0
    for line in lines:
        if line[0] == '#' or line[0] == '@':
            m +=1 #number of parameter lines

    for line in lines[m:m+N_max]: #read in data starting from (m+1)-th line and read in N_max lines in total
        if line[0] != '#' and line[0] != '@':
            tokens = line.split()
            d = float(tokens[1]) # ion-pair distance
            d_kn[k,n] = d
            
            n += 1
    N_k[k] = n

    if (DifferentTemperatures):  # if different temperatures are specified the metadata file, 
                                 # then we need the energies to compute the PMF
        # Read energies
        filename = 'data/umbrella%d_energies.xvg' % k
        print("Reading %s..." % filename)
        infile = open(filename, 'r')
        lines = infile.readlines()
        infile.close()
        # Parse data.
        n = 0
        m = 0
        for line in lines:
            if line[0] == '#' or line[0] == '@':
                m +=1 #number of parameter lines

        for line in lines[m:m+N_max]: #read in data starting from (m+1)-th line and read in N_max lines in total
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()            
                u_kn[k,n] = beta_k[k] * (float(tokens[2]) - float(tokens[1])) # reduced potential energy without umbrella restraint
                n += 1

    # Compute correlation times for potential energy and d (ion-pair distance)
    # timeseries.  If the temperatures differ, use energies to determine samples
            
    if (DifferentTemperatures):        
        g_k[k] = timeseries.statisticalInefficiency(u_kn[k,:], u_kn[k,0:N_k[k]])
        print("Correlation time for set %5d is %10.3f" % (k,g_k[k]))
        indices = timeseries.subsampleCorrelatedData(u_kn[k,0:N_k[k]])
    else:
        d_temp = d_kn[k,0:N_k[k]]
        g_k[k] = timeseries.statisticalInefficiency(d_temp)
        print("Correlation time for set %5d is %10.3f" % (k,g_k[k]))
        indices = timeseries.subsampleCorrelatedData(d_temp, g=g_k[k]) 
    # Subsample data.
    N_k[k] = len(indices)
    u_kn[k,0:N_k[k]] = u_kn[k,indices]
    d_kn[k,0:N_k[k]] = d_kn[k,indices]

N_max = numpy.max(N_k) # shorten the array size
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

# Construct torsion bins
print("Binning data...")
delta = (d_max - d_min) / float(nbins)
# compute bin centers
bin_center_i = numpy.zeros([nbins], numpy.float64)
for i in range(nbins):
    bin_center_i[i] = d_min + delta/2 + delta * i
# Bin data
bin_kn = numpy.zeros([K,N_max], numpy.int32)
for k in range(K):
    for n in range(N_k[k]):
        # Compute bin assignment.
        bin_kn[k,n] = int((d_kn[k,n] - d_min) / delta)

# Evaluate reduced energies in all umbrellas
print("Evaluating reduced potential energies...")
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image ion-pair distance deviation from umbrella center l
        dd = d_kn[k,n] - d0_k

        # Compute energy of snapshot n from simulation k in umbrella potential l
        u_kln[k,:,n] = u_kn[k,n] + beta_k[k] * (K_k/2.0) * dd**2

# Initialize MBAR.
print("Running MBAR...")
mbar = pymbar.MBAR(u_kln, N_k, verbose = True)

# Compute PMF in unbiased potential (in units of kT).
results = mbar.computePMF(u_kn, bin_kn, nbins)
f_i = results[0] #original: f_i=results['f_i']
df_i = results[1] #original: f_i=results['df_i']

# Write out PMF
print("PMF (in units of kT)")
print("%8s %8s %8s" % ('bin', 'f', 'df'))
for i in range(nbins):
    print("%8.3f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i]))

