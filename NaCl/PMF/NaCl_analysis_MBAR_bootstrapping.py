# Adaption applicable for the case of NaCl system taking ion-pair distance as the reaction coordinate

from __future__ import print_function
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
import random
import scipy.stats
import matplotlib.pyplot as plt
# Constants.
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = 47 # number of umbrellas
N_max = 2501 # maximum number of snapshots/simulation (time: 10 ns)
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol)) 
d_min = 0.25 # min for PMF (nm)
d_max = 1.00 # max for PMF (nm)
nbins = 200 # number of bins for 1D PMF
RunBootstrapping = True
nboot = 200 # number of bootstrap samples

# Allocate storage for simulation data
N_k = numpy.zeros([K], dtype = int) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K]) # K_k[k] is the spring constant (in kJ/mol/nm**2) for umbrella simulation k
d0_k = numpy.zeros([K]) # d0_k[k] is the spring center location (in nm) for umbrella simulation k
d_kn = numpy.zeros([K,N_max]) # d_kn[k,n] is the ion-pair distance (in nm) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max]) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = numpy.zeros([K]) #scaling factor for the k-th simulation (statistical inefficiency, g=1+2tau)

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
            
    if (DifferentTemperatures):     # for different temperatures   
        g_k[k] = timeseries.statisticalInefficiency(u_kn[k,:], u_kn[k,0:N_k[k]])
        print("Correlation time for set %5d is %10.3f" % (k,g_k[k]))
        indices = timeseries.subsampleCorrelatedData(u_kn[k,0:N_k[k]])
    else:                          # for constant temperature
        d_temp = d_kn[k,0:N_k[k]]    # all the distance data (from first one the N_max-th one)
        g_k[k] = timeseries.statisticalInefficiency(d_temp)     
        print("Correlation time for set %5d is %10.3f" % (k,g_k[k]))   # actually, g_k[k] is not correlation time...
        indices = timeseries.subsampleCorrelatedData(d_temp, g=g_k[k]) # list of int
    # Note: g is the statistical inefficiency, g = 1 + 2 tau, where tau is the correlation time
    # statisticalInefficiency(A_n, B_n=None, fast=False, mintime=3, fft=False)
    # timeseries.subsampleCorrelatedData determines the indices of an uncorrelated subsample of the data.
    
    # Subsample data.
    N_k[k] = len(indices)    # At this point, N_k contains the number of uncorrelated samples for each state k                
    u_kn[k,0:N_k[k]] = u_kn[k,indices]    # now the array only contains data of uncorrelated samples
    d_kn[k,0:N_k[k]] = d_kn[k,indices]    # now the array only contains data of uncorrelated samples
    # note there is no zero entries in d_kn, the elements with indices bigger than the number of subsamples
    # are just the original values assign to d_kn before we subsample the data.

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
        bin_kn[k,n] = int((d_kn[k,n] - d_min) / delta) # contains the information about which bin was the data assigned to
        
# Evaluate reduced energies in all umbrellas
print("Evaluating reduced potential energies...")
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image ion-pair distance deviation from umbrella center l
        dd = d_kn[k,n] - d0_k          # delta d

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

plt.plot(bin_center_i, f_i)
plt.xlabel('Ion-pair distance (nm)')
plt.ylabel('Potential of mean force (PMF) (kT)')
    
    
"""The following section is used for the bootstrapping of MBAR"""   
if RunBootstrapping == True:
    print('Start bootstrapping!')
    print('Number of bootstrap sample:', nboot)
    u_kln_b = numpy.zeros([K, K, N_max, nboot], numpy.float64)  
    u_kn_b = numpy.zeros([K,N_max, nboot])
    d_kn_b = numpy.zeros([K,N_max, nboot])


    index = 0
    for j in range(nboot):
        for i in range(K):
            # first randomly select N_k[0] (number of uncorrelated samples) samples from the data of i-th simulation
            d_kn_rand = random.choices(d_kn[i,:], k=N_k[i])    # length = number of uncorrelate samples
            u_kn_rand = random.choices(u_kn[i,:], k=N_k[i])    # length = number of uncorrelate samples
            # then assign values to bootstrapping array
            d_kn_b[i,:N_k[i],j] = d_kn_rand
            u_kn_b[i,:N_k[i],j] = u_kn_rand

    # Set zero of u_kn -- this is arbitrary.
    u_kn_b -= u_kn_b.min() # might be a problem if not constant temp (to be modified)

    bin_kn_b = numpy.zeros([K,N_max,nboot], numpy.int32)
    for j in range(nboot):
        for k in range(K):
            for n in range(N_k[k]):
                # Compute bin assignment.
                bin_kn_b[k,n,j] = int((d_kn_b[k,n,j] - d_min) / delta)  # contains the information about which bin was the data assigned to
                
                # Compute minimum-image ion-pair distance deviation from umbrella center l
                dd_b = d_kn_b[k,n,j] - d0_k      # delta d (bootstrapping)

                # Compute energy of snapshot n from simulation k in umbrella potential l
                u_kln_b[k,:,n,j] = u_kn_b[k,n,j] + beta_k[k] * (K_k/2.0) * dd_b**2   

    MBAR_b = []          # A list of object of pymbar.MBAR(u_kln_b[j], N_k, verbose = True)
    bootdata = numpy.zeros([nboot, nbins]) 
    # bootdata is an array storing bootstrapping data for each value of reaction coordinate (data of f_i)

    for j in range(nboot):
        # Initialize MBAR.
        print("Running MBAR...")
        print("Calculating bootstrap sample", j+1)
        MBAR_b.append(pymbar.MBAR(u_kln_b[:,:,:,j], N_k, verbose = True, initial_f_k = mbar.f_k))
        # initial_f_k = m_bar.f_k is use the mbar.f_k calculated previously as the initial guess, which can speed up the process (default = None, that is f_k = [0, 0, ...,0])
        # take a lookt a __init__ of pymbar and see how to output f_k data --> just use mbar.f_k !
        
        # Compute PMF in unbiased potential (in units of kT).
        results_b = MBAR_b[j].computePMF(u_kn_b[:,:,j], bin_kn_b[:,:,j], nbins)
        bootdata[j,:] = results_b[0]      # results_b[0]: data of f_i
        #f_i_b = results_b[0]
        #df_i_b = results_b[1]   

    # Now let's calculate the mean, variance and the confidence interval of bootstrap samples!
    m = numpy.zeros(nbins)       # mean 
    var = numpy.zeros(nbins)     # variance (sample variance)
    sem = numpy.zeros(nbins)     # standard error of means
    h = numpy.zeros(nbins)       # for the range of confidence interval
    confidence = 0.95            # confidence level
    # note: as a rule of thumb, it requires about 50-200 bootstrap samples to get a good variance and 1000-5000 bootstrap samples to get a good 95CI 
    
    for i in range(nbins):
        m[i] = numpy.mean(bootdata[:,i])
        var[i] = numpy.var(bootdata[:,i])
        sem[i] = scipy.stats.sem(bootdata[:,i])
        # note that the default of ddof of numpy.var is 0, we have to make it as one to calculate sample variance instead of the population variance
        h[i] = sem[i] * scipy.stats.t.ppf((1 + confidence) / 2., len(bootdata[:,i])-1)
        # then the confidence interval is from m-h to m+h
        
    # error bar based on variance or CI
    plt.figure()
    plt.title('Errorbar: variance')
    plt.errorbar(bin_center_i, f_i, yerr = var)
    plt.xlabel('Ion-pair distance (nm)')
    plt.ylabel('Potential of mean force (PMF) (kT)')
    
    plt.figure()
    plt.title('95% of confidence level')
    plt.plot(bin_center_i, f_i, color = 'blue')
    plt.fill_between(bin_center_i, m-h, m+h, color='yellow')
    plt.xlabel('Ion-pair distance (nm)')
    plt.ylabel('Potential of mean force (PMF) (kT)')

