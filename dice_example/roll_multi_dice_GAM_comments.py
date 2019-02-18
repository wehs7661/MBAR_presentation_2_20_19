#!python

#GAM = Comments added by Garrett A. Meek on 2/18/2019

import numpy as np
import matplotlib.pyplot as plt
import pymbar
import scipy.stats as sps

#GAM Adding the random module
import random

#Goal: Use MBAR to estimate "weights" needed to bias a dice towards a certain expectation value. Find FE between a fair dice and a loaded dice.

# Defining Functions
def roll_dice(n):
   return(np.random.randint(1,7,n))

def roll_loaded_dice(n,*args):
   if len(args) % 2 != 0:          # Function can take in multiple changes to the probability of dice rolls
      arg.apend(1)                 # Input of roll_loaded_dice(n, number, odds) --> Where roll_loaded_dice(10, 4,1) should give the distribution of a fair dice
   choices = list(range(1,7))
   for i in range(int(len(args)/2)):
      choices.extend([args[2*i]]*(args[2*i+1]-1))           # Interesting (bad) logic to append number (odds - 1) times to the fair dice distribution
   return(np.random.choice(choices, n))

#Main


plot = True

   # Samples for a fair dice

N_samples = 100000
k_sims = 2
n_dice = 10
bins = dict()
samples =list()
for i in range(N_samples): 
 roll = roll_dice(n_dice)
 samples.append(roll)   # Writing samples for fair dice
      #if roll not in bins.keys():
      #   bins[roll] = 1
      #else:
      #   bins[roll] += 1
   #freq = list(map(lambda x: bins[x], bins.keys()))


fair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1}
unfair = {1:-2.0986, 2:-1, 3:-1, 4:-1, 5: -1, 6:-1}
u_fair = [sum(list(map(lambda y: fair[y], x))) for x in samples]
u_unfair = [sum(list(map(lambda y: unfair[y], x))) for x in samples]
#print(type(u_fair))
#   print(samples[0:100])

#GAM Suppose we define two additional u_kln where, instead of assigning 'states'
#GAM whose definitions are identical to the distribution a sample was drawn from,
#GAM we define states which span different total values for the 10 dice that are rolled.
#GAM In other words, for both u_fair and u_unfair we define two states:

#GAM State 1: Totals from -10 to -12
#GAM State 2: Totals < -12

#GAM So now we have a total of 4 'states' in u_kln, but we are still using the same total number of samples,
#GAM and, importantly, in MBAR, we don't care which distributions these samples come from.  (Stated differently,
#GAM in the language of MD simulations, we don't care which simulation a sample came from, since we throw them
#GAM all into the mixture distribution anyway and reweight them.

#GAM To illustrate these properties we can perform three tests:

#GAM 1) Keep the samples for the "fair" and "unfair" distributions in separate 'states'

#GAM 2) Randomly swap half of the samples between corresponding states in the "fair" and "unfair" distributions

#GAM 3) Randomly swap half of the samples between non-corresponding states in the "fair" and unfair" distributions

#GAM First we will generate "pure" states for the fair and unfair distributions:
u_fair_state_1 = []
u_fair_state_2 = []
u_unfair_state_1 = []
u_unfair_state_2 = []
for i in u_fair:
 if i >= -12.0:
  u_fair_state_1.extend([i])
 if i < -12. :
  u_fair_state_2.extend([i])
for i in u_unfair:
 if i >= -12.0:
  u_unfair_state_1.extend([i])
 if i < -12.0:
  u_unfair_state_2.extend([i])
#GAM Now we are ready to construct the u_kln for test 1:
#GAM Now we fill in zeros for all missing states:
#GAM 
max_num_samples_per_state = max([len(u_fair_state_1),len(u_fair_state_2),len(u_unfair_state_1),len(u_unfair_state_2)])
if len(u_fair_state_1) < max_num_samples_per_state: 
 u_fair_state_1.extend([0. for i in range(0,max_num_samples_per_state-len(u_fair_state_1))])
if len(u_fair_state_2) < max_num_samples_per_state: 
 u_fair_state_2.extend([0. for i in range(0,max_num_samples_per_state-len(u_fair_state_2))])
if len(u_unfair_state_1) < max_num_samples_per_state: 
 u_unfair_state_1.extend([0. for i in range(0,max_num_samples_per_state-len(u_unfair_state_1))])
if len(u_unfair_state_2) < max_num_samples_per_state:
 u_unfair_state_2.extend([0. for i in range(0,max_num_samples_per_state-len(u_unfair_state_2))])

u_kln_test_1 = np.array([u_fair_state_1,u_fair_state_2,u_unfair_state_1,u_unfair_state_2])

#print(u_kln_test_1[0])

if True:
#GAM Make a copy of the distributions for the pure states
   pure_states = np.array([u_fair_state_1,u_fair_state_2,u_unfair_state_1,u_unfair_state_2])
#GAM Next we randomly swap half of the samples between corresponding states in each distribution:
   for i in range(0,max_num_samples_per_state/2):
    swap = random.randint(0,max_num_samples_per_state)
#    print(swap)
#    print(len(u_fair_state_1))
    temp = u_fair_state_1[swap]
    u_fair_state_1[swap] = u_unfair_state_1[swap]
    u_unfair_state_1[swap] = temp
    swap = random.randint(0,max_num_samples_per_state)
#    print(swap)
#    print(len(u_fair_state_2))
    temp = u_fair_state_2[swap]
    u_fair_state_2[swap] = u_unfair_state_2[swap]
    u_unfair_state_2[swap] = temp

#GAM Now we are ready to construct u_kln for test 2:
   u_kln_test_2 = np.array([u_fair_state_1,u_fair_state_2,u_unfair_state_1,u_unfair_state_2])

#GAM Next we use the pure states that we stored earlier, and randomly swap samples between
# NON-corresponding states from the fair and unfair distributions:
   for i in range(0,max_num_samples_per_state/2):
    swap = random.randint(0,max_num_samples_per_state)
    temp = pure_states[0][swap]
    pure_states[0][swap] = pure_states[3][swap]
    pure_states[3][swap] = temp
    swap = random.randint(0,max_num_samples_per_state)
    temp = pure_states[1][swap]
    pure_states[1][swap] = pure_states[2][swap]
    pure_states[2][swap] = temp
#GAM Now we are ready to construct u_kln for test 3:
   u_kln_test_3 = np.array([u_fair_state_1,u_fair_state_2,u_unfair_state_1,u_unfair_state_2])

   u_kln = np.array([u_fair, u_unfair])
   
   # Samples for a loaded dice 
   if True:
      u_loaded = list()
      # n_dice = 2
      bins_loaded = dict()
      samples_loaded =list()

      for i in range(N_samples):
         roll = roll_loaded_dice(n_dice,1,3)
         samples_loaded.append(roll)
         #if roll not in bins_loaded.keys():
         #   bins_loaded[roll] = 1
         #else:
         #   bins_loaded[roll] += 1
      #freq_loaded = list(map(lambda x: bins_loaded[x], bins_loaded.keys()))

   # MBAR stuff

   # u_kn = np.array([
   # state 1 [16, 18, 10, 34, 23, ...]
   # state 2 [18, 26, ...]
   # ])
   sums = list(map(sum, samples))
   sums_loaded = list(map(sum, samples_loaded))

   if True:
#GAM I'm wondering why you define N_k this way?  I think this means that all of your uncorrelated samples
#GAM will come from the first distribution, exclusively...
      N_k = np.array([N_samples, 0])
      results_1 = pymbar.mbar.MBAR(u_kln, N_k)
      N_k = np.array([N_samples/2, N_samples/2])
      results_2 = pymbar.mbar.MBAR(u_kln, N_k)
#GAM Now we're going to generate MBAR results for our three test cases:
      N_k = np.array([N_samples/4,N_samples/4,N_samples/4,N_samples/4])
      results_test_1 = pymbar.mbar.MBAR(u_kln_test_1, N_k)
      results_test_2 = pymbar.mbar.MBAR(u_kln_test_2, N_k)
      results_test_3 = pymbar.mbar.MBAR(u_kln_test_3, N_k)
      print('Free Energy Differences with two states (unevenly distributed): ')
      print(results_1.getFreeEnergyDifferences()[0])
#      results_test_3 = pymbar.mbar.MBAR(u_kln_test_3, N_k)
      print('Free Energy Differences with two states (evenly distributed): ')
      print(results_2.getFreeEnergyDifferences()[0])
#GAM Printing out the free energy differences for these three tests also, for comparison:
      print('Free Energy Differences with four states (Test #1): ')
      print(results_test_1.getFreeEnergyDifferences()[0])
      print('Free Energy Differences with four states and random, correspondent state swaps (Test #2): ')
      print(results_test_2.getFreeEnergyDifferences()[0])
      print('Free Energy Differences with four states and random, non-correspondent state swaps (Test #3): ')
      print(results_test_3.getFreeEnergyDifferences()[0])
      print('Weights with two states: ')
      print(results_1.getWeights())
      print('Weights with four states (Test #1): ')
      print(results_test_1.getWeights())
      print('Weights with four states and random, corespondent state swaps (Test #2): ')
      print(results_test_2.getWeights())
      print('Weights with four states and random, non-correspondent state swaps (Test #3): ')
      print(results_test_3.getWeights())
      print('Average fair dice = '+str(np.mean(sums)))
      print('Average loaded dice = '+str(np.mean(sums_loaded)))
      print('MBAR Average loaded dice = '+str(results_1.computeExpectations(sums)[0][1]))
   

   # Plotting

   if plot == True:
      # Histograms
      plt.hist(list(map(sum, samples)), alpha=0.5, bins=range(1,6*n_dice+1))
      plt.xlabel('Roll Energy')
      plt.ylabel('Frequency')
      plt.legend(['Fair Dice'])
      # plt.hist(list(map(sum, samples_loaded)), alpha=0.5, bins=range(1,6*n_dice))
      plt.show()

      # KDE of data

      density = sps.gaussian_kde(list(map(sum, samples)))
      density_MBAR = sps.gaussian_kde(list(map(sum, samples)), weights=results_1.getWeights()[ :, 1])
      density_loaded = sps.gaussian_kde(list(map(sum, samples_loaded)))


      x_range = np.linspace(1*n_dice, 6*n_dice, 200)
      plt.plot(x_range, density(x_range))
      plt.plot(x_range, density_loaded(x_range))
      plt.plot(x_range, density_MBAR(x_range))
      plt.show()


      # Scatter of weights

      plt.scatter(sums, results_1.getWeights()[:,1])
      plt.show()
