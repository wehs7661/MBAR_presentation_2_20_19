#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pymbar

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


if __name__ == "__main__":
# if True
   plot = True

   # Samples for a fair dice

   N_samples = 10000
   k_sims = 2
   n_dice = 1
   bins = dict()
   samples =list()
   for i in range(N_samples): 
      roll = int(roll_dice(n_dice))
      samples.append(roll)   # Writing samples for fair dice
      if roll not in bins.keys():
         bins[roll] = 1
      else:
         bins[roll] += 1
   freq = list(map(lambda x: bins[x], bins.keys()))


   fair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1}
   unfair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-3.3025}

   u_fair = [fair[x] for x in samples]
   u_unfair = [unfair[x] for x in samples]

   u_kln = np.array([u_fair, u_unfair])

   print(u_kln)
   
   # Samples for a loaded dice 
   if True:
      u_loaded = list()
      N_samples = 10000
      n_dice = 1
      bins_loaded = dict()
      samples_loaded =list()

      for i in range(N_samples):
         roll = int(roll_loaded_dice(n_dice,6,10))
         samples_loaded.append(roll)
         if roll not in bins_loaded.keys():
            bins_loaded[roll] = 1
         else:
            bins_loaded[roll] += 1
      freq_loaded = list(map(lambda x: bins_loaded[x], bins_loaded.keys()))

   # MBAR stuff

   # u_kn = np.array([
   # state 1 [16, 18, 10, 34, 23, ...]
   # state 2 [18, 26, ...]
   # ])
   print(N_samples)
   if True:
      N_k = np.array([N_samples, 0])
      results = pymbar.mbar.MBAR(u_kln, N_k)
      print(results.getFreeEnergyDifferences()[0])
      print(results.getWeights())
      print('Average fair dice = '+str(np.mean(samples)))
      print('Average loaded dice = '+str(np.mean(samples_loaded)))
      print('MBAR Average loaded dice = '+str(results.computeExpectations(samples)[0][1]))
   

   # Plotting

   if plot == True:
      plt.hist(samples)
      plt.hist(samples_loaded)
      plt.show()
