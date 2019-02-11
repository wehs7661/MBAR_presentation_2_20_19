#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pymbar
import scipy.stats as sps

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
# if True:
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

   # print(u_unfair)
   #print(samples)
   #exit()
   u_kln = np.array([u_fair, u_unfair])

   print(u_kln)
   
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
   print(N_samples)
   sums = list(map(sum, samples))
   sums_loaded = list(map(sum, samples_loaded))

   if True:
      N_k = np.array([N_samples, 0])
      results = pymbar.mbar.MBAR(u_kln, N_k)
      print(results.getFreeEnergyDifferences()[0])
      print(results.getWeights())
      print('Average fair dice = '+str(np.mean(sums)))
      print('Average loaded dice = '+str(np.mean(sums_loaded)))
      print('MBAR Average loaded dice = '+str(results.computeExpectations(sums)[0][1]))
   

   # Plotting

   if plot == True:
      # Histograms
      plt.hist(list(map(sum, samples)), alpha=0.5, bins=range(1,6*n_dice))
      plt.hist(list(map(sum, samples_loaded)), alpha=0.5, bins=range(1,6*n_dice))
      plt.show()

      # KDE of data

      density = sps.gaussian_kde(list(map(sum, samples)))
      density_loaded = sps.gaussian_kde(list(map(sum, samples_loaded)))

      x_range = np.linspace(1*n_dice, 6*n_dice, 200)
      plt.plot(x_range, density(x_range))
      plt.plot(x_range, density_loaded(x_range))
      plt.show()


      # Scatter of weights

      plt.scatter(sums, results.getWeights()[:,1])
      plt.show()
