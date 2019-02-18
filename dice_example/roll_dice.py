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
      args.apend(1)                 # Input of roll_loaded_dice(n, number, odds) --> Where roll_loaded_dice(10, 4,1) should give the distribution of a fair dice
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
   unfair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-3.3025} # -1 + ln(10) - ln(1)

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
      print('Variance fair dice = '+str(np.std(samples)))
      print('Average loaded dice = '+str(np.mean(samples_loaded)))
      print('Variance loaded dice = '+str(np.std(samples_loaded)))
      print('MBAR Average loaded dice = '+str(results.computeExpectations(samples)[0][1]))
      print('MBAR Variance loaded dice = '+str(np.sum(results.getWeights()[:,1]*np.power(samples,2)) - np.power(results.computeExpectations(samples)[0][1],2)))
      print('MBAR Partition Function :')
      print(np.exp(-1*results.getFreeEnergyDifferences()[0]))
   

   # Plotting

   if plot == True:
      plt.figure(figsize=[10,10])
      plt.hist(samples, density = True, bins=range(0,8), align='left', rwidth=0.9  )
      plt.xlabel('Dice Roll', fontsize=20)
      plt.ylabel('Density', fontsize=20)
      plt.xlim([0.5,6.5])
      # plt.hist(samples_loaded, bins=range(1,7))
      plt.gca()
      plt.savefig('outputs/1_fair_dice.png')
      plt.show()
      
      #plt.figure(figsize=[10,10])
      a = plt.hist(samples, density = True, bins=range(0,8), align='left', rwidth=0., alpha = 0.4)
      a_bins = np.array([(a[1][i] + a[1][i+1])/2 for i in range(len(a[1])-1)])
      b = plt.hist(samples, bins=range(0,8), density=True, weights=results.getWeights()[:,1], rwidth = 0.9, align='left', alpha=0.4)
      b_bins = np.array([(b[1][i] + b[1][i+1])/2 for i in range(len(b[1])-1)])
      c = plt.hist(samples_loaded, density = True, bins=range(0,8), align='left', rwidth=0.9, alpha = 0.4)
      c_bins = np.array([(c[1][i] + c[1][i+1])/2 for i in range(len(c[1])-1)])
      plt.show()


      #print(a[0])
      #print(a_bins)
      plt.figure(figsize=[10,10])
      plt.bar(a_bins-0.8, a[0], width=0.3, align='center')  # Needed to shift all bins over by -0.5, messy, but gives the plot I wanted
      plt.bar(c_bins-0.5, c[0], width=0.3, align='center')
      plt.bar(b_bins-0.2, b[0], width=0.3, align='center')
      plt.xlim([0.5,6.5])
      plt.legend(['Fair Dice', 'Weighted Dice', 'MBAR Weighted Dice' ])
      plt.xlabel('Dice Roll', fontsize=20)
      plt.ylabel('Density', fontsize=20)
      plt.gca()
      plt.savefig('outputs/1_weighted_dice.png')
      plt.show()



      plt.scatter(samples,results.getWeights()[:,1])
      plt.show()