from roll_multi_dice import *

class mbarDiceMixtureDist:
    def __init__(self, dists, N_k, energys):
        """
        mbarDiceMixtureDist will be used to 
        make multiple examples with multiple dice
        reweighting from a mixture distribution rather
        than a single dice distribution

        dists : List of each individual distribution {1:-1, 2,:-1, 3:-1 ... }
        energys : List of energy functions corresponding to given distributions + new
        energies wished to be calculated

        """
        self.mix_dists = np.array(dists).reshape(-1, np.array(dists).shape[-1])
        print(self.mix_dists)
        self.mix_dist_sum = [ sum(x) for x in self.mix_dists]
        self.energys = energys
        self.u_kln = []
        self.N_k = N_k
        self.u_kln = [[ sum(list(map(lambda y: f_energy[y], x))) for x in self.mix_dists] for f_energy in energys ]
        self.u_kln = np.array(self.u_kln)
        self.mbar = pymbar.mbar.MBAR(self.u_kln, self.N_k)
    
    def printResults(index):
        print('Free Energy Differences : ')
        print(results.getFreeEnergyDifferences()[0])
        print('Weights : ')
        print(results.getWeights()[:,index])
        print('Average fair dice = '+str(np.mean(sums)))
        print('Average loaded dice = '+str(np.mean(sums_loaded)))
        print('MBAR Average loaded dice = '+str(results.computeExpectations(sums)[0][1]))


        

        
if __name__ == "__main__":

    N_samples = 100000
    k_sims = 2
    n_dice = 10

    fair_samples = [roll_dice(n_dice) for _ in range(N_samples)]
    unfair_samples_high = [roll_loaded_dice(n_dice, 6, 10) for _ in range(N_samples)]
    unfair_samples_low = [roll_loaded_dice(n_dice, 1, 10) for _ in range(N_samples)]
    unfair_test = [roll_loaded_dice(n_dice, 6, 3) for _ in range(N_samples)]
    unfair_test_sums = [ sum(x) for x in unfair_test]
    fair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1}
    unfair_low = {1:-3.3026, 2:-1, 3:-1, 4:-1, 5: -1, 6:-1}
    unfair_high = {6:-3.3026, 2:-1, 3:-1, 4:-1, 5: -1, 1:-1}
    unfair_test = {6:-2.0986, 2:-1, 3:-1, 4:-1, 5: -1, 1:-1}
    low_roller = {1:- 3.3026, 2:-2.0986, 3:-1, 4: -1, 5: -2.0986, 6: -3.2026}
    dists = [fair_samples, unfair_samples_high, unfair_samples_low]
    N_k = [N_samples, N_samples, N_samples, 0, 0]
    energys = [fair, unfair_low, unfair_high, unfair_test, low_roller]
    woot = mbarDiceMixtureDist(dists, N_k, energys)
    print(woot.mbar.computeExpectations(woot.mix_dist_sum))
    

    # Generate Plots for mixture distribution made from a combination of 3
    # different states

    plt.figure(figsize=[10,10])
    plt.hist(woot.mix_dist_sum, alpha = 0.3, density=True, bins=range(1,61), color='k')
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,0], bins=range(1,61))
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,1], bins=range(1,61))
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,2], bins=range(1,61))
    # plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,2], bins=range(1,61))
    # plt.hist(unfair_test_sums, alpha = 0.3, bins=range(1,61),density=True)
    plt.legend(['Mixture Distribution','Fair Dice', 'High Roll Dice', 'Low Roll Dice'])
    plt.xlabel('Dice Value', fontsize=20)
    plt.ylabel('Density' ,fontsize=20)
    plt.xlim([n_dice, n_dice*6])
    plt.gca()
    plt.savefig('outputs/mixture_dist.png')
    plt.show()



    # Defining a KDE to show the mixture distribution
    # so as to not obscure the overlap of the estimated distributions

    bw = 0.1
    x_range = np.linspace(1*n_dice, 6*n_dice, 200)

    density_mix = sps.gaussian_kde(woot.mix_dist_sum, bw_method=bw)

    # These shows validation of intermediately loaded dice
    # See how edges of samples are much more in agreement between
    # the sampled and estimated distribution

    plt.figure(figsize=[10,10])
    plt.plot(x_range, density_mix(x_range),'k')
    plt.hist(unfair_test_sums, alpha = 0.3, bins=range(1,61),density=True)
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,3], bins=range(1,61))
    plt.legend(['Mixture Distribution', 'Loaded Dice','MBAR Weighted Dice'])
    plt.xlabel('Dice Value', fontsize=20)
    plt.ylabel('Density' ,fontsize=20)
    plt.xlim([n_dice, n_dice*6])
    plt.gca()
    plt.savefig('outputs/mixture_dist_estimates.png')
    plt.show()


    # Showing an estimated distribution for strange energy functions
    # inverse quadratic?

    plt.figure(figsize=[10,10])
    plt.plot(x_range, density_mix(x_range),'k')
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,4], bins=range(1,61))
    plt.legend(['Mixture Distribution', 'Loaded Dice','MBAR Weighted Dice'])
    plt.xlabel('Dice Value', fontsize=20)
    plt.ylabel('Density' ,fontsize=20)
    plt.xlim([n_dice, n_dice*6])
    plt.gca()
    plt.savefig('outputs/mixture_dist_prediction.png')
    plt.show()
        
    


    # Scatter of weights
      
    fig, ax1 = plt.subplots(figsize=(10,15))

    color = 'black'
    ax1.plot(x_range, density_mix(x_range), color=color)
    ax1.set_xlabel('Dice Value', fontsize=20)
    ax1.set_ylabel('Mixture Distribution', fontsize=20, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([0, max(density_mix(x_range))+np.std(density_mix(x_range))])
    ax2 = ax1.twinx()

    color = 'C0'
    ax2.scatter(woot.mix_dist_sum, woot.mbar.getWeights()[:,3])
    ax2.set_ylabel('Weight Values', fontsize=20, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim([min(woot.mbar.getWeights()[:,1]), max(woot.mbar.getWeights()[:,1]) + 10*np.std(woot.mbar.getWeights()[:,1])])
    
    ax1.autoscale(enable = True, axis='x', tight=True)

    plt.gca()
    plt.savefig('outputs/mixture_dist_weights.png')
    plt.show()

    # Plot low rollwer energy function as a series of step functions

    x_die = np.linspace(1,6.99999,500)
    x_eval = np.floor(x_die)
    energy = np.array([low_roller[x] for x in x_eval])
    plt.plot(x_die, energy)
    plt.xlabel('Die Value')
    plt.ylabel('Energy')
    plt.gca()
    plt.savefig('outputs/funky_energy.png')
    plt.show()