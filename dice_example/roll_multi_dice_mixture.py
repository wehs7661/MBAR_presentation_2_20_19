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
    unfair_samples = [roll_loaded_dice(n_dice, 6, 10) for _ in range(N_samples)]
    unfair_test = [roll_loaded_dice(n_dice, 6, 3) for _ in range(N_samples)]
    unfair_test_sums = [ sum(x) for x in unfair_test]
    fair = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1}
    unfair = {6:-3.3026, 2:-1, 3:-1, 4:-1, 5: -1, 1:-1}
    unfair_test = {6:-2.0986, 2:-1, 3:-1, 4:-1, 5: -1, 1:-1}
    #low_roller = {1:- 3.3026, 2:-3.3026, 3:-3.3026, 4: -1, 5: -1, 6: -1}
    dists = [fair_samples, unfair_samples]
    N_k = [N_samples, N_samples, 0]
    energys = [fair, unfair, unfair_test]
    woot = mbarDiceMixtureDist(dists, N_k, energys)
    print(woot.mbar.computeExpectations(woot.mix_dist_sum))
    
    plt.figure(figsize=[10,10])
    plt.hist(woot.mix_dist_sum, alpha = 0.3, density=True, bins=range(1,61))
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,0], bins=range(1,61))
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,1], bins=range(1,61))
    # plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,2], bins=range(1,61))
    # plt.hist(unfair_test_sums, alpha = 0.3, bins=range(1,61),density=True)
    plt.legend(['Mixture Distribution','Fair Dice', 'Loaded Dice'])
    plt.xlabel('Dice Value', fontsize=20)
    plt.ylabel('Density' ,fontsize=20)
    plt.xlim([n_dice, n_dice*6])
    plt.gca()
    plt.savefig('outputs/mixture_dist.png')
    plt.show()

    plt.figure(figsize=[10,10])
    plt.hist(woot.mix_dist_sum, alpha = 0.3, density=True, bins=range(1,61))
    # plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,0], bins=range(1,61))
    # plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,1], bins=range(1,61))
    plt.hist(unfair_test_sums, alpha = 0.3, bins=range(1,61),density=True)
    plt.hist(woot.mix_dist_sum, alpha =0.3, weights=woot.mbar.getWeights()[:,2], bins=range(1,61))
    plt.legend(['Mixture Distribution', 'Loaded Dice','MBAR Weighted Dice'])
    plt.xlabel('Dice Value', fontsize=20)
    plt.ylabel('Density' ,fontsize=20)
    plt.xlim([n_dice, n_dice*6])
    plt.gca()
    plt.savefig('outputs/mixture_dist_estimates.png')
    plt.show()

        
        