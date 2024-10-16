from scipy.stats import binomtest
import numpy as np
import matplotlib.pyplot as plt

n_stars = np.array([25, 50, 100, 200, 400])
alphas = [0.2, 0.4, 0.6, 0.8, 1]
fiducial_rate = 0.25


fig, ax = plt.subplots(1, 1, sharex=True)
plt.subplots_adjust(hspace=0)

for i, n_star in enumerate(n_stars):
    pvalues = []
    pvalues_low = []
    pvalues_high = []
    true_rates = np.arange(0, 2.1, 0.0001)
    int_list = [int(x) for x in n_star * fiducial_rate*true_rates]
    n_successes = np.unique(int_list)
    for n_obs in n_successes:
        result = binomtest(n_obs, n_star, p=fiducial_rate, alternative='two-sided')
        pvalues.append(result.pvalue)
        print(n_star, n_obs, n_obs/n_star, result.pvalue, result.proportion_ci(confidence_level=0.95))

        result = binomtest(int(n_obs-np.sqrt(n_obs)), n_star, p=fiducial_rate, alternative='two-sided')
        pvalues_low.append(result.pvalue)
        print(n_star, n_obs, int(n_obs-np.sqrt(n_obs)), int(n_obs-np.sqrt(n_obs))/n_star, result.pvalue, result.proportion_ci(confidence_level=0.95))

        result = binomtest(int(n_obs+np.sqrt(n_obs)), n_star, p=fiducial_rate, alternative='two-sided')
        pvalues_high.append(result.pvalue)
        print(n_star, n_obs, int(n_obs+np.sqrt(n_obs)), int(n_obs+np.sqrt(n_obs))/n_star, result.pvalue, result.proportion_ci(confidence_level=0.95))

    plt.fill_between(n_successes/n_star/fiducial_rate, pvalues_low, pvalues_high, alpha=alphas[i]-0.1, color='rebeccapurple', label='{} stars'.format(n_star))
    #plt.plot(n_successes/n_star/fiducial_rate, pvalues, alpha=alphas[i], c='rebeccapurple', label='{} stars'.format(n_star))

plt.axhline(0.05, linestyle='--')
plt.yscale('log')
plt.ylim(1e-4, 1)
plt.xlim(0,2)
plt.text(0.05, 0.1, 'Not significant', ha='left', va='center', fontsize=10, color='k')
plt.text(0.05, 0.025, 'Significant', ha='left', va='center', fontsize=10, color='k')
plt.xlabel('suppression/enhancement factor')
plt.ylabel(r'p-value compared to $\eta=0.25$')
plt.legend(loc='lower right')
plt.savefig('binomial_sim.png', dpi=250)