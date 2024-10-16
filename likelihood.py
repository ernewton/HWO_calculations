import numpy as np
import matplotlib.pyplot as plt
#from scipy.stats import norm
#import emcee
import os
import sys
from scipy.stats import truncnorm
random_seed = int(sys.argv[3])

np.random.seed(random_seed)

## variables: sigma_true, n_stars, age_unc
## make one plot for sigma_true = 1, one plot for sigma_true = 0.1

"""
assume true values for the inner and outer semi-major axis where binary suppression matters
"""

influencer = sys.argv[1]

# set up parameters for the binary star distribution of a
logmu = np.log10(40) # AU
logsigma = 1.5 # log10(AU)
trunc_low_loga = 1 # 10 AU, assume nothing closer is observed
trunc_high_loga = 3 # 10^4 AU is max a shown in Offner plot
trunc_high = (trunc_high_loga-logmu)/logsigma # truncates at no. of sigma from loc
trunc_low = (trunc_low_loga-logmu)/logsigma # truncates at no. of sigma from loc

single_detection_prob = 0.25 # prob of whatever planet I am simulating existing around a single star

"""
draw a dataset
"""

n_stars = int(sys.argv[2])
MF = 1 #0.47 # multiplicity fraction
n_binaries = int(n_stars*MF)
#n_singles = n_stars - n_binaries


if influencer == 'suppression':
    a_inner_true = 10  # AU (suppression 100%)
    a_outer_true = 200  # AU (suppression 0%)

    working_dir = "results/suppression_p{}_ai{}_ao{}_nstars{}".format(
    single_detection_prob, int(a_inner_true), int(a_outer_true), n_stars)

    def suppression_factor(a_values):
        results = (a_values - np.log10(a_inner_true)) / (np.log10(a_outer_true) - np.log10(a_inner_true))
        return np.clip(results, a_min=0, a_max=1)

if influencer == 'extrasuppression':
    a_inner_true = 10  # AU (suppression 100%)
    a_outer_true = 1000  # AU (suppression 0%)

    working_dir = "results/extrasuppression_p{}_ai{}_ao{}_nstars{}".format(
    single_detection_prob, int(a_inner_true), int(a_outer_true), n_stars)

    def suppression_factor(a_values):
        results = (a_values - np.log10(a_inner_true)) / (np.log10(a_outer_true) - np.log10(a_inner_true))
        return np.clip(results, a_min=0, a_max=1)

elif influencer == 'encouragement':
    a_max_true = 10  # AU (suppression 100%)
    a_outer_true = 200  # AU (suppression 0%)

    working_dir = "results/encouragement_p{}_ai{}_ao{}_nstars{}".format(
    single_detection_prob, int(a_max_true), int(a_outer_true), n_stars)

    def suppression_factor(a_values):
        results = np.ones_like(a_values)
        results[a_values<np.log10(a_outer_true)] = 1.6
        return results

        #results = 1.6 - np.abs( (a_values - np.log10(a_max_true)) )
        #results[(a_values>1) & (results<1)] = 1
        #return np.clip(results, a_min=0, a_max=2)

else:
    raise

if not os.path.exists(working_dir):
    os.mkdir(working_dir)

# draw a sample of binaries, table 2 from Offner+2023 PP proceedings
a_values = truncnorm.rvs(trunc_low, trunc_high,
                                loc=logmu, scale=logsigma,
                                size=n_binaries)

# assign each star a "detection" or "nondetection" of planet
# with presence of planets suppressed by a function of semi-major axis
my_factor = suppression_factor(a_values)
binary_detection_prob = single_detection_prob*my_factor
binary_detections = np.random.rand(n_binaries) < binary_detection_prob
#single_detections = np.random.rand(n_singles) < single_detection_prob


fig, ax = plt.subplots()
plt.plot(a_values, binary_detection_prob,'.')
plt.ylim(-0.01,0.51)
plt.xlabel("log(a) [au]")
plt.ylabel("suppressed/enhanced occurrence rate")
plt.savefig("{}/suppression_factor.png".format(working_dir), dpi=250)
plt.close()


fig, ax = plt.subplots(1, 1)
plt.hist(a_values, range=(1,3),
         histtype='step', linestyle='-', color='k',
         label='simulated targets',
         zorder=1)
plt.hist(a_values[binary_detections],range=(1,3),
         alpha=0.8, color='indianred',
         label='planet detected',
         zorder=1)
plt.legend()
plt.ylabel("no. of stars")
plt.xlabel(r"$log_{10}(a)$ [au]")

ax2 = ax.twinx()
bin_edges = [1.,1.5, 2., 2.5, 3.]
bin_centers = bin_edges[:-1] + np.diff(bin_edges)/2.
bin_widths = np.diff(bin_edges)/2.
count_all, _ = np.histogram(a_values, bins=bin_edges)
count_det, _ = np.histogram(a_values[binary_detections], bins=bin_edges)
f_obs = count_det/count_all
ax2.errorbar(bin_centers, f_obs,
             yerr=np.sqrt(count_det)/count_all,
             xerr=bin_widths,
             capsize=3,
             linestyle='None',
             zorder=2)
ax2.set_ylabel('detection fraction')
plt.ylim(0,1)

# add line at fiducial occurrence rate
fid_err = np.sqrt(25)/100 # assume look at 100 stars and find 25 earths
ax2.axhspan(0.25-fid_err, 0.25+fid_err, ec='k', fc='lightgray', alpha=0.3, ls=':', zorder=0)

plt.tight_layout()

#ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
#ax.patch.set_visible(False) # hide the 'canvas'

plt.savefig("{}/seed{}_sample.png".format(working_dir, random_seed), dpi=250)

with open("{}/seed{}_data.txt".format(working_dir, random_seed), 'w') as f:
    for aa, det, all in zip(bin_centers, count_det, count_all):
        f.write("{:.3f}\t{}\t{}\n".format(aa, det, all))

with open("{}/seed{}_data.txt".format(working_dir, random_seed), 'w') as f:
    for aa, det, all in zip(bin_centers, count_det, count_all):
        f.write("{:.3f}\t{}\t{}\n".format(aa, det, all))

from scipy.stats import binomtest

with open("{}/seed{}_stats.txt".format(working_dir, random_seed), 'w') as f:
    # chi2 = chisquare(f_obs=count_det, f_exp=count_all*0.25) # doesn't work, must have all f_obs = f_exp
    # Perform the binomial test
    for i in range(len(f_obs)):
        print(count_det[i], count_all[i], count_all[i]*single_detection_prob )
        result = binomtest(count_det[i], count_all[i], single_detection_prob, alternative='two-sided')
        f.write("{}\t{}\n".format(bin_centers[i], result.pvalue))





"""
define the likelihood
"""


# def single_star_loglike(params_arr, star_age, star_age_unc, star_meas_o3):
#     mu_pop = params_arr[0]
#     sigma_pop = params_arr[1]
#
#     # O3 measured in star: we want the prob. that we would measure O3
#     # in this star given the proposal distribution and the age unc.
#     if star_meas_o3:
#         # likelihood = prob. that age (characterized by measurement unc. distribution)
#         # is greater than population distribution (characterized by proposal distribution).
#         # This is the same as the prob. that the difference of the two Gaussians (which is
#         # also Gaussian) is greater than 0
#         loglike = norm.logcdf(
#             0, loc=mu_pop - star_age, scale=np.sqrt(sigma_pop**2 + star_age_unc**2)
#         )
#
#     # no O3 measured (opposite of above)
#     else:
#         loglike = norm.logcdf(
#             0, loc=star_age - mu_pop, scale=np.sqrt(sigma_pop**2 + star_age_unc**2)
#         )
#
#     return loglike
#
#
# def loglike(params_arr, star_ages_arr, star_age_unc_arr, star_meas_o3_arr):
#     n_stars = len(star_ages_arr)
#     loglike = 0
#     for i in range(n_stars):
#         loglike += single_star_loglike(
#             params_arr, star_ages_arr[i], star_age_unc_arr[i], star_meas_o3_arr[i]
#         )
#
#     return loglike
#
#
# def logposterior(params_arr, star_ages_arr, star_age_unc_arr, star_meas_o3_arr):
#     mu_pop = params_arr[0]
#     sigma_pop = params_arr[1]
#
#     # set priors on pop-level parameters
#     if mu_pop < 0 or mu_pop > 13:
#         return -np.inf
#     if sigma_pop < 0 or sigma_pop > 13:
#         return -np.inf
#
#     return loglike(params_arr, star_ages_arr, star_age_unc_arr, star_meas_o3_arr)
#
#
# """
# run mcmc!
# """
#
# ndim, nwalkers = 2, 100
#
# p0 = np.random.uniform(0, 13, size=(nwalkers, ndim))
# sampler = emcee.EnsembleSampler(
#     nwalkers, ndim, logposterior, args=(sampled_ages, age_uncertainties, o3_detections)
# )
#
# print("running burn in!")
# state = sampler.run_mcmc(p0, 200)
# sampler.reset()
# print("running production chain!")
# sampler.run_mcmc(p0, 500)
# print("done")
#
# samples = sampler.flatchain
# mu_samples = samples[:, 0]
# sigma_samples = samples[:, 1]
# np.save(f"{working_dir}/seed{random_seed}_chains.npy", samples)
#
# """
# plot results
# """
#
# fig, ax = plt.subplots(2, 1)
# ax[0].hist(
#     mu_samples, bins=50, color="rebeccapurple", alpha=0.5, label="recovered posterior"
# )
# ax[0].axvline(mu_true, color="k", label="true value")
# ax[0].legend()
# ax[0].set_xlabel("$\\mu_{\\mathrm{{population}}}$ [Gyr]")
# ax[1].set_xlabel(("$\\sigma_{\\mathrm{{population}}}$ [Gyr]"))
# ax[1].hist(sigma_samples, bins=50, color="rebeccapurple", alpha=0.5)
# ax[1].axvline(sigma_true, color="k")
# plt.tight_layout()
# plt.savefig("{}/seed{}_recovery.png".format(working_dir, random_seed), dpi=250)
