import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import emcee
import os
import sys
from scipy.stats import truncnorm
random_seed = 1

np.random.seed(random_seed)

## variables: sigma_true, n_stars, age_unc
## make one plot for sigma_true = 1, one plot for sigma_true = 0.1

"""
assume true values for the inner and outer semi-major axis where binary suppression matters
"""

influencer = sys.argv[1]

mu = 40 # AU
sigma = 1.5 # log10(AU)
max_loga = 5 # 10^5 AU is max a shown in Offner plot
max_value = (max_loga-np.log10(mu))/1.5 # truncates at no. of sigma from loc

min_loga = np.log10(10) # assume we won't observe anything with companion closer than 10 AU
min_value = (min_loga-np.log10(mu))/1.5 # truncates at no. of sigma from loc

single_detection_prob = 0.25 # prob of whatever planet I am simulating existing around a single star

"""
draw a dataset
"""

n_stars = int(sys.argv[2])
MF = 1 #0.47 # multiplicity fraction
n_binaries = int(n_stars*MF)
#n_singles = n_stars - n_binaries


if influencer == 'suppression':
    a_inner_true = 100  # AU (suppression 100%)
    a_outer_true = 1000  # AU (suppression 0%)

    working_dir = "results/suppression_ai{}_ao{}_nstars{}".format(
    int(a_inner_true), int(a_outer_true), n_stars)

    def suppression_factor(a_values):
        results = (a_values - np.log10(a_inner_true)) / (np.log10(a_outer_true) - np.log10(a_inner_true))
        return np.clip(results, a_min=0, a_max=1)

elif influencer == 'encouragement':
    a_max_true = 500  # AU (suppression 100%)
    a_outer_true = 1500  # AU (suppression 0%)

    working_dir = "results/encouragement_ai{}_ao{}_nstars{}".format(
    int(a_max_true), int(a_outer_true), n_stars)

    def suppression_factor(a_values):
        results = 1.75 - np.abs( (a_values - np.log10(a_max_true)) )
        results[(a_values>3) & (results<1)] = 1
        return np.clip(results, a_min=0, a_max=2)

else:
    raise

if not os.path.exists(working_dir):
    os.mkdir(working_dir)

# draw a sample of binaries, table 2 from Offner+2023 PP proceedings
a_values = truncnorm.rvs(min_value, max_value,
                                loc=np.log10(mu), scale=sigma,
                                size=n_binaries)

# assign each star a "detection" or "nondetection" of planet
# with presence of planets suppressed by a function of semi-major axis
my_factor = suppression_factor(a_values)
fig, ax = plt.subplots()
plt.plot(a_values, my_factor,'.')
plt.xlabel("log(a) [au]")
plt.ylabel("suppression factor")
plt.savefig("{}/suppression_factor.png".format(working_dir), dpi=250)
plt.close()

binary_detection_prob = single_detection_prob*my_factor
binary_detections = np.random.rand(n_binaries) < binary_detection_prob

#single_detections = np.random.rand(n_singles) < single_detection_prob

fig, ax = plt.subplots(1, 1)
plt.hist(a_values, range=(1,5),
         histtype='step', linestyle='-', color='k',
         label='simulated targets')
plt.hist(a_values[binary_detections],range=(1,5),
         alpha=0.5, color='indianred',
         label='planet detected')
plt.legend()
plt.ylabel("no. of stars")
plt.xlabel(r"$log_{10}(a)$ [au]")

ax2 = ax.twinx()
bin_edges = [1., 2., 3., 4., 5.]
bin_centers = bin_edges[:-1] + np.diff(bin_edges)/2.
count_all, _ = np.histogram(a_values, bins=bin_edges)
count_det, _ = np.histogram(a_values[binary_detections], bins=bin_edges)
ax2.errorbar(bin_centers, count_det/count_all, np.sqrt(count_det)/count_all, linestyle='--')
ax2.set_ylabel('detection fraction')
plt.ylim(0,1)

plt.tight_layout()
plt.savefig("{}/seed{}_sample.png".format(working_dir, random_seed), dpi=250)






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
