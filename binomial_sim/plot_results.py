
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib as mpl

mpl.rc('lines',linewidth = 1.5)
mpl.rc('font',size = 12)
mpl.rc('axes',labelsize = 14, linewidth=1.25)
mpl.rc('xtick',labelsize = 14)
mpl.rc('ytick',labelsize = 14)
# enable math fonts
mpl.rc('mathtext', default = 'regular')

n_stars = np.array([50, 100, 200, 400])
alphas = [0.2, 0.4, 0.6, 0.8]
bin_centers = [1.25, 1.75, 2.25, 2.75]
influencer = 'extrasuppression_p0.25_ai10_ao1000'

fig, ax = plt.subplots(1, 1, sharex=True)
plt.subplots_adjust(hspace=0)

for j, n_star in enumerate(n_stars):

    files = glob.glob('results/{}_nstars{}/seed*_stats.txt'.format(influencer, n_star))
    pvalues = np.zeros((len(files),len(bin_centers)))
    for i, myfile in enumerate(files):
        output = np.genfromtxt(myfile).transpose()
        pvalues[i] = output[1]

        ax.plot(
                bin_centers,
                pvalues[i],
                #label=f"{age} Gyr age unc.",
                alpha=alphas[j],
                color="rebeccapurple")

        print(myfile, pvalues[i])

    print(np.mean(pvalues, axis=0))
    plt.yscale('log')
    plt.ylim(1e-9,5)
    plt.axhline(0.05, linestyle='--', c='gray', zorder=0)
    plt.xlabel(r"$log_{10}(a)$ [au]")
    plt.ylabel('binomial statistic p-value')


plt.tight_layout()
plt.savefig('results/{}_summary.png'.format(influencer, n_star), dpi=250)
plt.close()


