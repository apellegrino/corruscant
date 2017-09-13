import numpy as np
import matplotlib.pyplot as plt
import richardsplot
from scipy.stats import linregress

def analytic(n, r):
    return n * (1 - 2*r)**3 * (4. / 3. * np.pi * n * r**3 + 1)

#with open("analytic.csv", 'r') as f:
data = np.loadtxt("analytic.csv", delimiter=',')
n = data[:,0]
trials = data[:,1:]

fig = plt.figure(figsize=(5,5))
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
#scale_err = np.std(trials, axis=1)/analytic(n, 0.05)
scale_err = (np.average(trials, axis=1) - analytic(n, 0.05))/analytic(n, 0.05)
scale_err = np.abs(scale_err)
ax.scatter(n, scale_err,
        color='b')

#mask = np.array([False, True, False, True, False, False, False, False, True, True, True])
#regression = linregress(np.log10(n[mask]), np.log10(scale_err[mask]))

regression = linregress(np.log10(n), np.log10(scale_err))
print regression
print "R^2 = {:f}".format(regression.rvalue**2)
ax.plot(n, 10**regression.intercept * n ** regression.slope, color='k', zorder=-99)
#ax.plot(n, 1.2 * n ** -0.5, color='k', zorder=-99)
ax.set_ylabel("$\\sigma_{numerical}/f(N)_{analytic}$")
ax.set_xlabel("$N$")
plt.tight_layout()
plt.savefig("errcomplexity.png", dpi=200)
#plt.show()

fig = plt.figure(figsize=(5,5))
ax = plt.gca()
#ax.plot(n, analytic(n, 0.05))
ax.errorbar(n, np.average(trials, axis=1)/analytic(n, 0.05),
            yerr=np.std(trials, axis=1)/analytic(n, 0.05),
            color='b', marker='o', capsize=2.)
ax.set_xscale("log")
#ax.set_yscale("log")
#ax.set_ylim(0,ax.get_ylim()[1])

xl = ax.get_xlim()
ax.plot(xl, (1,1), ls='dashed', color='#cccccc', zorder=-99)
ax.set_xlim(xl)

ax.set_xlabel("$N$")
ax.set_ylabel("Numerical / Analytic Ratio")

plt.tight_layout()
plt.savefig("convergence.png", dpi=200)
#plt.show()
