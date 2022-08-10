###
# Minimal code to plot the sample of free eccentricities
###
import numpy as np
import matplotlib.pyplot as plt

npbar = 1000
nw = 5000
[nej, efree] = np.loadtxt("RandomGrowth_FreeEccentricity_Levi_npbar%s_nw%s.dat" % (npbar, nw))
# nej = Np in the paper, efree is in the unit of the final forced eccentricity

print(f"{len(nej)} chains")

plt.scatter(nej, efree, s=0.1)

plt.axvline(x=npbar, color='k', label=r"$\bar{N}_{\rm p}$")

plt.xlabel(r"$N_{\rm p}$")
plt.ylabel(r"$e_{\rm free} / e_{\rm forced}$")
plt.legend()
plt.ylim(ymin=0)
plt.margins(0)

plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)
plt.show()
