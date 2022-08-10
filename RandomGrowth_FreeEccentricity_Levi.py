###
# Creates a sample of planet-planet scattering evolutions of different sizes following a Levi distribution.
###
from numpy import exp, cos, sin, sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc, erfcinv

eforced = 1 # final forced eccentricity is set to 1 as the unit
npbar = 1000
s = sqrt(eforced*eforced/(2*npbar)) # scale of the kicks
print(f"Walk step s = {s:.2f}")

nw = 5000 # number of kicks within one precession period Nw
wdt = 2*np.pi/nw

nwalks = 10000 # number of generated walks

# Distribution of the lengths of the walks
npmax = int(10*npbar)
npmin = int(np.floor(eforced/s)+1)
print(f"Minimum length of the chains: {npmin}")
print(f"Maximum length of the chains: {npmax}")
def f(n, b): #Levy ditribution
    return(b/sqrt(2*np.pi*(n**3)) * exp(-b*b/(2*n)))
b = sqrt(3*npbar) # b is defined such that npbar is the most likely value
nrange = np.arange(npmin, npmax+1)

# See Appendix B for the sampling of the Levi distribution
pmin = erfc(b/sqrt(2*npmin))
pmax = erfc(b/sqrt(2*npmax))
p = pmin + np.random.random(nwalks)*(pmax-pmin) #uniform sampling
n = np.round(b*b/(2*erfcinv(p)**2),0).astype("int") #n follows f
nwalks_per_n = np.array([np.sum(n==nn) for nn in range(0, npmax+1)])

nwalks = int(np.sum(nwalks_per_n))
print(f"New number of chains = {nwalks}") # can be a bit different from the desired value

x = np.zeros(nwalks) # ep cos(wp)
y = np.zeros(nwalks) # ep sin(wp)

nej = np.zeros(nwalks, dtype='int')
efree = np.zeros(nwalks)
ctry = 0
c = 0
cpern = np.zeros(npmax+1, dtype='int')
while np.any(cpern < nwalks_per_n):
    x = np.zeros(npmax+1)
    y = np.zeros(npmax+1)
    for n in range(1, npmax+1):

        dx = np.random.normal(scale=s)
        dy = np.random.normal(scale=s)

        x[n] = x[n-1] + dx
        y[n] = y[n-1] + dy

        ep = sqrt(x[n]**2 + y[n]**2)

        if (cpern[n] < nwalks_per_n[n]) and (ep > 1-s/2) and (ep < 1+s/2):
            nej[c] = n

            # Rotation so that the walks starts at x+iy = 0 and ends at x+iy = 1
            theta = np.arctan2(y[n], x[n])
            xrot = (cos(theta)*x + sin(theta)*y)[:(n+1)]
            yrot = (-sin(theta)*x + cos(theta)*y)[:(n+1)]

            wtsteps = np.arange(n+1)*wdt
            ksteps = np.zeros(n+1) # e cos(w)
            hsteps = np.zeros(n+1) # e sin(w)

            for j in range(n+1):
                if j == 0:
                    ksteps[j] = 0.
                    hsteps[j] = 0.
                else:
                    wt0 = wtsteps[j-1]
                    k0 = ksteps[j-1]
                    h0 = hsteps[j-1]
                    wt = wtsteps[j]-wt0
                    ksteps[j] = k0*cos(wt) - h0*sin(wt) + xrot[j-1]*(1-cos(wt)) + yrot[j-1]*sin(wt)
                    hsteps[j] = k0*sin(wt) + h0*cos(wt) - xrot[j-1]*sin(wt) + yrot[j-1]*(1-cos(wt))

            kf = ksteps[-1]
            hf = hsteps[-1]

            efree[c] = sqrt((kf-xrot[-1])**2+(hf-yrot[-1])**2)

            c += 1
            cpern[n] += 1
            if c % int(nwalks/10) == 0:
                print(f"{c} walks")
    ctry += 1

print(f"Successful walks: {c}/{ctry}")

np.savetxt("RandomGrowth_FreeEccentricity_Levi_npbar%s_nw%s.dat" % (npbar, nw), [nej,efree])

### Plot the distribution of lengths and the Levi distribution it derives from
plt.hist(nej, histtype="step", bins='auto', density=True)
plt.plot(nrange, f(nrange,b))
plt.plot(nrange, f(nrange,b)/erfc(np.sqrt(3/5)/2)) # renormalized Levi distribution with nmax = 10 npbar
plt.axvline(x=npbar, color='k', label=r"$\bar{N}_{\rm p}$")
plt.margins(x=0)
plt.legend()
plt.xlabel(r"$N_{\rm p}$")
plt.ylabel("Number of walks")

plt.tight_layout()
plt.show()
