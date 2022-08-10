###
# Creates a planet-planet scattering evolution and the consequence on an outer test particle (Figure 5)
###
from numpy import cos, sin
import numpy as np
import matplotlib.pyplot as plt

# eforced = 1 Final forced eccentricity is set to 1 as the unit

# First let us generate the trajectory of the complex eccentricity of the planet times nu/omega (or equivalently the trajectory of the forced eccentricity)
nej = 1000 # Number of steps Np
s = np.sqrt(1/(2*nej)) # Kick scale so that a chain of length nej has a likely chance of reaching 1 (\bar{Np} = Np)
print(f"Walk step s = {s:.2f}")

x = [0] # ep cos(wp)
y = [0] # ep sin(wp)
ep = 0
# We generate random walks of size nej starting at e=0 until one ends with e within s of epf:
while (ep < 1-s/2) or (ep > 1+s/2):

    dx = np.random.normal(scale=s, size=nej)
    dy = np.random.normal(scale=s, size=nej)

    x = np.array([np.sum(dx[:j]) for j in range(nej+1)])
    y = np.array([np.sum(dy[:j]) for j in range(nej+1)])
    ep = np.sqrt(x[-1]**2 + y[-1]**2)

# Rotation so that the walks starts at x+iy = 0 and ends at x+iy = 1
theta = np.arctan2(y[-1], x[-1])
xrot = cos(theta)*x + sin(theta)*y
yrot = -sin(theta)*x + cos(theta)*y
x = xrot
y = yrot

# Consequence on an outer test particle
nw = 1000 # number of kicks within one precession period Nw
wdt = 2*np.pi/nw

# The test particle's complex eccentricity follows a circular trajectory around center x+iy
nt = 10*nej # Number of snapshot in the evolution of the test particle
wtej = wdt*nej
wtrange = np.linspace(0, wtej, nt)

k = np.zeros(nt) # e cos(w)
h = np.zeros(nt) # e sin(w)
for j, wtt in enumerate(wtrange):
    if j == 0:
        k[j] = 0
        h[j] = 0
    else:
        n = int(wtt/wdt)
        wt0 = wtrange[j-1]
        k0 = k[j-1]
        h0 = h[j-1]

        wt = wtt-wt0
        k[j] = k0*cos(wt) - h0*sin(wt) + x[n]*(1-cos(wt)) + y[n]*sin(wt)
        h[j] = k0*sin(wt) + h0*cos(wt) - x[n]*sin(wt) + y[n]*(1-cos(wt))

# Representation on the complex plane

# Final free eccentricity and representation
efree = np.sqrt((k[-1]-x[-1])**2+(h[-1]-y[-1])**2)
print(f"efree={efree}")
plt.plot([1, x[-1]+efree], [0, 0], "C1")
plt.text(0.5*(1+x[-1]+efree), 0, r"$e_{\rm free}$", color="C1", ha="center", va="bottom", fontsize=10)

# Representation of the final circular trajectory of the test particle eccentricity
psi = np.arctan2(h[-1]-y[-1],k[-1]-x[-1])
phi = np.linspace(0, 2*np.pi, 1000)
keq = efree*np.cos(phi+psi) + x[-1]
heq = efree*np.sin(phi+psi) + y[-1]
plt.plot(keq, heq, "C1--")

# Representation of the trajectory of the test particle eccentricity before reaching the final circle
plt.plot(k, h, "C0.-", label=r"$\mathcal{E}(t)$")

# Representation of the trajectory of the test particle eccentricity when reaching the final circle
phi = np.linspace(0, np.pi/2, 1000)
k = efree*np.cos(phi+psi) + x[-1]
h = efree*np.sin(phi+psi) + y[-1]
plt.plot(k, h, "C0.-")

# Representation of the trajectory of the planet eccentricity
plt.plot(x, y, "k.-", label=r"$\frac{\nu}{\omega} \mathcal{E}_{\rm p}(t)$")
plt.plot(1, 0, "C1o") # Forced eccentricity

plt.axhline(y=0, lw=0.1, color="k")
plt.axvline(x=0, lw=0.1, color="k")
plt.xlabel(r"$e~\cos(\varpi-\varpi_{\rm p})$ [$e_\mathrm{forced}$]")
plt.ylabel(r"$e~\sin(\varpi-\varpi_{\rm p})$ [$e_\mathrm{forced}$]")
plt.legend()
plt.axis("equal")

plt.tight_layout()
plt.show()
