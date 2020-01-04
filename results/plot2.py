#!/usr/bin/env python3

import matplotlib

#matplotlib.use('Qt4Agg')

import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import make_interp_spline, BSpline, splev, splrep

fig = plt.figure()
ax = fig.add_subplot(111)
s = ['beta = 50', 'beta = 100', 'beta = 250', 'beta = 500', 'beta = 1000']
print(sys.argv)
for f in range(len(sys.argv) - 1):
    data = np.genfromtxt(str(sys.argv[f + 1]) + '', delimiter=' ')

    x = data[:, 0].astype(float)
    col = 'bgrcmyk'

    xnew = np.linspace(x.min(),x.max(),300)

    #y = data[:, 2].astype(float)
    ##spl = make_interp_spline(x, y, k = 3)
    #spl = splrep(x, y, s = 5)
    #spl_y = splev(x, spl)
    #poly = np.polyfit(x, y, 5)
    #poly_y = np.poly1d(poly)(x)
    #ax.plot(x, poly_y, color=col[1], label='Max')

    y = data[:, 1].astype(float)
    spl = make_interp_spline(x, y, k = 3)
    ax.plot(xnew, spl(xnew), color=col[f], label=s[f])

    #i, l = 0, 0
    #while l < len(x) - 2:
    #    i += 1
    #    beta, t = x[l], y[l]
    #    l = l + 1
    #    z = int(t / 100) + 2
    #    ax.plot(x[l: l + z], y[l: l + z], color=col[i%len(col)], label=str(beta))
    #    l += z
    #

ax.set_title('time vs err')
ax.set_xlabel('time (in s)')
ax.set_ylabel('err')
ax.legend()
#plt.show()
plt.savefig('all.png')
plt.close()
