#!/usr/bin/env python3

import matplotlib
matplotlib.use('Qt4Agg')

import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import make_interp_spline, BSpline

data = np.genfromtxt(str(sys.argv[1]) + '', delimiter=', ')
fig = plt.figure()
ax = fig.add_subplot(111)

x = data[:, 0].astype(int)
col = 'bgrcmyk'

xnew = np.linspace(x.min(),x.max(),300)

#s = ['k=2', 'k=4', 'k=8', 'k=16', 'k=32', 'k=50', 'k=62']
s = ['n=4941', 'n=5000', 'n=8000', 'n=10680', 'n=15000']
for i in range(1, len(s) + 1):
    if i==1: continue
    y = data[:, i].astype(float)
    print(x, y)
    spl = make_interp_spline(x, y, k = 3)
    #ax.plot(xnew, spl(xnew), color=col[i], label=s[i - 1])
    ax.plot(x, y, color=col[i%len(col)], label=s[i-1])

#i, l = 0, 0
#while l < len(x) - 2:
#    i += 1
#    beta, t = x[l], y[l]
#    l = l + 1
#    z = int(t / 100) + 2
#    ax.plot(x[l: l + z], y[l: l + z], color=col[i%len(col)], label=str(beta))
#    l += z
#

ax.set_title('k vs time')
ax.set_ylabel('time (in s)')
ax.set_xlabel('k (steps)')
ax.legend()
#plt.show()
plt.savefig(sys.argv[1] + '.png')
plt.close()
