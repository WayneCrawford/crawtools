"""
Compare three methods for calculating tanh(x)
"""
import numpy as np
from matplotlib import pyplot as plt

def tanh_mosher(x):
    xsq = np.power(x, 2.)
    y = x*(27 + xsq)/(27 + 9*xsq)
    y[x>2.96] = 1
    return y
    
def dtanh(x):
    a = np.exp(x*(x <= 50))
    one = np.ones(np.shape(x))

    y = (abs(x) > 50) * (abs(x)/x) + (abs(x) <= 50)*((a-one/a) / (a+one/a))
    return y

x = np.power(10., np.arange(-2, 2, 0.01))
fig, ax = plt.subplots(2, 1, sharex=True)
dta = dtanh(x)
mos = tanh_mosher(x)
npt = np.tanh(x)
ax[0].semilogx(x, mos, label='mosher')
ax[0].semilogx(x, dta, label='dtanh')
ax[0].semilogx(x, npt, label='np.tanh')
ax[0].set_ylabel('tanh(x)')
ax[0].legend()
ax[1].loglog(x, abs(dta-npt)/dta, label='abs(dtanh-np.tanh)/dtanh')
ax[1].loglog(x, abs(dta-mos)/dta, label='abs(dtanh-mosher)/dtanh')
ax[1].set_ylabel('diff')
ax[1].set_xlabel('x')
ax[1].legend()
plt.show()