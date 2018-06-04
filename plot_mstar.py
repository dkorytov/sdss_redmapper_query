#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo,get_current
import astropy
from background import *

zs = np.linspace(0,.35,100)

plt.figure()
plt.plot(zs,mstar(zs)-cosmo.distmod(zs).value)
plt.xlabel('redshift')
plt.ylabel('Absolute Mag')
plt.grid()


plt.figure()
plt.plot(zs,mstar(zs))
plt.xlabel('redshift')
plt.ylabel('Apparent Mag')
plt.grid()


plt.figure()
plt.plot(zs,-cosmo.distmod(zs).value)
plt.xlabel('redshift')
plt.ylabel('Distance Modulus')
plt.grid()


plt.show()


