#!/usr/bin/env python

import numpy as np
import matplotlib as plt

def get_mstar(z):
    return 12.27+62.36*z-289.79*z**2+729.69*z**3-709.42*z**4


def phi(m,z):
    mstar = get_mstar(z)
    return 10.0**(-.4*(m-mstar)*(1.8))*np.exp(-10.0**(-.4*(m-mstar)))
