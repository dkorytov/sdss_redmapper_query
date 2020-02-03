#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
# import pyfits
import astropy.io.fits as fits
import dtk 

def create_miscentering_distribution(fname):
    hdulist = fits.open(fname)
    print(hdulist)
    tdata = hdulist[1].data
    cdata = hdulist[1].columns
    print(tdata)
    print(hdulist.info())
    print(cdata)

if __name__ == '__main__':
    create_miscentering_distribution(sys.argv[1])
    plt.show()
