#! /usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess, sys, os
import numpy as np
from astropy.io import fits
from lacosmic import lacosmic
from utils import cp_p

def run_lacosmic(filename, FWHM, flim=2, sigma=5, outLevel=1):
    """Run lacosmic to remove cosmic rays from the input image"""

    imagelist = np.atleast_1d(filename)

    for i, ima in enumerate(imagelist):
        path, filename_ext = os.path.split(ima)
        if path:
            folder = path + '/'
        else:
            folder = ''

        filename2 = filename_ext.split('.')[0]

        # Make copy of original image
        if outLevel == 2:
            cp_p(ima, folder+filename2+'_CR_notcleaned.fits')

        hdulist = fits.open(ima)
        hdr = hdulist[0].header
        try:
            gain = hdr['GAIN']
        except:
            gain = 1
        try:
            RN = hdr['RN']
        except:
            RN = 10

        if FWHM[i] > 2:
            flim = 2
        else:
            flim = 5
        data = np.asarray(hdulist[0].data, dtype=float)
        lacosmic_res = lacosmic(data,flim,sigma,sigma,effective_gain=gain,readnoise=RN)

        # Create image cleaned from cosmic rays
        hdulist[0].data = lacosmic_res[0]
        hdulist.writeto(ima, overwrite=True)

        if outLevel == 2:
            # Create mask of cosmic rays
            hdulist[0].data = np.asarray(lacosmic_res[1],dtype=int)
            hdulist.writeto(folder+filename2+'_CRmask.fits', overwrite=True)

