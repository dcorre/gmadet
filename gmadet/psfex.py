#! /usr/bin/env python
# -*- coding: utf-8 -*-
  
"""Detection of sources."""

import errno, glob, os, shutil, subprocess, sys
import numpy as np
from astropy.io import fits
from utils import mv_p, mkdir_p
import xmltodict

def psfex(filename, config, useweight=False):
    """Compute PSF in astronomical images"""

    FWHM_list = []

    #imagelist=glob.glob(path+'/*.fits')
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
        print ('\nRunning psfex to estimate FWHM in %s'  % ima)
        root = os.path.splitext(ima)[0]
        if useweight:
            weight = root + '.weight.fits'
            subprocess.call(['sex', '-c', config['psfex']['sextractor'], \
                    ima, \
                    '-WEIGHT_IMAGE', weight, \
                    '-PARAMETERS_NAME', config['psfex']['param']])
        else:
            subprocess.call(['sex', '-c', config['psfex']['sextractor'], \
                    ima, \
                    '-PARAMETERS_NAME', config['psfex']['param']])

        cat = 'prepsfex.cat'
        subprocess.call(['psfex', cat, '-c', config['psfex']['conf'] ])
        mv_p('snap_prepsfex.fits', root + '.psf.fits')

    # Get the mean PSF FWHM in pixels
    with open('psfex.xml') as fd:
        doc = xmltodict.parse(fd.read())
        FWHM_stats = doc['VOTABLE']['RESOURCE']['RESOURCE']['TABLE'][0]['DATA']['TABLEDATA']['TR']['TD'][20:23]
        FHWM_min = float(FWHM_stats[0])
        FHWM_mean = float(FWHM_stats[1])
        FHWM_max = float(FWHM_stats[2])
   
        print ('\nFWHM min: %.2f pixels' % FHWM_min)
        print ('FWHM mean: %.2f pixels' % FHWM_mean) 
        print ('FWHM max: %.2f pixels\n' % FHWM_max)
    
        os.remove('psfex.xml')
        FWHM_list.append(FHWM_mean)

    return FWHM_list

