#! /usr/bin/env python
# -*- coding: utf-8 -*-
  
"""Detection of sources."""

import errno, glob, os, shutil, subprocess, sys
import numpy as np
from astropy.io import fits
from .utils import mv_p, mkdir_p


def psfex(filename, useweight=False):
    """Compute PSF in astronomical images"""
    
    #imagelist=glob.glob(path+'/*.fits')
    imagelist = np.atleast_1d(filename)
    print (imagelist)
    for ima in imagelist:
         root = os.path.splitext(ima)[0]
         if useweight:
             weight = root + '.weight.fits'
             #print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima, '-WEIGHT_IMAGE', weight])
         else:
             print ('Sextractor')
             #print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima])

         cat = 'prepsfex.cat'
         print ('psfex')
         #subprocess.call(['psfex', cat, '-c', 'config.psfex'])
         subprocess.call(['psfex', cat])
         mv_p('snap_prepsfex.fits', root + '.psf.fits')

def sextractor(path, outdir='sources_cat/', useweight=False):
    """Detect sources in astronomical images"""
    
    imagelist=glob.glob(path+'/*.fits')
    mkdir_p(outdir)

    for ima in imagelist:
        if '.psf' not in ima:
            root = os.path.splitext(ima)[0]
            if useweight:
                weight = root + '.weight.fits'
                #print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima, '-WEIGHT_IMAGE', weight])
            else:
                print (ima.split('/')[-1].split('.')[0])
                #print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima,
                                 '-CATALOG_NAME', outdir+ima.split('/')[-1].split('.')[0]+'.cat'])

def psfex(filename, telescope):
    """Run psfex to estimate PSF FWHM"""

    subprocess.call(['sex', '-c', 'config/%s/sourcesdet.sex' % telescope, filename, '-SEEING_FWHM', str(fwhmpsf), '-DETECT_THRESH', str(thresh), '-PARAMETERS_NAME', 'config/%s/sourcesdet.param' % telescope])

def run_psfex(filename):
    psfex(filename)
    psffilename = filename.split('.')[0] + '.psf.fits'
