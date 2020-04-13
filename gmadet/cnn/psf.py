#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:36:04 2019

@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, math, os, shutil, subprocess, sys
import argparse
from gmadet.utils import getpath,load_config, mv_p, rm_p
from astropy.io import fits
import xmltodict

def compute_psf(path, telescope, useweight):
    """ Compute psf """

    path_gmadet = getpath()

    config = load_config(telescope)
    imagelist = sorted(glob.glob(path + '/**/*.fit*', recursive=True))

    useweight = bool(useweight)

    filtername = path_gmadet + '/config/conv_kernels/default.conv'

    for ima in imagelist:
         root = os.path.splitext(ima)[0]
         if useweight:
             weight = root + '.weight.fits'
             print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', ima, \
                              '-c', config['psfex']['sextractor'], \
                              '-PARAMETERS_NAME', config['psfex']['param'], \
                              '-FILTER_NAME', filtername, \
                              '-WEIGHT_IMAGE', weight])
         else:
             print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', ima, \
                              '-c', config['psfex']['sextractor'], \
                              '-PARAMETERS_NAME', config['psfex']['param'], \
                              '-FILTER_NAME', filtername])
         cat = 'preppsfex.cat'
         subprocess.call(['psfex', cat, '-c', config['psfex']['conf'] ])
         mv_p('snap_preppsfex.fits', root + '.psf.fits')
         # Get the mean PSF FWHM in pixels
         with open('psfex.xml') as fd:
            doc = xmltodict.parse(fd.read())
            FWHM_stats = doc['VOTABLE']['RESOURCE']['RESOURCE']['TABLE'][0]['DATA']['TABLEDATA']['TR']['TD'][20:23]
            FHWM_min = float(FWHM_stats[0])
            FHWM_mean = float(FWHM_stats[1])
            FHWM_max = float(FWHM_stats[2])
         # Get number of psf snapshot per axis
         psf_snaps = os.popen("sed -n '/PSFVAR_NSNAP/p' %s" % config['psfex']['conf']).read()
         nb_snaps = psf_snaps.split()[1]

         hdulist = fits.open(root+'.psf.fits')
         hdr = hdulist[0].header
         hdr['FWHMMIN'] = str(FHWM_min)
         hdr['FWHMMEA'] = str(FHWM_mean)
         hdr['FWHMMAX'] = str(FHWM_max)
         hdr['PSF_NB'] = str(nb_snaps)
         hdulist.writeto(root+'.psf.fits',overwrite=True)

         rm_p(cat)
         rm_p('preppsf.psf')
         rm_p('psfex.xml')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Compute PSF for images in a given repository.')

    parser.add_argument('--path',
                        dest='path',
                        required=True,
                        type=str,
                        help='Path to images')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--useweight',
                        dest='useweight',
                        action='store_true',
                        help='If set, use weight map')


    args = parser.parse_args()

    compute_psf(args.path, args.telescope, useweight=args.useweight)

