#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import numpy as np
from astropy.io import fits
from lacosmic import lacosmic
from astroscrappy import detect_cosmics
from gmadet.utils import cp_p


# Can try to automatise the contrast value with the estimated PSF FWHM
# Not used at the moment.

def run_lacosmic(filename, FWHM, contrast=5, cr_threshold=5,
                 neighbor_threshold=5.0, niter=4, outLevel=1):
    """Run lacosmic to remove cosmic rays from the input image"""

    imagelist = np.atleast_1d(filename)

    for i, ima in enumerate(imagelist):
        path, filename_ext = os.path.split(ima)
        if path:
            folder = path + "/"
        else:
            folder = ""

        filename2 = os.path.splitext(filename_ext)[0]

        # Make copy of original image
        if outLevel == 2:
            cp_p(ima, folder + filename2 + "_CR_notcleaned.fits")

        hdulist = fits.open(ima)
        hdr = hdulist[0].header
        gain = hdr.get('GAIN', 1)
        RN = hdr.get('RN',10)

        """
        if FWHM[i] > 2:
            contrast = 2
        else:
            contrast = contrast
        """
        data = np.asarray(hdulist[0].data, dtype=float)
        lacosmic_res = lacosmic(
            data, contrast, cr_threshold, neighbor_threshold,
            effective_gain=gain, readnoise=RN, maxiter=niter
        )

        # Create image cleaned from cosmic rays
        hdulist[0].data = lacosmic_res[0]
        hdulist.writeto(ima, overwrite=True)

        if outLevel == 2:
            # Create mask of cosmic rays
            hdulist[0].data = np.asarray(lacosmic_res[1], dtype=np.int16)
            hdulist.writeto(
                folder +
                filename2 +
                "_CRmask.fits",
                overwrite=True)


def run_astroscrappy(
        filename,
        psffwhm,
        contrast=5,
        cr_threshold=4.5,
        niter=4,
        inmask=None,
        sigfrac=0.3,
        readnoise=10,
        satlevel=50000,
        pssl=0.0,
        sepmed=True,
        cleantype='medmask',
        fsmode='median',
        psfmodel='gauss',
        psfsize=7,
        psfk=None,
        psfbeta=4.765,
        verbose=False,
        outLevel=1
        ):
    """Run astroscrappy to remove cosmics"""

    imagelist = np.atleast_1d(filename)
    psffwhm = np.atleast_1d(psffwhm)

    for i, ima in enumerate(imagelist):
        path, filename_ext = os.path.split(ima)
        if path:
            folder = path + "/"
        else:
            folder = ""

        filename2 = os.path.splitext(filename_ext)[0]

        # Make copy of original image
        if outLevel == 2:
            cp_p(ima, folder + filename2 + "_CR_notcleaned.fits")

        hdulist = fits.open(ima)
        hdr = hdulist[0].header

        saturate = np.min([satlevel, hdr.get('SATURATE', 50000)])
        gain = hdr.get('GAIN', 1)
        readnoise = hdr.get('RN', readnoise)

        crmask, cleanarr = detect_cosmics(
                hdulist[0].data,
                inmask=inmask,
                sigclip=cr_threshold,
                sigfrac=sigfrac,
                objlim=contrast,
                gain=gain,
                readnoise=readnoise,
                satlevel=saturate,
                pssl=pssl,
                niter=niter,
                sepmed=sepmed,
                cleantype=cleantype,
                fsmode=fsmode,
                psfmodel=psfmodel,
                psffwhm=psffwhm[i],
                psfsize=psfsize,
                psfk=psfk,
                psfbeta=psfbeta,
                verbose=verbose)

        # Create image cleaned from cosmic rays
        hdulist[0].data = cleanarr
        hdulist.writeto(ima, overwrite=True)

        if outLevel == 2:
            # Create mask of cosmic rays
            hdulist[0].data = np.asarray(crmask, dtype=np.int16)
            hdulist.writeto(
                folder +
                filename2 +
                "_CRmask.fits",
                overwrite=True)
