#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import numpy as np
from astropy.io import fits
from lacosmic import lacosmic
from gmadet.utils import cp_p


# Can try to automatise the contrast value with the estimated PSF FWHM
# Not used at the moment.

def run_lacosmic(filename, FWHM, contrast=5, cr_threshold=5,
                 neighbor_threshold=5.0, maxiter=4, outLevel=1):
    """Run lacosmic to remove cosmic rays from the input image"""

    imagelist = np.atleast_1d(filename)

    for i, ima in enumerate(imagelist):
        path, filename_ext = os.path.split(ima)
        if path:
            folder = path + "/"
        else:
            folder = ""

        filename2 = os.path.splitext(filename_ext)[0]

        #  Make copy of original image
        if outLevel == 2:
            cp_p(ima, folder + filename2 + "_CR_notcleaned.fits")

        hdulist = fits.open(ima)
        hdr = hdulist[0].header
        try:
            gain = hdr["GAIN"]
        except BaseException:
            gain = 1
        try:
            RN = hdr["RN"]
        except BaseException:
            RN = 10

        """
        if FWHM[i] > 2:
            contrast = 2
        else:
            contrast = contrast
        """
        data = np.asarray(hdulist[0].data, dtype=float)
        lacosmic_res = lacosmic(
            data, contrast, cr_threshold, neighbor_threshold,
            effective_gain=gain, readnoise=RN, maxiter=maxiter
        )

        #  Create image cleaned from cosmic rays
        hdulist[0].data = lacosmic_res[0]
        hdulist.writeto(ima, overwrite=True)

        if outLevel == 2:
            #  Create mask of cosmic rays
            hdulist[0].data = np.asarray(lacosmic_res[1], dtype=int)
            hdulist.writeto(
                folder +
                filename2 +
                "_CRmask.fits",
                overwrite=True)
