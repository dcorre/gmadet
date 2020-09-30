#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.stats import SigmaClip
from photutils import (
        Background2D,
        SExtractorBackground,
        MMMBackground,
        ModeEstimatorBackground,
        MedianBackground,
        MeanBackground
)
from astropy.io import fits


def bkg_estimation(filename, box=(20, 20), filter_size=(3, 3),
                   bkg_estimator='SExtractor', sigma=3.,
                   sigma_lower=None, sigma_upper=None,
                   maxiters=10, outLevel=1):
    
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
        print("\nEstimate background in %s" % ima)
        root = os.path.splitext(ima)[0]

        hdulist = fits.open(ima)
        data = hdulist[0].data

        sigma_clip = SigmaClip(sigma=sigma,
                               sigma_lower=sigma_lower,
                               sigma_upper=sigma_upper,
                               maxiters=maxiters)
        
        if bkg_estimator == 'SExtractor':
            bkg_estimator = SExtractorBackground()
        elif bkg_estimator == 'MMM':
            bkg_estimator = MMMBackground()
        elif bkg_estimator == 'ModeEstimator':
            bkg_estimator = ModeEstimatorBackground()
        elif bkg_estimator == 'Median':
            bkg_estimator = MedianBackground()
        elif bkg_estimator == 'Mean':
            bkg_estimator = MeanBackground()

        bkg = Background2D(
            data,
            box,
            filter_size=filter_size,
            sigma_clip=sigma_clip,
            bkg_estimator=bkg_estimator,
        )

        #  Create image with background substracted
        hdulist[0].data = data - bkg.background
        hdulist.writeto(ima, overwrite=True)

        if outLevel == 2:
            #  Create image with 2D background
            bkg_file = root + "_bkg_map.fits"
            hdulist[0].data = bkg.background
            hdulist.writeto(bkg_file, overwrite=True)
