#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run psfex to compute PSF."""

import errno
import glob
import os
import shutil
import subprocess
import sys
import numpy as np
from astropy.io import fits
from gmadet.utils import (rm_p, mv_p, mkdir_p, getpath)
import xmltodict


def psfex(filename, config, useweight=False,
          verbose="NORMAL", outLevel=0, outDir=''):
    """Compute PSF in astronomical images"""

    FWHM_list = []

    # imagelist=glob.glob(path+'/*.fits')
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
        print("\nRunning psfex to estimate FWHM in %s" % ima)
        root = os.path.splitext(ima)[0]
        if useweight:
            weight = root + ".weight.fits"
            subprocess.call(
                [
                    "sex",
                    ima,
                    "-c", config["psfex"]["sextractor"],
                    "-WEIGHT_IMAGE", weight,
                    "-VERBOSE_TYPE", verbose,
                    "-PARAMETERS_NAME", config["psfex"]["param"],
                    "-FILTER_NAME", config['sextractor']['convFilter'],
                ]
            )
        else:
            subprocess.call(
                [
                    "sex",
                    ima,
                    "-c", config["psfex"]["sextractor"],
                    "-VERBOSE_TYPE", verbose, 
                    "-PARAMETERS_NAME", config["psfex"]["param"],
                    "-FILTER_NAME", config['sextractor']['convFilter'],
                ]
            )
        cat = "preppsfex.cat"
        subprocess.call(
            ["psfex", cat, "-c", config["psfex"]["conf"],
             "-VERBOSE_TYPE", verbose]
        )
        rm_p(cat)

        #  Delete files depending on the required level of output files
        mv_p("snap_preppsfex.fits", root + "_psf.fits")

        # Get the mean PSF FWHM in pixels
        with open("psfex.xml") as fd:
            doc = xmltodict.parse(fd.read())
            FWHM_stats = doc["VOTABLE"]["RESOURCE"]["RESOURCE"]["TABLE"][0]["DATA"][
                "TABLEDATA"
            ]["TR"]["TD"][20:23]
            FHWM_min = float(FWHM_stats[0])
            FHWM_mean = float(FWHM_stats[1])
            FHWM_max = float(FWHM_stats[2])

            print("\nFWHM min: %.2f pixels" % FHWM_min)
            print("FWHM mean: %.2f pixels" % FHWM_mean)
            print("FWHM max: %.2f pixels\n" % FHWM_max)

        #  Get number of psf snapshot per axis
        psf_snaps = os.popen(
            "sed -n '/PSFVAR_NSNAP/p' %s" % config["psfex"]["conf"]
        ).read()
        nb_snaps = psf_snaps.split()[1]
        # Add info to the header
        hdulist = fits.open(root + "_psf.fits")
        hdr = hdulist[0].header
        hdr["FWHMMIN"] = str(FHWM_min)
        hdr["FWHMMEA"] = str(FHWM_mean)
        hdr["FWHMMAX"] = str(FHWM_max)
        hdr["PSF_NB"] = str(nb_snaps)
        hdulist.writeto(root + "_psf.fits", overwrite=True)

        FWHM_list.append(FHWM_mean)

        mv_p("preppsfex.psf", root + ".psf")
        rm_p("preppsfex.psf")
        rm_p("psfex.xml")
    return FWHM_list
