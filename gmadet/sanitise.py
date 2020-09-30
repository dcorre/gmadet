#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import numpy as np
from astropy.io import fits


def sanitise_headers(filename):
    """
    Keep only important keywords.
    Helps not to crash when headers are not well formated.
    """

    #  Need to keep minimum information about astrometric
    #  calibration for initialising scamp
    #  List of official keywords:
    #  https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
    # Other keywords commonly used:
    #  https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html
    #  Keywords not present in this list are discarded.
    keywords_to_keep = [
        "SIMPLE",
        "BITPIX",
        "NAXIS",
        "NAXIS1",
        "NAXIS2",
        "EXTEND",
        "DATE-OBS",
        "TELESCOP",
        "INSTRUME",
        "OBJECT",
        "EXPTIME",
        "FILTER",
        "GAIN",
        "SATURATE",
        "EQUINOX",
        "EPOCH",
        "RADESYS",
        "CTYPE1",
        "CTYPE2",
        "CUNIT1",
        "CUNIT2",
        "CRVAL1",
        "CRVAL2",
        "CRPIX1",
        "CRPIX2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
        "CDELT1",
        "CDELT2",
        "CROTA1",
        "CROTA2",
        "BSCALE",
        "BZERO",
        "BUNIT",
        "AIRMASS",
        "END",
    ]
    hdulist = fits.open(filename)
    # Verify and try to fix issue with fits standard
    hdulist.verify("fix")
    hdr = hdulist[0].header

    keywords2delete = []
    for key, value in hdr.items():
        #  if EXPOSURE keyword present, rename it EXPTIME
        if key == "EXPOSURE":
            #  check if EXPTIME exists, otherwise create it
            try:
                if hdr["EXPTIME"]:
                    pass
            except BaseException:
                hdr["EXPTIME"] = value
                #  if EXPOSURE keyword present, rename it EXPTIME
        if key == "FILTERS":
            #  check if FILTER exists, otherwise create it
            try:
                if hdr["FILTER"]:
                    pass
            except BaseException:
                hdr["FILTER"] = value

        if key not in keywords_to_keep:
            keywords2delete.append(key)
    keywords2delete = np.unique(keywords2delete)
    for key in keywords2delete:
        del hdr[key]
    hdulist.writeto(filename, overwrite=True)


def sanitise_data(filename):
    """
    Make sure that the data are in the primary hdu.
    Otherwise some astromatic softs are not working properly.
    Other solution would be to check astromatic soft config.
    """
    # Make sure to use only the Primary hdu.
    # sextractor, scamp seems to crash otherwise.
    hdul = fits.open(filename)
    # hdul.verify('fix')
    print(hdul.info())
    if len(hdul) > 1:
        print("More than one HDU in fits file.")
        print("Keeping only the PrimaryHDU.")
        newhdu = fits.PrimaryHDU()
        newhdu.data = hdul[0].data
        newhdu.header = hdul[0].header
        newhdulist = fits.HDUList([newhdu])
        newhdulist.writeto(filename, overwrite=True)
        hdul.close()


def sanitise_fits(filename):
    """Call function to sanitise fits headers and data"""
    sanitise_headers(filename)
    sanitise_data(filename)
