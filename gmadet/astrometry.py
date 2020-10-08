#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Scripts to perfrom astrometrci calibration.
Can use scamp or astrometry.net.

For astrometry.net You need to install locally astrometry.net and to
download the correct indexes files if needed.

"""

import subprocess
import sys
import os
import numpy as np
from astropy.io import fits
from gmadet.utils import (mv_p, mkdir_p, cp_p, getpath)
import xmltodict


def clean_tmp_files(filename, soft="scamp"):
    """Clean tempory files created by either astrometry.net or scamp"""

    fileroot = os.path.splitext(filename)[0]
    if soft == "astrometrynet":
        os.remove(fileroot[0] + "-indx.xyls")
        os.remove(fileroot[0] + ".axy")
        os.remove(fileroot[0] + ".corr")
        os.remove(fileroot[0] + ".match")
        os.remove(fileroot[0] + ".rdls")
        os.remove(fileroot[0] + ".solved")
        os.remove(fileroot[0] + ".wcs")
        os.remove(fileroot[0] + ".fits")
        os.rename(fileroot[0] + ".new", fileroot[0] + ".fits")
    elif soft == "scamp":
        os.remove("prepscamp.cat")
        os.remove("prepscamp.head")
        os.remove("scamp.xml")


def remove_astro_keywords(header):
    """ Remove from header all keywords related to astrometric solution  """

    #  First remove old keywords related to astrometric calibration
    keywords_to_remove = [
        "CRPIX1",
        "CRPIX2",
        "CRVAL1",
        "CRVAL2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
        "CDELT1",
        "CDELT2",
        "RADESYS",
        "CTYPE1",
        "CTYPE2",
        "EQUINOX",
        "CROTA1",
        "CROTA2",
    ]

    #  List of the first letters for standard astrometric coefficients.
    #  Such as TR, PV, SIA
    # coeff_astro = ['TR', 'SIA', 'A_', 'B_', 'AP_', 'BP_', 'LT', 'PV', 'PC']
    #  For PS1 skycell, though not sure what it means
    coeff_astro = ["PV", "PC"]

    for key, value in header.items():
        for coeff in coeff_astro:
            _len = len(coeff)
            if key[:_len] in [coeff]:
                keywords_to_remove.append(key)

    for keyword in keywords_to_remove:
        if keyword in header:
            del header[keyword]

    return header


def header_from_string(scamphead):
    """Create fits header from scamp .head output file"""

    s = ""
    with open(scamphead) as f:
        for line in f:
            s = s + line + "\n"
    header = fits.Header.fromstring(s, sep="\n")

    return header


def update_headers_scamp(filename, scamphead, pixelscale):
    """Modify the header after running scamp"""

    hdulist = fits.open(filename)
    # Verify and try to fix issue with fits standard
    hdulist.verify("fix")
    hdr = hdulist[0].header
    #  First remove old keywords related to astrometric calibration
    hdr = remove_astro_keywords(hdr)

    newheader = header_from_string(scamphead)

    for key, value in newheader.items():
        # print (key, value)
        # truncate long keys
        if len(key) > 8:
            key = key[:7]
        try:
            hdr.set(key.upper(), value)
        except BaseException:
            try:
                hdr.set(key.upper(), str(value))
            except BaseException:
                pass
    hdr.set("CDELT1", pixelscale[0])
    hdr.set("CDELT2", pixelscale[1])

    if (hdr["CTYPE1"] != "RA---TPV") and (hdr["CTYPE2"] != "DEC--TPV"):
        print("\nWARNING: scamp did not set CTYPE1, CTYPE2 to RA---TPV and DEC--TPV.")
        print(
            "Set them to to these values by hand, otherwise it can not be read by astropy.wcs"
        )
        print(
            "Likely due to some SIP distortions parameters already present in headers."
        )
        print("One might check the astrometry to be safe.\n")
        hdr["CTYPE1"] = "RA---TPV"
        hdr["CTYPE2"] = "DEC--TPV"

    hdulist.writeto(filename, overwrite=True)


def astrometrynet(filename, radius=4, scaleLow=3.5, scaleHigh=4,
                  scaleUnits="arcsecperpix"):
    """Run astrometry.net on the input image"""

    # Get the RA and DEC before removing the keywords as required
    # by astrometry.net.
    # Is it really required? need to check
    hdul = fits.open(filename)
    hdr = hdul[0].header

    ra = str(hdr["CRVAL1"])
    dec = str(hdr["CRVAL2"])

    del hdr["CRVAL1"]
    del hdr["CRVAL2"]
    del hdr["CRPIX1"]
    del hdr["CRPIX2"]
    hdul[0].header = hdr
    hdul.writeto(filename, overwrite=True)
    hdul.close()

    #  Run astrometry.net
    subprocess.call(
        [
            "solve-field", filename,
            "--ra", ra,
            "--dec", dec,
            "--radius", str(radius),
            "--scale-units", str(scaleUnits),
            "--scale-low", str(scaleLow),
            "--scale-high", str(scaleHigh),
            "--crpix-center", "--tweak-order", "5",
            "--uniformize", "0",
            "--no-plots",
            "--overwrite",
        ]
    )

    #  Delete temporary files
    clean_tmp_files(filename, soft="astrometrynet")


def scamp(filename, config, accuracy=0.5, itermax=10,
          band=None, useweight=False, CheckPlot=False,
          verbose="NORMAL"):
    """Compute astrometric solution of astronomical image using scamp"""
    
    path = os.path.dirname(filename)
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
        print("Performing astrometric calibration on %s using SCAMP." % ima)
        # print ('You required astrometric precision of %.3f arcsec.' % accuracy)

        root = os.path.splitext(ima)[0]
        _name = root.split("/")[-1]

        if CheckPlot:
            plot_fmt = "PNG"
            plotnames = (
                "%s_fgroups,%s_distort,%s_astr_interror2d,%s_astr_interror1d,%s_astr_referror2d,%s_astr_referror1d,%s_astr_chi2,%s_psphot_error"
                % (root, root, root, root, root, root, root, root)
            )
            plottypes = "FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR"
        else:
            plot_fmt = "NULL"
            plotnames = " "
            plottypes = " "

        if config["telescope"] == "PS1" and band is not None:
            astref_band = band
        else:
            astref_band = "DEFAULT"

        #  Initialise the while loop
        i = 0
        #  Dummy offset
        mean_offset = 100
        while (mean_offset >= accuracy) and (i <= itermax):
            i += 1
            #  Create catalog using sextractor
            # print ('Create FITS-LDAC file from SExtractor')
            subprocess.call(
                [
                    "sex",
                    "-c", config["scamp"]["sextractor"],
                    "-PARAMETERS_NAME", config["scamp"]["param"],
                    "-VERBOSE_TYPE", verbose,
                    "-FILTER_NAME", config['sextractor']['convFilter'],
                    ima,
                ]
            )
            #  Run SCAMP
            subprocess.call(
                [
                    "scamp",
                    "prepscamp.cat",
                    "-c", config["scamp"]["conf"],
                    "-ASTREF_BAND", astref_band,
                    "-CHECKPLOT_DEV", plot_fmt,
                    "-CHECKPLOT_NAME", plotnames,
                    "-CHECKPLOT_TYPE", plottypes,
                    "-VERBOSE_TYPE", verbose,
                ]
            )
            #  Check astrometry offset
            with open("scamp.xml") as fd:
                doc = xmltodict.parse(fd.read())
            """
            offset = doc['VOTABLE']['RESOURCE']['RESOURCE']['TABLE'][0]['DATA']['TABLEDATA']['TR']['TD'][34]
            offset = offset.split(' ')
            offset_axis1 = float(offset[0])
            offset_axis2 = float(offset[1])
            """
            header = header_from_string("prepscamp.head")
            offset_axis1 = float(header["ASTRRMS1"] * 3600)
            offset_axis2 = float(header["ASTRRMS2"] * 3600)

            mean_offset = np.mean([offset_axis1, offset_axis2])
            print(
                "Astrometric precision after run %d: %.2f arcseconds. Required: %.2f."
                % (i, mean_offset, accuracy)
            )

            pixelscale = doc["VOTABLE"]["RESOURCE"]["RESOURCE"]["TABLE"][0]["DATA"][
                "TABLEDATA"
            ]["TR"]["TD"][18].split("  ")
            pixelscale = [float(pixelscale[0]) / 3600,
                          float(pixelscale[1]) / 3600]
            # Update header of input fits file
            update_headers_scamp(ima, "prepscamp.head", pixelscale)
        # cp_p('prepscamp.cat', _name.split('.')[0]+'.cat')
        #  Delete temporary files
        clean_tmp_files(ima, soft="scamp")
    print("\n")


def astrometric_calib(filenames, config, soft="scamp",
                      accuracy=0.5, itermax=10, verbose="NORMAL"):
    """perform astrometric calibration"""

    imagelist = np.atleast_1d(filenames)
    for ima in imagelist:
        #  Use scamp for astrometric calibration
        if soft == "scamp":
            scamp(ima, config, accuracy=accuracy,
                  itermax=itermax, verbose=verbose)

        #  Use astrometry.net for astrometric calibration
        elif soft == "astrometrynet":
            # Get pixel scale in degrees
            header = fits.getheader(ima)
            try:
                pixScale = abs(header["CDELT1"])
            except Exception:
                try:
                    pixScale = abs(header["CD1_1"])
                except Exception:
                    print(
                        "Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1"
                    )
            #  Set up boundaries for plate scale for astrometry.net
            scaleLow = 0.7 * pixScale * 3600
            scaleHigh = 1.3 * pixScale * 3600
            radius = max(header["NAXIS1"] * pixScale,
                         header["NAXIS2"] * pixScale)
            asrometrynet(ima, radius=radius, scaleLow=scaleLow,
                         scaleHigh=scaleHigh)
