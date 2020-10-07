#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

"""
Contain scripts to extract sources in  astronomical images
using sextractor.

"""
import os
import subprocess
import numpy as np

from astropy.io import ascii, fits
from astropy import wcs
from astropy.table import Table

from gmadet.utils import mv_p, getpath


def run_sextractor(filelist, FWHM_list, thresh, telescope, config,
                   verbose="NORMAL", subFiles=None, outLevel=1):
    """Run sextractor """
    # if substraction have been performed
    # Run sextractor on input image to get photometry calibrated and
    # on substracted image to get interesting sources
    if subFiles:
        subFiles = np.array(subFiles)

        """
        # No mask on iput data
        mask = ["None"] * len(filelist)
        mask.extend([im for im in subFiles[:, 3]])
        weight_type = ["NONE"] * len(filelist)
        weight_type.extend(["MAP_WEIGHT"] * len(mask))

        psfs = [os.path.splitext(im)[0] + ".psf" for im in filelist]
        #  assume we take the psf of the original file
        psfs.extend([psfs[0]] * len(subFiles[:, 2]))
        # Rather take the original before registration as it introduces
        # artefact. Note: there is a crossmatch before calibration so no pb.
        filelist = [im for im in filelist]
        """
        # No mask on iput data
        mask = ['None'] * len(subFiles[:, 0])
        mask.extend([im for im in subFiles[:, 3]])
        weight_type = ['NONE'] * len(subFiles[:, 0])
        weight_type.extend(["MAP_WEIGHT"] * len(mask))
        psfs = [im.split(".")[0] + ".psf" for im in filelist]
        # Assume PSF is same for the input image cutouts
        # as it will take time to recompute it.
        psfs.extend([psfs[0]] * (len(subFiles[:, 0])-1))
        #  assume we take the psf of the original file for substracted im.
        psfs.extend([psfs[0]] * len(subFiles[:, 2]))

        filelist = [im for im in subFiles[:, 0]]
        filelist.extend([im for im in subFiles[:, 2]])
        # Duplicate FWHM list
        # FWHM_list.extend(FWHM_list)
        FWHM_list = np.ravel([[i] * len(filelist) for i in FWHM_list])
    else:
        mask = ["None"] * len(filelist)
        weight_type = ["NONE"] * len(filelist)

        psfs = [os.path.splitext(im)[0] + ".psf" for im in filelist]

    for i, filename in enumerate(filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + "/"
        else:
            folder = ""
        #  Get rid of the extension to keep only the name
        filename2 = os.path.splitext(filename_ext)[0]
        if outLevel == 2:
            checkimage_type = "BACKGROUND, SEGMENTATION"
            checkimage_name = (
                folder
                + filename2
                + "_background.fits"
                + ", "
                + folder
                + filename2
                + "_segmentation.fits"
            )

        else:
            checkimage_type = "NONE"
            checkimage_name = " "

        subprocess.call(
            [
                "sex",
                "-c", config["sextractor"]["conf"],
                filename,
                "-WEIGHT_TYPE", str(weight_type[i]),
                "-WEIGHT_IMAGE", str(mask[i]),
                "-SEEING_FWHM", str(FWHM_list[i]),
                "-DETECT_THRESH", str(thresh),
                "-PARAMETERS_NAME", config["sextractor"]["param"],
                # '-FILTER_NAME', config['sextractor']['default_conv'],
                "-CHECKIMAGE_TYPE", checkimage_type,
                "-CHECKIMAGE_NAME", checkimage_name,
                "-VERBOSE_TYPE", verbose,
                "-PSF_NAME", psfs[i],
                "-CATALOG_NAME", folder + filename2 + "_SourcesDet.cat",
                "-FILTER_NAME", config["sextractor"]["convFilter"]
            ]
        )


def filter_sources(filelist, soft, edge_cut=32, sigma=1, subFiles=None):
    # if substraction have been performed
    if subFiles:
        subfiles = np.array(subFiles)
        filelist = [im for im in subfiles[:, 0]]
        filelist.extend([im for im in subfiles[:, 2]])
        """
        # Rather take the original before registration as it introduces
        # artefact
        original_filelist = list(filelist)
        original_filelist.extend([im for im in filelist] * len(subFiles))
        # originallist = [im for im in filelist]
        # filelist = []
        # for im, sub in zip(originallist, subFiles[:,2]):
        #    filelist.append(im)
        #    filelist.append(sub)
        subfiles = np.array(subFiles)
        filelist = [im for im in filelist]
        filelist.extend([im for im in subfiles[:, 2]])
        """
    for filename in filelist:
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + "/"
        else:
            folder = ""
        #  Get rid of the extension to keep only the name
        filename2,fileext = os.path.splitext(filename_ext)

        if soft == "sextractor":
            sources = ascii.read(
                folder + filename2 + "_SourcesDet.cat", format="sextractor"
            )
            mv_p(
                folder + filename2 + "_SourcesDet.cat",
                folder + filename2 + "_SourcesDetnoFilter.cat",
            )
            #  only if there is at least one detection
            if sources:
                #  Remove sources too close to the imge edges
                header = fits.getheader(folder + filename2 + fileext)
                imsize = [int(header["NAXIS1"]), int(header["NAXIS2"])]
                mask_edge = (
                    (sources["X_IMAGE"] > edge_cut)
                    & (sources["Y_IMAGE"] > edge_cut)
                    & (sources["X_IMAGE"] < imsize[0] - edge_cut)
                    & (sources["Y_IMAGE"] < imsize[1] - edge_cut)
                )

                """
                # Remove sources that are likely cosmic rays
                # Compute flux ratio as total_flux / nb_pixels for each source
                # Compute the median and std for sources with:
                # 10 < nb_pixels < 100
                # Higher than 10 to discard cosmics
                # Smaller than 100 to discard saturated sources
                # Remove sources with:
                # nb_pixels < 10 and fluxratio > fluxratio_med + fluxratio_std
                flux = sources['FLUX_AUTO']
                nbpix = sources['ISOAREA_IMAGE']
                fluxratio = flux / nbpix
                # Compute the median and std on the original image only
                # This assumes that the substracted image follows the original
                # image in the for loop
                if '_sub' not in filename2:
                    mask = (nbpix > 10) & (nbpix < 100)
                    fluxratio_med = np.median(fluxratio[mask])
                    fluxratio_std = np.std(fluxratio[mask])

                mask_cosmics = (nbpix < 10) & (
                    fluxratio > fluxratio_med + sigma * fluxratio_std)
                mask_tot = np.bitwise_and(mask_edge, np.invert(mask_cosmics))
                """

                # Remove sources too close to the edges
                # sources_filt = sources[mask_edge]
                # Flag sources too close to the edges
                #  Only if there is at least one detection
                edge_flag = np.array(["N"] * len(sources))
                edge_flag[~mask_edge] = "Y"
                sources["edge"] = edge_flag
            else:
                #  Add something to make script not crashing
                sources["edge"] = []
            sources.write(
                folder + filename2 + "_SourcesDet.cat",
                format="ascii.commented_header",
                overwrite=True,
            )


def convert_xy_radec(filelist, soft="sextractor", subFiles=None):
    """
    Performs pyraf transformation of x-y into RA-DEC coordinates
    filename is WITHOUT suffix .fits
    Input is the *.magfiltered file from select_good_stars() function
    Output is the *.magwcs file
    """
    # If substraction has been performed
    if subFiles:
        subfiles = np.array(subFiles)
        original_filelist = [im for im in subfiles[:, 0]]
        original_filelist.extend([im for im in original_filelist] * 2)

        reference_filelist = ["None"] * len(subfiles[:, 0])
        reference_filelist.extend([im for im in subfiles[:, 1]])

        filelist = [im for im in subfiles[:, 0]]
        filelist.extend([im for im in subfiles[:, 2]])
        
        """
        original_filelist = list(filelist)
        original_filelist.extend([im for im in filelist] * len(subFiles))
        subfiles = np.array(subFiles)
        reference_filelist = ["None"] * len(filelist)
        reference_filelist.extend([im for im in subfiles[:, 1]])
        # filelist = [im for im in subfiles[:, 0]]
        # Rather take the original before registration as it introduces
        # artefact
        filelist = [im for im in filelist]
        filelist.extend([im for im in subfiles[:, 2]])
        """
    else:
        original_filelist = filelist
        reference_filelist = ["None"] * len(filelist)

    for filename, original_filename, refimage in zip(
        filelist, original_filelist, reference_filelist
    ):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + "/"
        else:
            folder = ""

        #  Get rid of the extension to keep only the name
        filename2 = os.path.splitext(filename_ext)[0]

        magfilewcs = folder + filename2 + ".magwcs"

        if soft == "sextractor":
            sources = ascii.read(
                folder + filename2 + "_SourcesDet.cat",
                format="commented_header"
            )
            #  If there is at least one detection
            if sources:
                header = fits.getheader(filename)
                w = wcs.WCS(header)
                ra, dec = w.wcs_pix2world(
                    sources["X_IMAGE"], sources["Y_IMAGE"], 1)

                filenames = [filename] * len(ra)
            else:
                ra = []
                dec = []
                filenames = []
            data = Table(
                [
                    sources["X_IMAGE"],
                    sources["Y_IMAGE"],
                    ra,
                    dec,
                    sources["MAG_AUTO"],
                    sources["MAGERR_AUTO"],
                    sources["edge"],
                    sources["CHI2_PSF"],
                    sources["MAG_PSF"],
                    sources["MAGERR_PSF"],
                    sources["FWHM_IMAGE"],
                    sources["FWHMPSF_IMAGE"],
                    filenames,
                ],
                names=[
                    "Xpos",
                    "Ypos",
                    "RA",
                    "DEC",
                    "Mag_inst",
                    "Magerr_inst",
                    "edge",
                    "psf_chi2",
                    "psf_mag",
                    "psf_magerr",
                    "FWHM",
                    "FWHMPSF",
                    "filenames",
                ],
            )

        # Flag to identify substraction image
        if "_sub" in magfilewcs:
            data["FlagSub"] = ["Y"] * len(data)
            data["OriginalIma"] = [original_filename] * len(data)
            data["RefIma"] = [refimage] * len(data)
        else:
            data["FlagSub"] = ["N"] * len(data)
            data["OriginalIma"] = [original_filename] * len(data)
            data["RefIma"] = [refimage] * len(data)

        data.write(magfilewcs, format="ascii.commented_header", overwrite=True)
        check_RADEC = data["RA", "DEC"]
        check_RADEC.write(
            magfilewcs + "2", format="ascii.commented_header", overwrite=True
        )
