#!usr/bin/env python
# -*- coding: utf-8 -*-

# Python module for performing photometric calibration
# David Corre, corre@lal.in2p3.fr
#

import os
import numpy as np
import matplotlib.pyplot as plt
import warnings

from gmadet.crossmatch import run_xmatch
from gmadet.utils import get_phot_cat, filter_catalog_data
from gmadet.phot_conversion import *

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import vstack, Table, Column, join
from astropy.stats import sigma_clip

from copy import deepcopy

warnings.simplefilter(action="ignore")


def crossmatch(detected_sources, radius, pixScale, catalog, nb_threads=4):
    """
    Crossmatch the output file of Sextractor / pyRAF
    Remove duplicate source during crossmatch
    Keep only the closest match
    """

    # Add units
    detected_sources["_RAJ2000"] *= u.deg
    detected_sources["_DEJ2000"] *= u.deg
    # Add index for each source
    detected_sources["idx"] = np.arange(len(detected_sources))

    # Do not consider sources close to the image edges
    # Do not add sources from substracted image
    mask = (detected_sources["edge"] == "N") & \
            (detected_sources["FlagSub"] == "N")

    # Run xmatch to crossmatch detected sources with available catalog
    # reduce size of the input catalog to its minimum
    cat = deepcopy(detected_sources["_RAJ2000", "_DEJ2000", "idx"])

    crossmatch = run_xmatch(
        cat[mask],
        catalog,
        radius * pixScale * 3600,
        nb_threads)
    # Do not consider duplicates
    # Consider only closest crossmatch if multiple association.
    # Assume that the first occurence is the closest one.
    _, referenced_star_idx = np.unique(crossmatch["idx"],
                                       return_index=True)
    ref_sources = crossmatch[referenced_star_idx]
    #keep_idx = np.isin(detected_sources['idx'], ref_sources['idx'])
    ref_sources = join(ref_sources, detected_sources, join_type='left')
    return ref_sources


def conv_mag_sys(data, band, catalog):
    """
    For a gieven telescope magnitude and refence star catalog, perform
    conversion from catalog bands to telescope bands
    """
    if band in ["g", "r", "i", "z", "y", "g+r"]:
        if catalog == "II/349/ps1":
            # No transformation
            bands = band.split("+")
            if len(bands) > 1:
                # jansky = []
                jansky = np.zeros(len(data))
                mag_err = np.zeros(len(data))
                for filt in bands:
                    jansky = jansky + 3631 * \
                        10 ** (-0.4 * (data["%smag" % filt]))
                    # Add error in quadrature. Might be too simple
                    mag_err = mag_err + data["e_%smag" % filt] ** 2
                newmag = -2.5 * np.log10(jansky / (3631))
                newmag_err = np.sqrt(mag_err)
            else:
                newmag = data["%smag" % bands[0]]
                newmag_err = data["e_%smag" % bands[0]]

            data["mag_cat"] = newmag
            data["magerr_cat"] = newmag_err
            data["magsys"] = "AB"
            catalogName = "Pan-STARRS DR1"

        elif catalog == "V/147/sdss12":
            # No transformation
            catalogName = "SDSS DR12"
            pass
        elif catalog == "I/345/gaia2":
            bands = band.split("+")
            if len(bands) > 1:
                # jansky = []
                jansky = np.zeros(len(data))
                mag_err = np.zeros(len(data))
                for filt in bands:
                    data2 = gaia2SDSS(filt, data)
                    jansky = jansky + 3631 * \
                        10 ** (-0.4 * (data2["%s_SDSSMag" % filt]))
                    # Add error in quadrature. Might be too simple
                    mag_err = mag_err + data2["calib_err"] ** 2
                newmag = -2.5 * np.log10(jansky / (3631))
                newmag_err = np.sqrt(mag_err)
                data["mag_cat"] = newmag
                data["magerr_cat"] = newmag_err
                data["magsys"] = "AB"
            else:
                data = gaia2SDSS(band, data)
                data.rename_column("%s_SDSSMag" % band, "mag_cat")
                data.rename_column("calib_err", "magerr_cat")
                data["magsys"] = "AB"
            catalogName = "Gaia DR2"

    elif band in ["B", "V", "R", "I"]:
        if catalog == "II/349/ps1":
            data = PS2Johnson(band, data)
            data['mag_cat'] = data["%sMag" % band]
            data['magerr_cat'] = data["calib_err"]
            #data.rename_column("%sMag" % band, "mag_cat")
            #data.rename_column("calib_err", "magerr_cat")
            data["magsys"] = "AB"
            catalogName = "Pan-STARRS DR1"
        elif catalog == "V/147/sdss12":
            data = SDSS2Johnson(band, data)
            data['mag_cat'] = data["%sMag" % band]
            data['magerr_cat'] = data["calib_err"]
            #data.rename_column("%sMag" % band, "mag_cat")
            #data.rename_column("calib_err", "magerr_cat")
            data["magsys"] = "AB"
            catalogName = "SDSS DR12"
        elif catalog == "I/345/gaia2":
            data = gaia2Johnson(band, data)
            data['mag_cat'] = data["%sMag" % band]
            data['magerr_cat'] = data["calib_err"]
            #data.rename_column("%sMag" % band, "mag_cat")
            #data.rename_column("calib_err", "magerr_cat")
            data["magsys"] = "Vega"
            catalogName = "Gaia DR2"
        elif catalog == "I/284/out":
            # Need to add
            catalogName = "USNO-B1"
            pass

    return data, catalogName


def zeropoint(data, sigma, quadrant, folder,
              fname, band, catalog, doPlot=False):
    """"Compute zeropoints"""
    # data.show_in_browser()
    # Sigma clipping for zeropoints
    # clip = sigma_clip(a['Delta_Mag'], sigma = 1.5, masked=True)
    # clip_mask = np.invert(clip.recordmask)
    # a = a[clip_mask]
    # Remove objects to get all objects with delte_mag < 1 sigma
    delta_mag = data["mag_inst"] - data["mag_cat"]
    clip = sigma_clip(delta_mag, sigma=sigma)
    clip_mask = np.invert(clip.recordmask)

    newdata = data[clip_mask]
    delta_mag = newdata["mag_inst"] - newdata["mag_cat"]
    delta_mag_median = np.median(delta_mag)
    delta_mag_std = np.std(delta_mag)

    newdata.write(
        folder + fname + "_ZP_%d.dat" % quadrant,
        format="ascii.commented_header",
        overwrite=True,
    )

    if doPlot:
        # newdata.show_in_browser(jsviewer=True)
        plt.figure()
        plt.scatter(data["mag_inst"], data["mag_cat"], color="blue")
        median_nosigmaclip = np.median(data["mag_inst"] - data["mag_cat"])
        std_nosigmaclip = np.std(data["mag_inst"] - data["mag_cat"])
        plt.plot(
            data["mag_inst"],
            data["mag_inst"] -
            median_nosigmaclip,
            color="green")
        plt.xlabel("Instrumental magnitude")
        plt.ylabel("%s band %s" % (band, catalog))
        plt.title(
            "zeropoints without sigma clipping.\nMedian: %.2f, std: %.2f"
            % (median_nosigmaclip, std_nosigmaclip)
        )
        plt.savefig(folder + fname + "_ZP_%d_nosigmaClipping.png" % quadrant)

        plt.figure()
        plt.scatter(newdata["mag_inst"], newdata["mag_cat"], color="blue")
        plt.plot(
            newdata["mag_inst"], newdata["mag_inst"] - delta_mag_median,
            color="green"
        )
        plt.xlabel("Instrumental magnitude")
        plt.ylabel("%s band %s" % (band, catalog))
        plt.title(
            "zeropoints after sigma clipping (sigma=1.5).\n"
            "Median: %.2f, std: %.2f"
            % (-delta_mag_median, delta_mag_std)
        )
        plt.savefig(folder + fname + "_ZP_%d.png" % quadrant)
        # plt.show()

    return newdata, delta_mag_median, delta_mag_std


def phot_calib(detected_sources, telescope, radius=3,
               doPlot=True, subFiles=None):
    """Perform photometric calibration using catalogs"""

    # Assume that the pixelscale is the same for all images.
    # It is the case as all the images are part of the same.
    # If the script is used outside of the gamdet_run, this need
    # to be changed.
    # Benefit is to make only one crossmatch with all sources instead
    # of multiple crossmatch with small number of sources.
    # Get pixel scale in degrees
    header = fits.getheader(detected_sources['OriginalIma'][0])
    try:
        pixScale = abs(header["CDELT1"])
    except Exception:
        try:
            pixScale = abs(header["CD1_1"])
        except Exception:
            print(
                "Pixel scale could not be found in fits header.\n"
                "Expected keyword: CDELT1 or CD1_1"
            )
    # Get filter and catalog to perform photometric calibration
    band_DB, band_cat, catalog = get_phot_cat(detected_sources['OriginalIma'][0],
            telescope)

    print("Crossmatching with catalog.")
    # Crossmtach sources in the non substracted images.
    # Filtering sources in substracted images is done inside.
    ref_sources = crossmatch(detected_sources, radius, pixScale, catalog)

    # Remove extended sources and bad measurements from reference
    # stars catalog
    # ref_sources.show_in_browser(jsviewer=True)
    print("Removed extended objects or with band quality flags.")
    good_ref_sources = filter_catalog_data(ref_sources, catalog)
    # good_ref_sources.show_in_browser(jsviewer=True)

    # create columns
    good_ref_sources['mag_cat'] = np.zeros(len(good_ref_sources))
    good_ref_sources['magerr_cat'] = np.zeros(len(good_ref_sources)) 
    good_ref_sources['magsys'] = ['None'] * len(good_ref_sources)
    

    deltaMagMedianlist = []
    deltaMagStdlist = []
    filename_list = []
    band_cat_list = []
    band_DB_list = []
    # Compute zeropoints
    for i, key in enumerate(
            detected_sources.group_by("OriginalIma").groups.keys):
        print("Processing photometric calibration for ", key[0])

        mask_key = good_ref_sources['OriginalIma'] == key[0]
        # Get path and filename to images
        path, fname_ext = os.path.split(key[0])
        if path:
            folder = path + "/"
        else:
            folder = ""

        #  Get rid of the extension to keep only the name
        fname2 = fname_ext.split(".")[0]
        extension = ""
        for ext in fname_ext.split(".")[1:]:
            extension = extension + "." + ext

        # Get filter and catalog to perform photometric calibration
        band_DB, band_cat, catalog = get_phot_cat(key[0], telescope)

        # Transform filter bands in catalog to telescope ones
        print("Convert magnitude system.")
        ref_cat, catalogName = conv_mag_sys(
            good_ref_sources[mask_key], band_cat, catalog)
        print("Compute zeropoint.")
        ref_cat_calibrated, deltaMagMedian, deltaMagStd = zeropoint(
            ref_cat, 1.5, i, folder, fname2, band_cat, catalogName, doPlot=True
        )

        deltaMagMedianlist.append(deltaMagMedian)
        deltaMagStdlist.append(deltaMagStd)
        filename_list.append(key[0])
        band_cat_list.append(band_cat)
        band_DB_list.append(band_DB)
    
    print("Add magnitude to candidates.")
    # Apply photmetric calibration to candidates
    mag_calib_col = Column(np.zeros(len(detected_sources)),
                           name="mag_calib")
    mag_calib_err_col = Column(np.zeros(len(detected_sources)),
                               name="mag_calib_err")
    ZP_col = Column(np.zeros(len(detected_sources)),
                           name="ZP")
    ZP_err_col = Column(np.zeros(len(detected_sources)),
                           name="ZP_err")
    magsys_col = Column(["None"] * len(detected_sources),
                        name="magsys")
    filter_cat_col = Column(["None"] * len(detected_sources),
                            name="filter_cat")
    filter_DB_col = Column(["NoFilterFound"] * len(detected_sources),
                           name="filter_DB")

    detected_sources.add_columns([mag_calib_col, mag_calib_err_col,
                                 ZP_col, ZP_err_col, magsys_col, filter_cat_col,
                                 filter_DB_col])

    for j, filename in enumerate(filename_list):
        #  Compute calibrated magnitudes for transient candidates only
        mask = detected_sources["OriginalIma"] == filename
        detected_sources["mag_calib"][mask] = (
            detected_sources["mag_inst"][mask] - deltaMagMedianlist[j]
        )
        # Quadratic sum of statistics and calibration errors.
        detected_sources["mag_calib_err"][mask] = np.sqrt(
            detected_sources["mag_inst_err"][mask] ** 2 +
            deltaMagStdlist[j] ** 2
        )

        detected_sources["ZP"][mask] = deltaMagMedianlist[j] 
        detected_sources["ZP_err"][mask] = deltaMagStdlist[j] 

        # Define magnitude system
        if catalog == "II/349/ps1":
            detected_sources["magsys"][mask] = "AB"
        else:
            pass

        detected_sources["filter_cat"][mask] = band_cat_list[j] 
        detected_sources["filter_DB"][mask] = band_DB_list[j] 

    fname = detected_sources['OriginalIma'][0].split("_reg_")[0] + ".alldetections"
    detected_sources.write(
        fname,
        format="ascii.commented_header",
        overwrite=True,
    )

    return detected_sources
