#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Scripts to perform crossmatch with astronomical catalogs
and with solar moving objects.

"""
import os
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astroquery import xmatch
from astroquery.imcce import Skybot
from astroML.crossmatch import crossmatch_angular
from copy import deepcopy
import multiprocessing as mp

def _run_xmatch(coordinates, catalog, radius):
    """
    Perform cross-match with a catalog using the CDS XMatch
    parameters: coordinates, catalog, radius:
                coordinates: astropy table with RA, DEC of all detected sources
                catalog: Vizier identifier of the catalog
                radius in arcsecond
    returns: astropy.table object

    Vizier catalog identifiers:
    Gaia DR2: I/345/gaia2
    SDSS DR12: V/147/sdss12
    2MASS: II/246/out
    USNO B1: I/284/out
    USNO A2: I/252/out
    GLADE 2: VII/281/glade2
    Panstarrs DR1: II/349/ps1
    """
    xmatch.XMatch.TIMEOUT = 3600

    matched_stars = xmatch.XMatch.query(
        coordinates,
        cat2="vizier:%s" % catalog,
        max_distance=radius * u.arcsec,
        colRA1="_RAJ2000",
        colDec1="_DEJ2000",
    )

    return matched_stars

def run_xmatch(coordinates, catalog, radius, nb_threads):
    """Run xmatch in parallel"""

    catalog_list = []
    idx_stop = []
    Ncat = len(coordinates)
    Ncut = int(Ncat / nb_threads)
    for i in range(nb_threads):
        if i == 0:
            catalog_list.append(coordinates[i * Ncut : (i+1) * Ncut])
            idx_stop.append((i+1) * Ncut)
        elif i > 0 and i < nb_threads-1:
            catalog_list.append(coordinates[i * Ncut : (i+1) * Ncut])
            idx_stop.append((i+1) * Ncut)
        elif i == nb_threads-1:
            catalog_list.append(coordinates[i * Ncut : Ncat])
            idx_stop.append(Ncat)
        if idx_stop[-1] >= Ncat:
            break
    pool = mp.Pool(nb_threads)
    # call apply_async() without callback
    result_objects = [pool.apply_async(_run_xmatch,
                                args=(cat, catalog, radius))
                          for cat in catalog_list]

    # result_objects is a list of pool.ApplyResult objects
    results = [r.get() for r in result_objects]

    # Don't forget to close
    pool.close()
    pool.join()

    # If one table is empty and one returns something,
    # there will be a conflict type, str vs something.
    # So keep only the one with data
    res2keep = []
    c = 0
    for i in range(len(results)):
        if len(results[i]) > 0:
            res2keep.append(results[i])
            c += 1
    # If all are empty select the first not to crash the code
    if c == 0:
        res2keep.append(results[0])

    crossmatch = vstack(res2keep)
    return crossmatch


def catalogs(image_table, radius,
             catalogs=["I/345/gaia2", "II/349/ps1", "I/271/out", "I/284/out"],
             Nb_cuts=(1, 1), subFiles=None, nb_threads=4):
    """
    Input file is *.magwcs and the output is the list of the stars *.oc
    which were not identified in the catalogue
    radius is expressed in pixels
    """

    cat_dict = {
        "I/284/out": "USNO-B1",
        "I/345/gaia2": "GAIA DR2",
        "II/349/ps1": "PS1 DR1",
        "I/271/out": "GSC",
    }
    counter = 0
    _filename_list = []

    if subFiles is not None:
        subfiles = np.array(subFiles)
        # filelist = [im for im in subfiles[:, 0]]
        # Rather take the original data as the registartion introduces artefact
        filelist = [im for im in image_table["filenames"]]
        filelist.extend([im for im in subfiles[:, 2]])
    else:
        filelist = image_table["filenames"]
    for i, filename in enumerate(filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + "/"
        else:
            folder = ""
        #  Get rid of the extension to keep only the name
        filename2,extension = os.path.splitext(filename_ext)

        magfilewcs = folder + filename2 + ".magwcs"

        if Nb_cuts == (1, 1):
            original_name = folder + filename2
            quadrant = 1
        else:
            split_file = filename2.split("_Q")
            original_name = folder + split_file[0]
            for name in split_file[1].split("_")[1:]:
                original_name = original_name + "_" + name
            original_name = os.path.splitext(original_name)[0]
            quadrant = split_file[1].split("_")[0]
            """
            if len(split_file[-1]) == 2:
                idx = 3
            elif len(split_file[-1]) == 3:
                idx = 4
            """
        _filename_list.append(original_name + ".oc")

        header = fits.getheader(filename)
        # Get pixel scale in degrees
        try:
            pixScale = abs(header["CDELT1"])
        except Exception:
            try:
                pixScale = abs(header["CD1_1"])
            except Exception:
                print(
                    "Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1"
                )

        # Load detected sources in astropy table
        detected_sources = ascii.read(
            magfilewcs,
            names=[
                "Xpos",
                "Ypos",
                "_RAJ2000",
                "_DEJ2000",
                "mag_inst",
                "mag_inst_err",
                "edge",
                "psf_chi2",
                "psf_mag",
                "psf_magerr",
                "FWHM",
                "FWHMPSF",
                "filenames",
                "FlagSub",
                "OriginalIma",
                "RefIma",
            ],
            format="commented_header",
        )
        if detected_sources:
            detected_sources["quadrant"] = [quadrant] * len(detected_sources)
            # Do not need it as the astrometric calibration
            # is performed on each quadrant now.
            """
            # Transform X_pos and Y_pos to original image in case it was split
            header2 = fits.getheader(filename)
            Naxis1 = float(header2['NAXIS1'])
            Naxis2 = float(header2['NAXIS2'])
            Naxis11 = int(Naxis1/Nb_cuts[0])
            Naxis22 = int(Naxis2/Nb_cuts[1])
            if Nb_cuts == (1,1):
                quad = None
                index_i = 0
                index_j = 0
            else:
                quad, index_i, index_j = image_table['quadrant'][i].split('_')
                quad = quad[1:]
            detected_sources['Xpos'] = detected_sources['Xpos'] + Naxis22 * int(index_j)
            detected_sources['Ypos'] = detected_sources['Ypos'] + Naxis11 * int(index_i)
            """
            detected_sources["Xpos_quad"] = detected_sources["Xpos"]
            detected_sources["Ypos_quad"] = detected_sources["Ypos"]

            if i == 0:
                detected_sources_tot = deepcopy(detected_sources)
            else:
                detected_sources_tot = vstack(
                    [detected_sources_tot, detected_sources])
    # Add units
    detected_sources_tot["_RAJ2000"] *= u.deg
    detected_sources_tot["_DEJ2000"] *= u.deg
    # Add index for each source
    detected_sources_tot["idx"] = np.arange(len(detected_sources_tot))

    #  Add a flag to indicate whether a source is crossmatched with
    #  a referenced source
    detected_sources_tot["Match"] = ["N"] * len(detected_sources_tot)

    #  Initialise candidates with all detected sources
    #candidates = deepcopy(detected_sources_tot)
    candidates = deepcopy(detected_sources_tot["_RAJ2000", "_DEJ2000",
                                               "idx", "Match", "FlagSub"])
    # candidates.write('test0.dat', format='ascii.commented_header', overwrite=True)
    print("\nCrossmatching sources with catalogs.")
    print(
        "Radius used for crossmatching with catalogs: %.2f arcseconds\n"
        % (radius * pixScale * 3600)
    )

    mask_matched = candidates["Match"] == "N"
    for catalog in catalogs:
        print(catalog, len(candidates[mask_matched]))
        # Use Xmatch to crossmatch with catalog
        crossmatch = run_xmatch(
            candidates[mask_matched], catalog,
            radius * pixScale * 3600, nb_threads
        )

        # crossmatch.write('test.dat', format='ascii.commented_header', overwrite=True)
        # Do not consider duplicates
        #  Meaning that if there are several sources
        #  we consider it as a crossmatch
        referenced_star_idx = np.unique(crossmatch["idx"])

        candidates["Match"][referenced_star_idx] = "Y"

        #  Update Match mask
        mask_matched = candidates["Match"] == "N"

        if subFiles is not None:
            mask = (candidates["FlagSub"] == "Y") & (mask_matched)
            mask2 = candidates["FlagSub"] == "Y"
            print(
                "%d/%d candidates left in substracted image after crossmatching with %s"
                % (len(candidates[mask]), len(candidates[mask2]), cat_dict[catalog])
            )
            mask = (candidates["FlagSub"] == "N") & (mask_matched)
            mask2 = np.invert(mask2)
            print(
                "%d/%d candidates left in original image (without substraction) after crossmatching with %s"
                % (len(candidates[mask]), len(candidates[mask2]), cat_dict[catalog])
            )
        else:
            print(
                "%d/%d candidates left after crossmatching with %s"
                % (len(candidates[mask_matched]), len(candidates), cat_dict[catalog])
            )
        if len(candidates[mask_matched]) == 0:
            break

    # Flag sources without any match in catalogs
    idx_no_match = np.isin(detected_sources_tot['idx'],
                           candidates['idx'][mask_matched])
    #candidates = deepcopy(detected_sources_tot[keep_idx])
    # Add the Match flag in original table data
    detected_sources_tot['Match'] = candidates['Match']

    #  Get filename
    _filename = np.unique(_filename_list)[0]
    # Write candidates file.
    # If substraction was performed, split transients into specific files
    if subFiles is not None:
        mask = (detected_sources_tot["FlagSub"] == "Y") & idx_no_match
        detected_sources_tot[mask].write(
            os.path.splitext(_filename)[0] + "_sub.oc",
            format="ascii.commented_header",
            overwrite=True,
        )
        oc = detected_sources_tot[mask]["_RAJ2000", "_DEJ2000"]
        oc.write(
            os.path.splitext(_filename)[0] + "_sub.oc_RADEC",
            format="ascii.commented_header",
            overwrite=True,
        )


    mask = (detected_sources_tot["FlagSub"] == "N") & idx_no_match
    detected_sources_tot[mask].write(
        _filename,
        format="ascii.commented_header",
        overwrite=True)
    oc = detected_sources_tot[mask]["_RAJ2000", "_DEJ2000"]
    oc.write(
        os.path.splitext(_filename)[0] + ".oc_RADEC",
        format="ascii.commented_header",
        overwrite=True,
    )

    #  Also write a file with all the sources detected to know
    detected_sources_tot.write(
        os.path.splitext(_filename)[0] + ".alldetections",
        format="ascii.commented_header",
        overwrite=True,
    )

    # oc=candidates['Xpos', 'Ypos']
    # oc.write(magfilewcs+'4',format='ascii.commented_header', overwrite=True)

    return detected_sources_tot


def skybot(ra_deg, dec_deg, date, radius, Texp):
    """
    query SkyBoT catalog using astroquery.skybot to search for moving objects
    parameters: ra_deg, dec_deg, date, Texp, radius:
                ra_deg, dec_deg: RA and DEC in degrees astropy.units
                date: date of observation (UTC)
                Texp: exposure time in astropy.unit
                radius in arcsecond
    returns: astropy.table object
    """
    #  if present, substitue T by a whitespace
    if "T" in date:
        date = date.replace("T", " ")

    field = SkyCoord(ra_deg, dec_deg)
    epoch = Time(str(date), format="iso")
    moving_objects_list = Skybot.cone_search(
        field, radius, epoch, position_error=120 * u.arcsecond
    )

    #  Add RA, DEC at end of exposure
    moving_objects_list["RA_Tend"] = (
        moving_objects_list["RA"] + moving_objects_list["RA_rate"] * Texp
    )
    moving_objects_list["DEC_Tend"] = (
        moving_objects_list["DEC"] + moving_objects_list["DEC_rate"] * Texp
    )
    #  Compute angular distance  between Tstart and Tend
    c1 = SkyCoord(
        moving_objects_list["RA"],
        moving_objects_list["DEC"],
        frame="icrs")
    c2 = SkyCoord(
        moving_objects_list["RA_Tend"], moving_objects_list["DEC_Tend"], frame="icrs"
    )
    sep = c1.separation(c2)
    moving_objects_list["RADEC_sep"] = sep

    return moving_objects_list


def crossmatch_skybot(sources, moving_objects, radius=5):
    """
    crossmatch list of detected sources with list of moving objects
    in this field using skybot

    parameters: sources, moving_objects, radius:
                sources: astropy.table containing list of unknown
                         sources detected
                moving_objects: astropy.table containing list of moving
                                objects from skybot in the same field of view
                radius in arcsecond
    returns: astropy.table object

    NOT WORKING WITH PYTHON2.7, seems ok WITH PYTHON3
    """
    # sources=ascii.read('/home/corre/codes/Tests/%s.fits.oc' % filename)
    # moving_objects=ascii.read('/home/corre/codes/gmadet/gmadet/moving_objects.dat')
    cat1 = np.empty((len(sources), 2), dtype=np.float64)
    cat2 = np.empty((len(moving_objects), 2), dtype=np.float64)
    cat1[:, 0] = sources["_RAJ2000"]
    cat1[:, 1] = sources["_DEJ2000"]

    cat2[:, 0] = moving_objects["RA"]
    cat2[:, 1] = moving_objects["DEC"]
    dist, ind = crossmatch_angular(cat1, cat2, radius / 3600)
    match = ~np.isinf(dist)
    # print (match)
    dist_match = dist[match]
    #  Convert in arcseconds
    dist_match *= 3600
    # print (dist_match)
    if dist_match:
        mov_match = sources[ind[np.unique(match)]]
        print(mov_match)
    candidates = sources[match == False]

    return candidates


def moving_objects(filelist):
    """
    Crossmatch the list of candidates with moving objects using SkyBoT

    """
    for filename in filelist:
        header = fits.getheader(filename)
        ra_deg = float(header["CRVAL1"]) * u.deg
        dec_deg = float(header["CRVAL2"]) * u.deg
        date = header["DATE-GPS"]
        radius = 2 * u.deg
        Texp = float(header["exposure"]) * u.second

        moving_objects = skybot(ra_deg, dec_deg, date, radius, Texp)

        moving_objects.write(
            "moving_objects.dat", format="ascii.commented_header", overwrite=True
        )
        moving_objects["RA", "DEC"].write(
            "moving_objects.reg", format="ascii.commented_header", overwrite=True
        )

        candidates_out = crossmatch_skybot(
            candidates, moving_objects, radius=10)

        print(
            "%d match with a moving object found"
            % (len(candidates) - len(candidates_out))
        )

    # return candidates_out
