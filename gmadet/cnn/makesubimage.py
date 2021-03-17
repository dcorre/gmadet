#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: David Corre (IJCLab/CNRS)
"""

import errno
import glob
import math
import os
import shutil
import subprocess
import sys
import argparse
from copy import deepcopy
import numpy as np
from gmadet.utils import getpath, make_sub_image, make_fits, mv_p, rm_p, mkdir_p
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord


def getCandPos(path, pattern=".alldetections", flag_notsub=False):
    """Combine all the .alldetections file in the provided path."""

    resfiles = glob.glob(os.path.join(path, "**", "*%s*" % pattern), recursive=True)
    filelist = []
    origfilelist = []
    Xpos_list = []
    Ypos_list = []
    RA_list = []
    Dec_list = []
    mag_list = []
    magerr_list = []
    band_list = []
    FWHM = []
    FWHM_PSF = []
    edge_list = []
    for resfile in resfiles:
        data = ascii.read(resfile)
        # Select only detection in the substracted images and not
        # close to the edge.
        if not flag_notsub:
            mask = data["FlagSub"] == "Y"  # & (data["edge"] == "N")
        else:
            mask = data["FlagSub"] == "N"  # & (data["edge"] == "N")
        RA_list.extend(data["_RAJ2000"][mask])
        Dec_list.extend(data["_DEJ2000"][mask])
        Xpos_list.extend(data["Xpos"][mask])
        Ypos_list.extend(data["Ypos"][mask])
        mag_list.extend(data["mag_calib"][mask])
        magerr_list.extend(data["mag_calib_err"][mask])
        band_list.extend(data["filter_cat"][mask])
        filelist.extend(data["filenames"][mask])
        FWHM.extend(data["FWHM"][mask])
        FWHM_PSF.extend(data["FWHMPSF"][mask])
        edge_list.extend(data["edge"][mask])
        origfile = os.path.abspath(data["OriginalIma"][mask][0])
        origfilelist.extend([origfile] * len(data[mask]))
    cand_id = np.arange(len(RA_list))
    table = Table(
        [
            cand_id,
            RA_list,
            Dec_list,
            Xpos_list,
            Ypos_list,
            edge_list,
            band_list,
            mag_list,
            magerr_list,
            FWHM,
            FWHM_PSF,
            filelist,
            origfilelist,
        ],
        names=[
            "ID",
            "RA",
            "Dec",
            "Xpos",
            "Ypos",
            "edge",
            "Band",
            "mag",
            "magerr",
            "FWHM",
            "FWHMPSF",
            "filename",
            "OriginalIma",
        ],
    )

    table.write(
        path + "/candidates_list.dat", format="ascii.commented_header", overwrite=True
    )

    return table


def crossmatch_detections(path, candidates_list, radius=1):
    """ Crossmatch detections with simulated events positions """

    # Combine all simulated_objects.list
    pattern = "simulated_objects.list"
    files = glob.glob(os.path.join(path, "**", "*%s*" % pattern), recursive=True)

    sim_list = None
    for f in files:
        data = ascii.read(f)
        if sim_list is None:
            sim_list = deepcopy(data)
        else:
            sim_list = vstack([sim_list, data])

    # Load simulated events summary file
    newfilename = []
    for sim in sim_list:
        filename = os.path.basename(sim["filename"])
        newfilename.append(filename)
    sim_list["filename2"] = newfilename

    # Check whether the simulated events have been detected
    Ndet = []
    candID_list = []
    closest_candID_list = []
    edge_list = []

    Nsim = len(sim_list)
    for i, event in enumerate(sim_list):
        print("processing simulated event  %d/%d ..." % (i, Nsim), end="\r", flush=True)

        event_coords = [event["RA"], event["Dec"]]
        # Check whether it is a simulated object
        # FIXME: it is currently broken, as candidates_list contain
        # _reg_ annotated filenames when simulated_objects.list does not.
        # Answer:
        # This mask was used to use only one image at a time to speed up
        # this loop. But now that the original image in .alldetection
        # is the registered image and not the true original image, it
        # does not work anymore as the filename in simulated_object.list
        # point to the true original image. Can either add a new column
        # in the .alldetection to have the true original image, or find
        # an another solution to speed up the loop.
        # For now we just ignore it.
        # mask1 = candidates_list["OriginalIma"] == event["filename2"]
        candidates = deepcopy(candidates_list)  # [mask1])

        # Compute the separation with detections and sort by ascending values
        offset = (
            SkyCoord(candidates["RA"], candidates["Dec"], unit="deg")
            .separation(SkyCoord(event["RA"], event["Dec"], unit="deg"))
            .degree
        )

        candidates["offset"] = offset
        candidates.sort("offset")
        mask2 = candidates["offset"] < radius / 3600
        Nmatch = len(candidates[mask2])
        Ndet.append(Nmatch)
        if len(candidates[mask2]) > 0:
            text = ""
            for c, j in enumerate(list(candidates["ID"][mask2])):
                text += "%s" % j
                if c + 1 < len(candidates[mask2]):
                    text += ","
            candID_list.append(text)
            closest_candID_list.append(candidates["ID"][mask2][0])
            edge_list.append(candidates["edge"][mask2][0])
        else:
            candID_list.append(None)
            closest_candID_list.append(None)
            edge_list.append(None)

    # Add a column 'detected' to the simulated events list table
    sim_list["Nmatches"] = Ndet
    sim_list["closest_candID"] = closest_candID_list
    sim_list["all_candIDs"] = candID_list
    sim_list["edge"] = edge_list
    sim_list.write(
        os.path.join(path, "crossmatch.dat"),
        format="ascii.commented_header",
        overwrite=True,
    )
    return sim_list


def subimage(path, training, size=32, radius=1, flag_notsub=False, false=False):
    """ Extract a sub-image centered on the candidate position """

    path_gmadet = getpath()

    # size of the extracted image
    cutsize = (size, size)
    print("Combine the detections from all simulated images.")
    candidates_list = getCandPos(path, flag_notsub=flag_notsub)

    if training:
        print("Crossmatch simulated events with detections")
        sim_list = crossmatch_detections(path, candidates_list, radius=radius)

        resdir = os.path.join(path, "candidates_training")
        mkdir_p(resdir)
        truedir = os.path.join(path, "candidates_training", "true")
        mkdir_p(truedir)
        falsedir = os.path.join(path, "candidates_training", "false")
        mkdir_p(falsedir)
    else:
        resdir = os.path.join(path, "candidates")
        mkdir_p(resdir)

    for i, cand in enumerate(candidates_list):
        print(
            "processing candidates %d/%d ..." % (i, len(candidates_list)),
            end="\r",
            flush=True,
        )

        OT_coords = [cand["RA"], cand["Dec"]]

        if training:
            # Check if corresponds to a simulated event
            mask = sim_list["closest_candID"] == cand["ID"]
            # Check whether it is a simulated object
            # mask1 = sim_list['filename2'] == cand['OriginalIma']
            # mask2 = (sim_list['RA'] - cand['RA'])**2 + \
            #        (sim_list['Dec'] - cand['Dec'])**2 < (radius/3600)**2
            # mask = np.bitwise_and(mask1,mask2)

            # Consider only sources with a single match.
            # Some real sources close to a cosmic or another will not be
            # consider. But with a visual inspection they can be easily
            # classified as true transient.
            if len(sim_list[mask]) == 1:
                outdir = truedir
            else:
                if false:
                    outdir = falsedir
                else:
                    outdir = resdir
            inputname = cand["filename"]
        else:
            outdir = resdir
            inputname = cand["filename"]

        outname = os.path.join(outdir, "candidate_%d.fits" % i)
        # If the inputname can not be found for some reasons
        # Make sure the code is not crashing
        try:
            # Get initial image size
            hdr_input = fits.getheader(inputname)
            Naxis1 = hdr_input["NAXIS1"]
            Naxis2 = hdr_input["NAXIS2"]

            # Extract small image
            subimage, header, size_list, pixref, origin = make_sub_image(
                inputname,
                OT_coords,
                coords_type="world",
                sizes=[size, size],
            )

            make_fits(subimage[0], outname, header[0], size_list[0], pixref[0])

            # add information to header
            hdus = fits.open(outname, memmap=False)
            hdr = hdus[0].header
            hdr["MAG"] = cand["mag"]
            hdr["MAGERR"] = cand["magerr"]
            hdr["FWHM"] = cand["FWHM"]
            hdr["FWHMPSF"] = cand["FWHMPSF"]
            hdr["FILTER"] = cand["Band"]
            hdr["RA"] = cand["RA"]
            hdr["Dec"] = cand["Dec"]
            hdr["Xpos"] = cand["Xpos"]
            hdr["Ypos"] = cand["Ypos"]
            hdr["FILE"] = cand["filename"]
            hdr["CANDID"] = cand["ID"]
            # Whether it is close to the edge of the image
            # If yes the image will not be size x size in pixels
            """
            if (
                (cand["Xpos"] > Naxis1 - size)
                or (cand["Xpos"] < size)
                or (cand["Ypos"] > Naxis2 - size)
                or (cand["Ypos"] < size)
            ):
                hdr["edge"] = "True"
                print("Edge ", outname)
            else:
                hdr["edge"] = "False"
            """
            # If source too close to the edge, the cutout has dimension
            # smaller than the reauired (size x size)
            # It will cause the CNN to crash so flag them.
            if hdus[0].data.shape != (size, size):
                hdr["edge"] = "True"
            else:
                hdr["edge"] = "False"
            hdus.writeto(outname, overwrite=True)

        except BaseException:
            print("Could not extract candidate in %s" % inputname)
