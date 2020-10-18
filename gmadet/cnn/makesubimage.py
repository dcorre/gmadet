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
import copy
import numpy as np
from gmadet.utils import (getpath, make_sub_image,
                          mv_p, rm_p, mkdir_p)
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.coordinates import SkyCoord

def getCandPos(path, pattern=".alldetections", flag_notsub=False):
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
    for resfile in resfiles:
        data = ascii.read(resfile)
        # Select only detection in the substracted images and not
        # close to the edge.
        if not flag_notsub:
            mask = (data["FlagSub"] == "Y") & (data["edge"] == "N")
        else:
            mask = (data["FlagSub"] == "N") & (data["edge"] == "N")
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
        origfile = os.path.basename(data["OriginalIma"][mask][0])
        origfilelist.extend([origfile] * len(data[mask]))
    cand_id = np.arange(len(RA_list))
    table = Table(
        [
            cand_id,
            RA_list,
            Dec_list,
            Xpos_list,
            Ypos_list,
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
            "Band",
            "mag",
            "magerr",
            "FWHM",
            "FWHMPSF",
            "filename",
            "OriginalIma",
        ],
    )

    table.write(path + "/candidates_list.dat",
                format="ascii.commented_header",
                overwrite=True)

    return table


def crossmatch_detections(path, candidates_list, radius=1):
    """ Crossmatch detections with simulated events positions """

    # Load simulated events summary file
    sim_list = ascii.read(os.path.join(path, "simulated_objects.list"))
    newfilename = []
    for sim in sim_list:
        filename = os.path.basename(sim["filename"])
        newfilename.append(filename)
    sim_list["filename2"] = newfilename

    # Check whether the simulated events have been detected
    Ndet = []
    candID_list = []
    closest_candID_list = []

    Nsim = len(sim_list)
    for i, event in enumerate(sim_list):
        print("processing simulated event  %d/%d ..." % (i, Nsim), end="\r", flush=True)

        event_coords = [event["RA"], event["Dec"]]
        # Check whether it is a simulated object
        # FIXME: it is currently broken, as candidates_list contain _reg_ annotated filenames
        # when simulated_objects.list does not
        mask1 = candidates_list["OriginalIma"] == event["filename2"]
        candidates = copy.deepcopy(candidates_list) # [mask1])

        # Compute the separation with detections and sort by ascending values
        offset = SkyCoord(candidates["RA"], candidates["Dec"], unit='deg').separation(
            SkyCoord(event["RA"], event["Dec"], unit='deg')).degree
        candidates["offset"] = offset
        candidates.sort("offset")
        mask2 = candidates["offset"] < radius / 3600
        Nmatch = len(candidates[mask2])
        Ndet.append(Nmatch)
        if len(candidates[mask2]) > 0:
            text = ""
            for j in list(candidates["ID"][mask2]):
                text += "%s" % j
                if j < len(candidates[mask2]):
                    text += ","
            candID_list.append(text)
            closest_candID_list.append(candidates["ID"][mask2][0])
        else:
            candID_list.append(None)
            closest_candID_list.append(None)

    # Add a column 'detected' to the simulated events list table
    sim_list["Nmatches"] = Ndet
    sim_list["closest_candID"] = closest_candID_list
    sim_list["all_candIDs"] = candID_list

    sim_list.write(os.path.join(path, "crossmatch.dat"),
                   format="ascii.commented_header",
                   overwrite=True)
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
            # mask2 = (sim_list['RA'] - cand['RA'])**2 + (sim_list['Dec'] - cand['Dec'])**2 < (radius/3600)**2
            # mask = np.bitwise_and(mask1,mask2)

            if len(sim_list[mask]) == 1:
                outdir = truedir
            else:
                if false:
                    outdir = falsedir
                else:
                    outdir = resdir
            """
            inputname = (
                path_gmadet
                + "/cnn/data/sim/%s/images/gmadet_results/substraction/" % telescope
                + cand["filename"].split("/")[-1]
            )
            """
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
            make_sub_image(
                inputname,
                OT_coords,
                coords_type="world",
                output_name=outname,
                size=[size, size],
                fmt="fits",
                addheader=False,
            )
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
            hdus.writeto(outname, overwrite=True)

        except BaseException:
            print("Could not extract candidate in %s" % inputname)
