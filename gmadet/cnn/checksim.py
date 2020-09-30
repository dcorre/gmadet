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
from astropy.io import ascii, fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from gmadet.utils import mkdir_p

def getCandPos(path, pattern=".alldetections"):

    #  Try to load if already exists
    #  Otherwise create it
    try:
        table = ascii.read(path + "/candidates_list.dat")
    except BaseException:
        resfiles = glob.glob(
            path + "/**/gmadet_results/*%s*" %
            pattern, recursive=True)
        filelist = []
        origfilelist = []
        Xpos_list = []
        Ypos_list = []
        RA_list = []
        Dec_list = []
        mag_list = []
        magerr_list = []
        band_list = []
        for resfile in resfiles:
            data = ascii.read(resfile)
            mask = (data["FlagSub"] == "Y") & (data["edge"] == "N")
            RA_list.extend(data["_RAJ2000"][mask])
            Dec_list.extend(data["_DEJ2000"][mask])
            Xpos_list.extend(data["Xpos"][mask])
            Ypos_list.extend(data["Ypos"][mask])
            mag_list.extend(data["mag_calib"][mask])
            magerr_list.extend(data["mag_calib_err"][mask])
            band_list.extend(data["filter_cat"][mask])
            filelist.extend(data["filenames"][mask])
            _, origfile = os.path.split(data["OriginalIma"][mask][0])
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
                "filename",
                "OriginalIma",
            ],
        )
        table.write(
            path + "/candidates_list.dat",
            format="ascii.commented_header",
            overwrite=True,
        )

    return table


def crossmatch_detections(path, candidates_list, radius=1, getFilter=False):
    """ Crossmatch detections with simulated events positions """

    #  Try to load file if already exists
    #  Otherwise create it
    try:
        sim_list = ascii.read(path + "/crossmatch.dat")
    except BaseException:

        #  Load simulated events summary file
        sim_list = ascii.read(path + "/simulated_objects.list")
        newfilename = []
        for sim in sim_list:
            filename = sim["filename"].split("/")[-1]
            newfilename.append(filename)
        sim_list["filename2"] = newfilename

        #  Check whether the simulated events have been detected
        Ndet = []
        candID_list = []
        closest_candID_list = []

        Nsim = len(sim_list)
        for i, event in enumerate(sim_list):
            print(
                "processing simulated event  %d/%d ..." % (i, Nsim),
                end="\r",
                flush=True,
            )
            event_coords = [event["RA"], event["Dec"]]
            #  Check whether it is a simulated object
            mask1 = candidates_list["OriginalIma"] == event["filename2"]
            candidates = copy.deepcopy(candidates_list[mask1])
            #  Compute the sepaation with detections and sort by ascending
            # values
            offset = (candidates["RA"] - event["RA"]) ** 2 + (
                candidates["Dec"] - event["Dec"]
            ) ** 2
            offset = np.sqrt(offset)
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

        #  Add a column 'detected' to the simulated events list table
        sim_list["Nmatches"] = Ndet
        sim_list["closest_candID"] = closest_candID_list
        sim_list["all_candIDs"] = candID_list
        sim_list.write(path + "/crossmatch.dat",
                       format="ascii.commented_header",
                       overwrite=True)
    # sim_list.show_in_browser()
    return sim_list

def fraction_detected():
    """Number of simulated events detected"""




def makestats(path, radius=2):
    """ Create some statistics on the simulated events """
    mkdir_p(path+'/CheckSim/')
    candidates_list = getCandPos(path)

    #  Load file crossmatching the detected candidates with simulated events
    crossmatch = crossmatch_detections(
        path, candidates_list, radius=radius, getFilter=True
    )

    mask_det = crossmatch["Nmatches"] == 1

    bands = crossmatch.group_by("filter").groups.keys

    #  plot histogram of simulated sources magnitude
    plt.figure()
    n1, bins1, patches = plt.hist(
        crossmatch["mag"],
        30,
        facecolor="C0",
        alpha=0.75,
        density=False,
        stacked=False,
        label="Sim",
    )
    n2, bins2, patches = plt.hist(
        crossmatch["mag"][mask_det],
        bins1,
        facecolor="C1",
        alpha=0.5,
        density=False,
        stacked=False,
        label="Det",
    )
    plt.xlabel("mag")
    plt.ylabel("N")
    plt.title("All filters")
    plt.grid(True)
    plt.legend()
    plt.savefig(path + "/CheckSim/sim_mag_distrib_allbands.png")

    # plot fraction of detection for a given magnitude bin
    plt.figure()
    x = (bins1[1:]+bins1[:-1])/2
    plt.plot(x, n2/n1)
    plt.xlabel('Simulated magnitude')
    plt.ylabel('Fraction (%)')
    plt.title('Fraction of simulated events detected.')
    plt.savefig(path + '/CheckSim/detection_fraction_allbands.png')

    for band in bands:
        mask_band = crossmatch["filter"] == band[0]
        mask = np.bitwise_and(mask_det, mask_band)
        plt.figure()
        n1, bins1, patches = plt.hist(
            crossmatch["mag"][mask_band],
            30,
            facecolor="C0",
            alpha=0.75,
            density=False,
            stacked=False,
            label="Sim",
        )
        n2, bins2, patches = plt.hist(
            crossmatch["mag"][mask],
            bins1,
            facecolor="C1",
            alpha=0.5,
            density=False,
            stacked=False,
            label="Det",
        )
        plt.xlabel("mag")
        plt.ylabel("N")
        plt.title("%s band" % band[0])
        plt.grid(True)
        plt.legend()
        plt.savefig(path + "/CheckSim/sim_mag_distrib_%s.png" % band[0])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute sub-images centered on candidates position."
    )

    parser.add_argument(
        "--path", dest="path", required=True, type=str, help="Path to images"
    )

    parser.add_argument(
        "--radius",
        dest="radius",
        required=False,
        default=2,
        type=float,
        help="Radius for crossmatching detected sources with simulated events."
             " Default: 2 arcseconds",
    )

    args = parser.parse_args()

    makestats(args.path, radius=args.radius)
