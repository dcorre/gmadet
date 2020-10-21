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
from gmadet.cnn.makesubimage import (
        getCandPos,
        crossmatch_detections)

def fraction_detected():
    """Number of simulated events detected"""


def makestats(path, radius=2):
    """ Create some statistics on the simulated events """
    mkdir_p(os.path.join(path, 'CheckSim'))
    candidates_list = getCandPos(path)

    # Load file crossmatching the detected candidates with simulated events
    # Try to load file if already exists
    # Otherwise create it
    try:
        crossmatch = ascii.read(os.path.join(path, "crossmatch.dat"))
    except BaseException:
        crossmatch = crossmatch_detections(
            path, candidates_list, radius=radius)

    # mask_det = crossmatch["Nmatches"] == 1
    # Actually take everything with a match.
    mask_det = crossmatch["Nmatches"] > 0

    bands = crossmatch.group_by("filter").groups.keys

    # plot histogram of simulated sources magnitude
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
    plt.savefig(os.path.join(path, "CheckSim/sim_mag_distrib_allbands.png"))

    # plot fraction of detection for a given magnitude bin
    plt.figure()
    x = (bins1[1:]+bins1[:-1])/2
    plt.plot(x, n2/n1)
    plt.xlabel('Simulated magnitude')
    plt.ylabel('Fraction (%)')
    plt.title('Fraction of simulated events detected.')
    plt.savefig(os.path.join(path, 'CheckSim/detection_fraction_allbands.png'))

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
        plt.savefig(os.path.join(
            path, "CheckSim/sim_mag_distrib_%s.png" % band[0]))
