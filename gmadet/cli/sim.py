#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import os
import argparse
import warnings

from gmadet.utils import (
    load_config,
    list_files,
    make_copy,
    getpath,
    getTel
)
from gmadet.sanitise import sanitise_fits
from gmadet.astrometry import astrometric_calib
from gmadet.psfex import psfex
from gmadet.cnn.sim import sim

warnings.simplefilter(action="ignore", category=FutureWarning)

def main():

    path_gmadet = getpath()

    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        description="Insert simulated point like sources in astronomical "
                    "images using the estimated PSFs of each image."

        )

    parser.add_argument(
        "--path_data",
        dest="path_data",
        required=True,
        type=str,
        help="Path to data"
    )

    parser.add_argument(
        "--telescope",
        dest="telescope",
        choices=telescope_list,
        required=True,
        type=str,
        help="Alias for the available telescopes.",
    )

    parser.add_argument(
        "--Ntrans",
        dest="Ntrans",
        required=False,
        type=int,
        default=100,
        help="Number of transients to insert in each image. "
             "(Defautl: 100)",
    )

    parser.add_argument(
        "--size",
        dest="size",
        required=False,
        type=int,
        default=48,
        help="Size of the transient image to insert in each image. "
             "Assumed to be a square, in pixels. (Defautl: 48)",
    )

    parser.add_argument(
        "--magrange",
        dest="magrange",
        required=False,
        type=int,
        nargs='+',
        default=[14,23],
        help="Magnitude range of simulated sources. (Default: 14 23)"
    )

    parser.add_argument(
        "--ZP",
        dest="ZP",
        required=False,
        type=int,
        default=30,
        help="Zeropoint used for the simulated sources. (Defautl: 30).",
    )

    parser.add_argument(
        "--gain",
        dest="gain",
        required=False,
        type=float,
        default=None,
        help="Gain to use, in e-/ADU. If None, take value from header. "
             " (Defautl: None).",
    )

    parser.add_argument(
        "--doAstrometry",
        dest="doAstrometry",
        required=False,
        default="scamp",
        choices=["No", "scamp"],
        type=str,
        help="Whether to perform astrometric calibration, with scamp. "
             "(Default: scamp)",
    )

    parser.add_argument(
        "--accuracy",
        dest="accuracy",
        required=False,
        type=float,
        default=0.15,
        help="Astrometric accuracy to reach, in arcseconds. "
             "(Defautl: 0.15 arcseconds)",
    )

    parser.add_argument(
        "--itermax",
        dest="itermax",
        required=False,
        type=float,
        default=5,
        help="Max number of iteration to reach the required accuracy. "
             "(Default: 5)",
    )

    parser.add_argument(
        "--convFilter",
        dest="convFilter",
        required=False,
        default="default",
        type=str,
        help="Corresponds to FILTER_NAME keyword for sextractor "
             "(without .conv)."
             "\nDifferent filter available listed here: %s" \
                     % path_gmadet + "/config/conv_kernels/"
             "\n(Default: default)"
        ,
    )

    parser.add_argument(
        "--verbose",
        dest="verbose",
        required=False,
        default="NORMAL",
        choices=["QUIET", "NORMAL", "FULL", "LOG"],
        type=str,
        help="Level of verbose, according to astromatic software. "
             "(Default: NORMAL)",
    )
    
    args = parser.parse_args()

    #  Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)

    filenames = list_files(args.path_data)
    #  copy original images
    #  Create list of the copy images
    filenames = make_copy(filenames, args.path_data,
                          outputDir="gmadet_sim/")

    for filename in filenames:

        if args.doAstrometry != "No": 
            print("Sanitise header and data of %s.\n" % filename)
            sanitise_fits(filename)
            astrometric_calib(
                filenames,
                config,
                soft=args.doAstrometry,
                verbose=args.verbose,
                accuracy=args.accuracy,
                itermax=args.itermax,
            )

        # Estimate the PSF FWHM for each image/quadrants using psfex
        FWHM_list = psfex(
                filenames,
                config,
                verbose=args.verbose,
                outLevel=2,
            )

        # Keep only the path if a file was provided
        if os.path.isdir(args.path_data):
            datapath = args.path_data
        else:
            text = os.path.split(args.path_data)
            if text[0]:
                datapath= text[0] + '/'
            else:
                datapath = ''
        sim(datapath,
            filenames,
            Ntrans=args.Ntrans,
            size=args.size,
            magrange=args.magrange,
            gain=args.gain,
            magzp=args.ZP
        )


if __name__ == "__main__":
    main()
