#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import os
import argparse
import warnings
import subprocess

from gmadet.utils import (
    load_config,
    list_files,
    make_results_dir,
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
        usage="usage: %(prog)s [options] data",
        description="Insert simulated point like sources in astronomical "
                    "images using the estimated PSFs of each image."

    )

    parser.add_argument(
        "--results",
        dest="path_results",
        required=False,
        type=str,
        default='gmadet_sim',
        help="Base path to store the results. "
             "(Default: gmadet_sim)"
    )

    parser.add_argument(
        "--keep-old",
        "--keep",
        dest="keep",
        required=False,
        action="store_true",
        help="Keep previous results"
    )

    parser.add_argument(
        "--skip-processed",
        "--skip",
        dest="skip",
        required=False,
        action="store_true",
        help="Skip already processed files"
    )

    parser.add_argument(
        "--preprocess",
        dest="preprocess",
        required=False,
        type=str,
        default=None,
        help="Pre-process the image using external program before analysing. "
             "The program should accept two positional arguments - original "
             "filename and new one. (Default: just copy the image)"
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
        "--ntrans",
        "--num-trans",
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
        "--mag-range",
        dest="magrange",
        required=False,
        type=int,
        nargs='+',
        default=[14, 23],
        help="Magnitude range of simulated sources. (Default: 14 23)"
    )

    parser.add_argument(
        "--zp",
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
        "--astrometry",
        dest="doAstrometry",
        required=False,
        default="scamp",
        choices=["no", "scamp"],
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
        "--iter-max",
        dest="itermax",
        required=False,
        type=float,
        default=5,
        help="Max number of iteration to reach the required accuracy. "
             "(Default: 5)",
    )

    parser.add_argument(
        "--conv-filter",
        dest="convFilter",
        required=False,
        default="default",
        type=str,
        help="Corresponds to FILTER_NAME keyword for sextractor "
             "(without .conv)."
             "\nDifferent filter available listed here: %s"
        % path_gmadet + "/config/conv_kernels/"
             "\n(Default: default)",
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

    args, filenames = parser.parse_known_args()

    # Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)

    filenames = list_files(filenames, exclude=args.path_results)

    for raw_filename in filenames:
        filename = make_results_dir(
            raw_filename,
            outputDir=args.path_results,
            keep=args.keep,
            skip=args.skip,
            copy=False if args.preprocess else True
        )

        if not filename:
            print("%s is already processed, skipping. \n" % raw_filename)
            continue

        if args.preprocess:
            # We need to call external code what will copy (processed)
            # image to results dir
            print("Pre-processing %s" % raw_filename)
            subprocess.call(args.preprocess.split() + [raw_filename,
                                                       filename])

            if not os.path.exists(filename):
                print("Pre-processing failed")
                continue

        if args.doAstrometry != "no":
            print("Sanitise header and data of %s.\n" % filename)
            sanitise_fits(filename)
            astrometric_calib(
                filename,
                config,
                soft=args.doAstrometry,
                verbose=args.verbose,
                accuracy=args.accuracy,
                itermax=args.itermax,
            )

        # Estimate the PSF FWHM for each image/quadrants using psfex
        FWHM_list = psfex(
            filename,
            config,
            verbose=args.verbose,
            outLevel=2,
        )

        datapath = os.path.split(filename)[0]
        sim(datapath,
            filename,
            Ntrans=args.Ntrans,
            size=args.size,
            magrange=args.magrange,
            gain=args.gain,
            magzp=args.ZP
            )


if __name__ == "__main__":
    main()
