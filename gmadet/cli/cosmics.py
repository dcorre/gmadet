#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings
import os
import subprocess

from gmadet.utils import (
    load_config,
    list_files,
    make_results_dir,
    getpath,
    getTel
)
# from gmadet.psfex import psfex
from gmadet.remove_cosmics import run_lacosmic

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    # path_gmadet = getpath()
    # telescope_list = getTel()

    parser = argparse.ArgumentParser(
        usage="usage: %(prog)s [options] data",
        description="Remove cosmics in astronomical images."
    )

    parser.add_argument(
        "--results",
        dest="path_results",
        required=False,
        type=str,
        default='gmadet_results',
        help="Base path to store the results. "
             "(Default: gmadet_results)"
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
        "--contrast",
        dest="contrast",
        required=False,
        default=5.0,
        type=float,
        help="Contrast threshold between the Laplacian image and the "
             "fine-structure image. Check "
             "https://lacosmic.readthedocs.io/en/latest/api/lacosmic.lacosmic.html#lacosmic.lacosmic "
             "for more details. (Default: 5.0)",
    )

    parser.add_argument(
        "--cr-threshold",
        dest="cr_threshold",
        required=False,
        default=5.0,
        type=float,
        help="The Laplacian signal-to-noise ratio threshold for cosmic-ray "
             "detection. (Default: 5.0)",
    )

    parser.add_argument(
        "--neighbor-threshold",
        dest="neighbor_threshold",
        required=False,
        default=5.0,
        type=float,
        help="The Laplacian signal-to-noise ratio threshold for detection of "
             "cosmic rays in pixels neighboring the initially-identified "
             "cosmic rays. (Default: 5.0)",
    )

    parser.add_argument(
        "--maxiter",
        "--max-iter",
        dest="maxiter",
        required=False,
        default=4,
        type=int,
        help="The maximum number of iterations. (Default: 4)",
    )
    # Commented as the FWHM is not used anymore at the moment to automatise
    # the contrast value.
    """
    parser.add_argument(
        "--telescope",
        dest="telescope",
        choices=telescope_list,
        required=True,
        type=str,
        help="Alias for the available telescopes.",
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
        "--useweight",
        dest="useweight",
        action="store_true",
        help="If set, use weight map. "
             "Must be same name as image with .weight.fits extension. "
             "(Default: False)",
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
    """
    args, filenames = parser.parse_known_args()

    # Load config files for a given telescope
    # config = load_config(args.telescope, args.convFilter)
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

        FWHM = None
        run_lacosmic(filename, FWHM, contrast=args.contrast,
                     cr_threshold=args.cr_threshold,
                     neighbor_threshold=args.neighbor_threshold,
                     maxiter=args.maxiter, outLevel=2)


if __name__ == "__main__":
    main()
