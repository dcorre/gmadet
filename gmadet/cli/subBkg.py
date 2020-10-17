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
    list_files,
    make_results_dir,
)

from gmadet.background import bkg_estimation

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        usage="usage: %(prog)s [options] data",
        description="Remove background in astronomical images with photutils."
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
        "--box",
        dest="box",
        required=False,
        type=int,
        nargs='+',
        default=[20, 20],
        help="Size of the background mesh. (Default: 20 20)"
    )

    parser.add_argument(
        "--filter-size",
        dest="filter_size",
        required=False,
        type=int,
        nargs='+',
        default=[3, 3],
        help="Size of the filter box. (Default: 3 3)"
    )

    parser.add_argument(
        "--bkg-estimator",
        dest="bkg_estimator",
        required=False,
        type=str,
        choices=["SExtractor", "MMM", "Median", "Mean"],
        default="SExtractor",
        help="Background estimator method among the available ones. Check "
             "https://photutils.readthedocs.io/en/stable/background.html for"
             " more info. (Default: SExtractor)."
    )

    parser.add_argument(
        "--sigma-clip",
        dest="sigma_clip",
        required=False,
        type=float,
        default=3.0,
        help="The number of standard deviations to use for both the lower "
             "and upper clipping limit. These limits are overridden by "
             "`sigma_lower` and `sigma_upper`, if input. (Default: 3.0)"
    )

    parser.add_argument(
        "--sigma-lower",
        dest="sigma_lower",
        required=False,
        type=float,
        default=None,
        help="The number of standard deviations to use as the lower bound "
             "for the clipping limit. (Default: None)"
    )

    parser.add_argument(
        "--sigma-upper",
        dest="sigma_upper",
        required=False,
        type=float,
        default=None,
        help="The number of standard deviations to use as the upper bound "
             "for the clipping limit. (Default: None)"
    )

    parser.add_argument(
        "--maxiters",
        "--max-iters",
        dest="maxiters",
        required=False,
        type=int,
        default=10,
        help="The maximum number of sigma-clipping iterations to perform or "
             "None to clip until convergence is achieved (i.e., iterate until "
             "the last iteration clips nothing). (Default: 10)"
    )

    args, filenames = parser.parse_known_args()

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

        bkg_estimation(filename, box=args.box, filter_size=args.filter_size,
                       bkg_estimator=args.bkg_estimator, sigma=args.sigma_clip,
                       sigma_lower=args.sigma_lower, sigma_upper=args.sigma_upper,
                       maxiters=args.maxiters, outLevel=2)


if __name__ == "__main__":
    main()
