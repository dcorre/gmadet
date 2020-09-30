#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

from gmadet.utils import (
    list_files,
    make_copy,
)

from gmadet.background import bkg_estimation

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        description="Remove background in astronomical images with photutils."
    )

    parser.add_argument(
        "--path_data",
        dest="path_data",
        required=True,
        type=str,
        help="Path to file"
    )

    parser.add_argument(
        "--box",
        dest="box",
        required=False,
        type=int,
        nargs='+',
        default=[20,20],
        help="Size of the background mesh. (Default: 20 20)"
    )

    parser.add_argument(
        "--filter_size",
        dest="filter_size",
        required=False,
        type=int,
        nargs='+',
        default=[3,3],
        help="Size of the filter box. (Default: 3 3)"
    )

    parser.add_argument(
        "--bkg_estimator",
        dest="bkg_estimator",
        required=False,
        type=str,
        choices=["SExtractor", "MMM", "Median", "Mean"],
        default="SExtractor",
        help="Background estimator method among the available ones. Check https://photutils.readthedocs.io/en/stable/background.html for more info. (Default: SExtractor)."
    )

    parser.add_argument(
        "--sigma_clip",
        dest="sigma_clip",
        required=False,
        type=float,
        default=3.0,
        help="The number of standard deviations to use for both the lower and upper clipping limit. These limits are overridden by `sigma_lower` and `sigma_upper`, if input. (Default: 3.0)"
    )

    parser.add_argument(
        "--sigma_lower",
        dest="sigma_lower",
        required=False,
        type=float,
        default=None,
        help="The number of standard deviations to use as the lower bound for the clipping limit. (Default: None)"
    )

    parser.add_argument(
        "--sigma_upper",
        dest="sigma_upper",
        required=False,
        type=float,
        default=None,
        help="The number of standard deviations to use as the upper bound for the clipping limit. (Default: None)"
    )

    parser.add_argument(
        "--maxiters",
        dest="maxiters",
        required=False,
        type=int,
        default=10,
        help="The maximum number of sigma-clipping iterations to perform or None to clip until convergence is achieved (i.e., iterate until the last iteration clips nothing). (Default: 10)"
    )

    args = parser.parse_args()

    #  Load config files for a given telescope
    filenames = list_files(args.path_data)
    filenames = make_copy(filenames, args.path_data,
                          outputDir="gmadet_subBkg/")
    bkg_estimation(filenames, box=args.box, filter_size=args.filter_size,
                   bkg_estimator=args.bkg_estimator, sigma=args.sigma_clip,
                   sigma_lower=args.sigma_lower, sigma_upper=args.sigma_upper,
                   maxiters=args.maxiters,outLevel=2)

if __name__ == "__main__":
    main()
