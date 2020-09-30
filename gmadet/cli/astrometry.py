#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

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

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    path_gmadet = getpath()
    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        description="Perform astrometric calibration of astronomical images."
    )

    parser.add_argument(
        "--path_data",
        dest="path_data",
        required=True,
        type=str,
        help="Path to file"
    )

    parser.add_argument(
        "--soft",
        dest="soft",
        required=False,
        choices=["scamp", "astrometry.net"],
        default="scamp",
        type=str,
        help="Soft to use for performing astrometric solution."
             "Not working with astrometry.net. (Default: scamp)",
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
    #  Create list of the copy images
    filenames = make_copy(filenames, args.path_data,
                          outputDir="gmadet_astrometry/")

    for filename in filenames:

        print("Sanitise header and data of %s.\n" % filename)
        sanitise_fits(filename)

        astrometric_calib(
            filename,
            config,
            soft=args.soft,
            verbose=args.verbose,
            accuracy=args.accuracy,
            itermax=args.itermax,
        )

if __name__ == "__main__":
    main()
