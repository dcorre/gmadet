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
from gmadet.psfex import psfex

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    path_gmadet = getpath()
    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        usage="usage: %(prog)s data [data2 ... dataN] [options]",
        description="Compute PSF of astronomical images."
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
        "--telescope",
        dest="telescope",
        choices=telescope_list,
        required=True,
        type=str,
        help="Alias for the available telescopes.",
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
        "--use-weight",
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

    args, filenames = parser.parse_known_args()

    # Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)

    filenames, subdirs = list_files(filenames, exclude=args.path_results)

    for raw_filename, subdir in zip(filenames, subdirs):
        filename = make_results_dir(
            raw_filename,
            outputDir=os.path.join(args.path_results, subdir),
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

        psfex(filename, config, useweight=args.useweight,
              verbose=args.verbose, outLevel=2)


if __name__ == "__main__":
    main()
