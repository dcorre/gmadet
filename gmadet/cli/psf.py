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
from gmadet.psfex import psfex

warnings.simplefilter(action="ignore", category=FutureWarning)

def main():

    path_gmadet = getpath()
    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        description="Compute PSF of astronomical images."
    )

    parser.add_argument(
        "--path_data",
        dest="path_data",
        required=True,
        type=str,
        help="Path to file"
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

    args = parser.parse_args()

    #  Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)
    filenames = list_files(args.path_data)
    filenames = make_copy(filenames, args.path_data,
                          outputDir="gmadet_psf/")
    psfex(filenames, config, useweight=args.useweight,
          verbose=args.verbose, outLevel=2)

if __name__ == "__main__":
    main()
