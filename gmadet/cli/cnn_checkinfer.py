#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

from gmadet.cnn.checkinfer import makestats

warnings.simplefilter(action="ignore", category=FutureWarning)

def main():

    parser = argparse.ArgumentParser(
        description="Do some tests to quantify the CNN training."
    )

    parser.add_argument(
        "--path_plots",
        dest="path_plots",
        required=True,
        type=str,
        help="Path where to store plots."
    )

    parser.add_argument(
        "--path_crossmatch",
        dest="path_crossmatch",
        required=True,
        type=str,
        help="Path where crossmatch.dat is stored."
    )

    parser.add_argument(
        "--path_infer",
        dest="path_infer",
        required=True,
        type=str,
        help="Path where infer.dat is stored."
    )

    parser.add_argument(
        "--maglim",
        dest="maglim",
        required=False,
        type=float,
        nargs='+',
        default=[12,18,21],
        help="Magnitudes used to split the magnitude range in the plots. "
             " (Default: 12 18 21)"
    )

    parser.add_argument(
        "--CNNproblim",
        dest="CNNproblim",
        required=False,
        type=float,
        nargs='+',
        default=[0.1,0.5,0.7],
        help="CNN proba used to split the results in the plots. "
             " (Default: 0.1 0.5 0.7)"
    )

    parser.add_argument(
        "--FWHM_ratio_lower",
        dest="FWHM_ratio_lower",
        required=False,
        type=float,
        default=0.5,
        help="Lower bound for the ratio FWHM / FWHM_PSF used for the plots. "
             " (Default: 0.5)"
    )
    
    parser.add_argument(
        "--FWHM_ratio_upper",
        dest="FWHM_ratio_upper",
        required=False,
        type=float,
        default=4.2,
        help="Upper bound for the ratio FWHM / FWHM_PSF used for the plots. "
             " (Default: 4.2)"
    )

    args = parser.parse_args()
    makestats(args.path_plots, args.path_crossmatch, args.path_infer,
              maglim=args.maglim, CNNproblim=args.CNNproblim,
              FWHM_ratio_lower=args.FWHM_ratio_lower,
              FWHM_ratio_upper=args.FWHM_ratio_upper)


if __name__ == "__main__":
    main()
