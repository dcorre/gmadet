#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

from gmadet.cnn.convert import convert

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        description="Convert simulated data before starting training."
    )

    parser.add_argument(
        "--path",
        dest="path_datacube",
        required=True,
        type=str,
        help="Path where to store the datacube."
    )

    parser.add_argument(
        "--cube",
        dest="cubename",
        required=True,
        type=str,
        help="Name of the datacube.",
    )

    parser.add_argument(
        "--cutouts",
        dest="path_cutouts",
        required=True,
        type=str,
        help="Path to the cutouts used for the training.",
    )

    parser.add_argument(
        "--frac_true",
        dest="frac_true",
        required=True,
        type=float,
        help="Fraction of true candidates to be included in the training set. Flat between 0 and 1.",
    )

    args = parser.parse_args()
    convert(args.path_datacube, cubename=args.cubename,
            path_cutouts=args.path_cutouts,
            frac_true=args.frac_true)


if __name__ == "__main__":
    main()
