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
        "--path_datacube",
        dest="path_datacube",
        required=True,
        type=str,
        help="Path where to store the datacube."
    )

    parser.add_argument(
        "--cubename",
        dest="cubename",
        required=True,
        type=str,
        help="Name of the datacube.",
    )

    parser.add_argument(
        "--path_cutouts",
        dest="path_cutouts",
        required=True,
        type=str,
        help="Path to the cutouts used for the training.",
    )

    args = parser.parse_args()
    convert(args.path_datacube, cubename=args.cubename,
            path_cutouts=args.path_cutouts)

if __name__ == "__main__":
    main()
