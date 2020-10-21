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
    getpath,
    getTel
)
from gmadet.cnn.makesubimage import subimage

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    path_gmadet = getpath()
    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        usage="usage: %(prog)s [options] data",
        description="Create cutouts centered on the optical candidates."
    )

    parser.add_argument(
        "--training",
        dest="training",
        action="store_true",
        help="If set, enters training mode. (Default: normal mode)"
    )

    parser.add_argument(
        "--false",
        dest="false",
        action="store_true",
        help="If set, put all unmatched objects to false folder "
             "in training mode. (Default: not applied)"
    )

    parser.add_argument(
        "--size",
        dest="size",
        required=False,
        default=32,
        type=int,
        help="Size in pixels of the extracted images. (Default: 32)",
    )

    parser.add_argument(
        "--radius",
        dest="radius",
        required=False,
        default=2,
        type=float,
        help="Radius for crossmatching detected sources with simulated events."
             " (Default: 2 arcseconds)",
    )

    parser.add_argument(
        "--flag-notsub",
        dest="flag_notsub",
        required=False,
        action="store_true",
        help="Whether the candidates are not the results of an image "
             "substraction. (Default: False)",
    )

    args, paths = parser.parse_known_args()

    for path in paths:
        subimage(
            path,
            args.training,
            size=args.size,
            radius=args.radius,
            flag_notsub=args.flag_notsub,
            false=args.false)


if __name__ == "__main__":
    main()
