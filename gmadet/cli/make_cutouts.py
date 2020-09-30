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
from gmadet.cnn.makesubimage import subimage

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    path_gmadet = getpath()
    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        description="Create cutouts centered on the optical candidates."
    )

    parser.add_argument(
        "--path_data",
        dest="path_data",
        required=True,
        type=str,
        help="Path to file"
    )

    parser.add_argument(
        "--training",
        dest="training",
        action="store_true",
        help="If set, enters training mode. (Default: normal mode)"
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
        "--flag_notsub",
        dest="flag_notsub",
        required=False,
        action="store_true",
        help="Whether the candidates are not the results of an image substraction."
             "(Default: False)",
    )


    args = parser.parse_args()
    subimage(args.path_data, args.training,
             size=args.size, radius=args.radius, flag_notsub=args.flag_notsub)

if __name__ == "__main__":
    main()
