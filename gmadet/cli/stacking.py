#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

#from gmadet.utils import list_files
from gmadet.stacking import stacking

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        description="Stack astronomical images."
    )

    parser.add_argument(
        "--path_data",
        "--data",
        dest="path_data",
        required=True,
        type=str,
        help="Path where the files to be stacked are.",
    )

    parser.add_argument(
        "--radius",
        dest="radius",
        required=False,
        type=float,
        default=5,
        help="Radius in arcmin to group fields. (Default: 5 arcmin)",
    )

    parser.add_argument(
        "--deltaT",
        "--deltat",
        dest="deltaT",
        required=False,
        type=float,
        default=1,
        help="Time interval in hours to group fields into same epoch. "
             "(Default: 1h)",
    )

    parser.add_argument(
        "--no_BackSub",
        "--no-backsub",
        dest="no_BackSub",
        action="store_false",
        help="If provided as argument, no background substraction is "
             "performed on each image prior to stacking. "
             "Default: Substraction is performed",
    )

    args = parser.parse_args()
    stacking(args.path_data, args.radius, args.deltaT,
             subBack=args.no_BackSub)

if __name__ == "__main__":
    main()
