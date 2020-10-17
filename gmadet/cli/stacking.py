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
        usage="usage: %(prog)s [options] data",
        description="Stack astronomical images."
    )

    parser.add_argument(
        "--results",
        dest="path_results",
        required=False,
        type=str,
        default='gmadet_stacking',
        help="Base path to store the results. "
             "(Default: gmadet_stacking)"
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
        "--radius",
        dest="radius",
        required=False,
        type=float,
        default=5,
        help="Radius in arcmin to group fields. (Default: 5 arcmin)",
    )

    parser.add_argument(
        "--deltat",
        dest="deltaT",
        required=False,
        type=float,
        default=1,
        help="Time interval in hours to group fields into same epoch. "
             "(Default: 1h)",
    )

    parser.add_argument(
        "--no-backsub",
        dest="no_BackSub",
        action="store_false",
        help="If provided as argument, no background substraction is "
             "performed on each image prior to stacking. "
             "Default: Substraction is performed",
    )

    args, paths = parser.parse_known_args()

    for path in paths:
        stacking(path, args.radius, args.deltaT,
                 subBack=args.no_BackSub, path_results=args.path_results, keep=args.keep)

if __name__ == "__main__":
    main()
