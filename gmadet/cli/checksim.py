#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

from gmadet.cnn.checksim import makestats

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        usage="usage: %(prog)s [options] data",
        description="Do some tests to quantify the simulation."
    )

    parser.add_argument(
        "--radius",
        dest="radius",
        required=False,
        default=2,
        type=float,
        help="Radius for crossmatching detected sources with simulated events."
             " Default: 2 arcseconds",
    )

    args, paths = parser.parse_known_args()

    for path in paths:
        makestats(path, radius=args.radius)


if __name__ == "__main__":
    main()
