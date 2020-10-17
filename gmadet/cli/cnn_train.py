#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

"""

import argparse
import warnings

from gmadet.cnn.train import train

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():

    parser = argparse.ArgumentParser(
        description="Train the CNN using the given datacube."
    )

    parser.add_argument(
        "--cube",
        dest="path_cubename",
        required=True,
        type=str,
        help="Path to the datacube, including its name and extension.",
    )

    parser.add_argument(
        "--model-path",
        dest="path_model",
        required=True,
        type=str,
        help="Path where to store the trained model.",
    )

    parser.add_argument(
        "--model-name",
        dest="modelname",
        required=True,
        type=str,
        help="Name of the trained model.",
    )

    parser.add_argument(
        "--epochs",
        dest="epochs",
        required=False,
        type=int,
        default=10,
        help="Nunmber of epochs. (Default: 10)",
    )

    parser.add_argument(
        "--frac",
        dest="frac",
        required=False,
        type=float,
        default=0.15,
        help="Fraction of the data used for the validation sample. "
             "(Default: 0.15)",
    )

    parser.add_argument(
        "--dropout",
        dest="dropout",
        required=False,
        type=float,
        default=0.3,
        help="Fraction used for each dropout. "
             "(Default: 0.3)",
    )

    args = parser.parse_args()
    train(args.path_cubename, args.path_model, args.modelname,
          args.epochs, frac=args.frac, dropout=args.dropout)


if __name__ == "__main__":
    main()
