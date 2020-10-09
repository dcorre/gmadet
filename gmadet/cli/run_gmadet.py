#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: David Corre, Orsay, France, corre@lal.in2p3.fr

Input arguments are filename, typical fwhm (or estimated by psfex),
sextracting threshold and maximal distance for catalog crosschecking in pixels

Example :
   python gmadet.py --filename /folder/image.fits --FWHM psfex
                    --threshold 4 --radius_crossmatch 2.5 --telescope TRE

"""

import sys
import subprocess
import glob
import math
import shutil
import os
import argparse
import warnings

from gmadet.phot_calibration import phot_calib
from gmadet.utils import (
    load_config,
    clean_folder,
    cut_image,
    list_files,
    mv_p,
    mkdir_p,
    make_copy,
    make_results_dir,
    clean_outputs,
    getpath,
    getTel
)
from gmadet.sanitise import sanitise_fits
from gmadet.remove_cosmics import run_lacosmic
from gmadet.astrometry import astrometric_calib
from gmadet.psfex import psfex
from gmadet.sources_extraction import (
    run_sextractor,
    filter_sources,
    convert_xy_radec
)
from gmadet.substraction import substraction
from gmadet.background import bkg_estimation
from gmadet.crossmatch import catalogs, moving_objects
from gmadet.database import send_data2DB
from gmadet.filter_candidates import filter_candidates

from astropy.io import ascii, fits
from astropy.table import vstack, Table, Column

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

from copy import deepcopy

warnings.simplefilter(action="ignore", category=FutureWarning)

def main():

    path_gmadet = getpath()

    telescope_list = getTel()

    parser = argparse.ArgumentParser(
        description="Finding unknown objects in astronomical images."
    )

    parser.add_argument(
        "--path_data",
        "--data",
        dest="path_data",
        action="append",
        required=False,
        type=str,
        help="Path to data"
    )

    parser.add_argument(
        "--path_results",
        "--results",
        dest="path_results",
        required=False,
        type=str,
        default='gmadet_results',
        help="Base path to store the results. "
             "(Default: gmadet_results)"
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
        "--skip-processed",
        "--skip",
        dest="skip",
        required=False,
        action="store_true",
        help="Skip already processed files"
    )

    parser.add_argument(
        "--FWHM",
        "--fwhm",
        dest="FWHM",
        required=True,
        help="Typical telescope FWHM"
    )

    parser.add_argument(
        "--radius_crossmatch",
        "--radius-crossmatch",
        dest="radius_crossmatch",
        required=False,
        type=float,
        default=3.0,
        help="Radius to use for crossmatching, in pixels. "
             "(Default: 3.0 pixels)",
    )

    parser.add_argument(
        "--threshold",
        dest="threshold",
        required=False,
        default=4.0,
        type=float,
        help="Consider only sources above this threshold. "
             "(Default: 4.0)",
    )

    parser.add_argument(
        "--soft",
        dest="soft",
        required=False,
        choices=["sextractor"],
        default="sextractor",
        type=str,
        help="Soft to use for detecting sources.\n (Default: sextractor)",
    )

    parser.add_argument(
        "--convFilter",
        "--conv-filter",
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
        "--telescope",
        dest="telescope",
        choices=telescope_list,
        required=True,
        type=str,
        help="Alias for the available telescopes.",
    )

    parser.add_argument(
        "--quadrants",
        dest="quadrants",
        required=False,
        default=1,
        type=int,
        help="Number of quadrants the image is divided. "
             "(Default: 1)",
    )

    parser.add_argument(
        "--doAstrometry",
        "--do-astrometry",
        dest="doAstrometry",
        required=False,
        default="scamp",
        choices=["No", "scamp"],
        type=str,
        help="Whether to perform astrometric calibration, with scamp. "
             "(Default: scamp)",
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

    parser.add_argument(
        "--doSub",
        "--do-sub",
        dest="doSub",
        required=False,
        type=str,
        help="Whether to perform astrometric calibration, with ps1 images "
             'or user provided reference image. Type "ps1" for PS1 reference '
             'image or provide the path to your reference image.',
    )

    parser.add_argument(
        "--ps1_method",
        "--ps1-method",
        dest="ps1_method",
        required=False,
        default="mosaic",
        choices=["mosaic", "individual"],
        type=str,
        help="When substracting images using Pan-STARRS reference images, "
             "there 2 options, either create a mosaic of all PS1 image and "
             "substract or do the substraction individually for each PS1 "
             "image. In the latter case, your image is cut to match the "
             "PS1 image. (Default: mosaic)",
    )

    parser.add_argument(
        "--doMosaic",
        "--do-mosaic",
        dest="doMosaic",
        action="store_true",
        help="Whether to combine the individual frames into a common mosaic "
             "when `ps1_method` is set to `individual`. (Default: not set)",
    )

    parser.add_argument(
        "--Remove_cosmics",
        "--remove-cosmics",
        dest="Remove_cosmics",
        action="store_true",
        help="Whether to remove cosmic rays using lacosmic. "
             " (Default: not set)",
    )

    parser.add_argument(
        "--sub_bkg",
        "--sub-bkg",
        dest="sub_bkg",
        action="store_true",
        help="Whether to substract background. (Default: not set)",
    )

    parser.add_argument(
        "--output_data_level",
        "--output-data-level",
        dest="outLevel",
        required=False,
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="Number of output files that are kept after the process. "
             "0: minimum, 2: maximum"
             "(Default: 2)",
    )

    parser.add_argument(
        "--owncloud_path",
        "--owncloud-path",
        dest="owncloud_path",
        required=False,
        type=str,
        help="Local path to the owncloud",
    )

    parser.add_argument(
        "--VOE_path",
        "--voe-path",
        dest="VOE_path",
        required=False,
        type=str,
        help="Path/filename of the VoEvent containing the observation plan.",
    )

    args,filenames = parser.parse_known_args()

    # Combine free-form arguments with path_data into single list of filenames
    if isinstance(args.path_data, str):
        filenames = [args.path_data] + filenames
    elif args.path_data:
        filenames = list(args.path_data) + filenames

    if not filenames:
        return

    Nb_cuts = (args.quadrants, args.quadrants)

    # Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)

    filenames = list_files(filenames, exclude=args.path_results)

    for raw_filename in filenames:
        # We should copy images to results dir one by one, while we are processing them
        filename = make_results_dir(raw_filename, outputDir=args.path_results, keep=args.keep, skip=args.skip)

        if not filename:
            print("%s is already processed, skipping. \n" % raw_filename)
            continue

        print("Sanitise header and data of %s.\n" % filename)
        sanitise_fits(filename)

        # Cut image into several quadrants if required
        # And create table with filename and quadrant ID
        image_table = cut_image(
            filename,
            config,
            Nb_cuts=Nb_cuts,
            doAstrometry=args.doAstrometry
        )

        if args.Remove_cosmics:
            print(
                "Running lacosmic on %s to remove cosmic rays. \n" %
                filename)
            # Clean cosmic rays
            # Not using FWHM anymore
            FWHM_list = [None] * len(image_table)
            run_lacosmic(
                image_table["filenames"],
                FWHM_list,
                contrast=5.0,
                cr_threshold=5.0,
                neighbor_threshold=5.0,
                maxiter=4,
                outLevel=args.outLevel
            )

        if args.sub_bkg:
            # Substract background
            bkg_estimation(
                image_table["filenames"],
                box=(20, 20),
                filter_size=(3, 3),
                bkg_estimator='SExtractor',
                sigma=3.,
                sigma_lower=None,
                sigma_upper=None,
                maxiters=10,
                outLevel=args.outLevel,
            )

        if args.FWHM == "psfex":
            # Estimate the PSF FWHM for each image/quadrants using psfex
            FWHM_list = psfex(
                image_table["filenames"],
                config,
                verbose=args.verbose,
                outLevel=args.outLevel,
            )
        else:
            FWHM_list = [args.FWHM] * len(image_table)

        if args.doAstrometry != "No":
            astrometric_calib(
                image_table["filenames"],
                config,
                soft=args.doAstrometry,
                verbose=args.verbose,
                accuracy=0.15,
                itermax=10
            )

        if args.doSub:
            substracted_files = substraction(
                image_table["filenames"],
                args.doSub,
                config,
                soft="hotpants",
                method=args.ps1_method,
                doMosaic=args.doMosaic,
                verbose=args.verbose,
                outLevel=args.outLevel,
            )
        else:
            substracted_files = None

        if args.soft == "sextractor":
            run_sextractor(
                image_table["filenames"],
                FWHM_list,
                args.threshold,
                args.telescope,
                config,
                verbose=args.verbose,
                subFiles=substracted_files,
                outLevel=args.outLevel,
            )

        filter_sources(
            image_table["filenames"],
            args.soft,
            sigma=1,
            subFiles=substracted_files,
        )
        convert_xy_radec(
            image_table["filenames"],
            soft=args.soft,
            subFiles=substracted_files
        )
        total_sources, candidates = catalogs(
            image_table,
            args.radius_crossmatch,
            Nb_cuts=Nb_cuts,
            subFiles=substracted_files,
        )
        #moving_objects(image_table["filenames"], candidates)

        sources_calib = phot_calib(
            total_sources,
            args.telescope,
            radius=args.radius_crossmatch,
            doPlot=True,
            subFiles=substracted_files,
        )
        # Apply filter to candidates
        # Remove candidates on the edge
        # Remove candidate depending the FWHM ratio
        # Apply the CNN model
        candidates_filtered = filter_candidates(
            sources_calib
        )


        # If both arguments VOE_path and owncloud_path are provided
        # Send candidates to database
        # Set the tile_id corresponding to your tile by hand at the moment
        if args.VOE_path and args.owncloud_path:
            send_data2DB(
                filename,
                candidates_filtered,
                Nb_cuts,
                args.owncloud_path,
                args.VOE_path,
                "utilsDB/usrpwd.json",
                debug=True,
                subFiles=substracted_files,
            )

        # clean output files
        clean_outputs(image_table["filenames"], args.outLevel)

if __name__ == "__main__":
    main()
