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
        dest="path_data",
        required=True,
        type=str,
        help="Path to data"
    )

    parser.add_argument(
        "--FWHM",
        dest="FWHM",
        required=True,
        help="Typical telescope FWHM"
    )

    parser.add_argument(
        "--radius_crossmatch",
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
        dest="doSub",
        required=False,
        type=str,
        help="Whether to perform astrometric calibration, with ps1 images "
             'or user provided reference image. Type "ps1" for PS1 reference '
             'image or provide the path to your reference image.',
    )

    parser.add_argument(
        "--ps1_method",
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
        "--Remove_cosmics",
        dest="Remove_cosmics",
        action="store_true",
        help="Whether to remove cosmic rays using lacosmic.",
    )

    parser.add_argument(
        "--sub_bkg",
        dest="sub_bkg",
        action="store_true",
        help="Whether to substract background.",
    )

    parser.add_argument(
        "--output_data_level",
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
        dest="owncloud_path",
        required=False,
        type=str,
        help="Local path to the owncloud",
    )

    parser.add_argument(
        "--VOE_path",
        dest="VOE_path",
        required=False,
        type=str,
        help="Path/filename of the VoEvent containing the observation plan.",
    )

    args = parser.parse_args()

    Nb_cuts = (args.quadrants, args.quadrants)

    #  Load config files for a given telescope
    config = load_config(args.telescope, args.convFilter)

    filenames = list_files(args.path_data)
    #  copy original images
    #  Create list of the copy images
    filenames = make_copy(filenames, args.path_data,
                          outputDir="gmadet_results/")

    for filename in filenames:

        print("Sanitise header and data of %s.\n" % filename)
        sanitise_fits(filename)

        #  Cut image into several quadrants if required
        #  And create table with filename and quadrant ID
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
            #  Clean cosmic rays
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
            #  Substract background
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
        total_candidates = catalogs(
            image_table,
            args.radius_crossmatch,
            Nb_cuts=Nb_cuts,
            subFiles=substracted_files,
        )
        # moving_objects(args.filename, total_candidates)

        # total_candidates = ascii.read('total_candidates.dat')
        total_candidates_calib = phot_calib(
            total_candidates,
            args.telescope,
            radius=args.radius_crossmatch,
            doPlot=True,
            subFiles=substracted_files,
        )

        # total_candidates_calib = ascii.read(
        #    'Test_sendDB/gmadet_results/jul1919-010r_sh_tot_cand2.dat')

        #  If both arguments VOE_path and owncloud_path are provided
        #  Send candidates to database
        # Set the tile_id corresponding to your tile by hand at the moment
        if args.VOE_path and args.owncloud_path:
            send_data2DB(
                filename,
                total_candidates_calib,
                Nb_cuts,
                args.owncloud_path,
                args.VOE_path,
                "utilsDB/usrpwd.json",
                debug=True,
                subFiles=substracted_files,
            )

        #  clean output files
        clean_outputs(image_table["filenames"], args.outLevel)

if __name__ == "__main__":
    main()
