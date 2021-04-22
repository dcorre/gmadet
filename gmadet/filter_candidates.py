#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import join
from astropy.time import Time
from astropy.io import ascii, fits
from gmadet.cnn.infer import infer
from gmadet.utils import (
    make_sub_image,
    make_fits,
    make_figure,
    combine_cutouts,
    mkdir_p,
)
import multiprocessing as mp


def filter_candidates(
    sources,
    FWHM_ratio_lower=0.5,
    FWHM_ratio_upper=5.0,
    CNN_model=None,
    CNN_thres=0.0,
    makecutout=True,
    size=32,
    size_cnn=32,
    fmt="png",
    outLevel=1,
    nb_threads=8,
    combined=False,
):
    """Filter transient candidates"""

    # if no sources skip the following
    if len(sources) == 0:
        print("No candidates, no need to filter.")
        return 0

    print("Filter candidates")
    # Take first candidate to extract the path where to store the file
    # No need to chack if substraction was performed,
    # as if it did only the ones from substracted files are with 'Match' == Y
    path, fname_ext = os.path.split(sources["filenames"][0])
    # Get rid of the extension to keep only the name
    fname2, extension = os.path.splitext(fname_ext)
    # Get rid of the _reg pattern
    fname2 = fname2.split("_ref")[0]

    # First get the sources not crossmatching with sources in catalogs
    mask_cat = sources["Match"] == "N"

    # Remove candidates on the edges
    mask_edge = sources["edge"] == "N"

    # Remove sources with FWHM ratio outside the desired range
    FWHM_ratio = sources["FWHM"] / sources["FWHMPSF"]
    mask_FWHM = (FWHM_ratio >= FWHM_ratio_lower) & (FWHM_ratio <= FWHM_ratio_upper)

    mask_tot = mask_cat & mask_edge & mask_FWHM
    # Create dictionary with fits info for cutouts to be be given to the CNN model
    # or simply for making fits cutouts
    if CNN_model is not None or fmt == "fits":

        if CNN_model is not None:
            path_CNN_cutouts = os.path.join(path, "CNN_cutouts")
            mkdir_p(path_CNN_cutouts)
            outnames = []
        args_data = []
        info_dicts = []
        for cand in sources:
            coords = [cand["_RAJ2000"], cand["_DEJ2000"]]
            if CNN_model is not None:
                outname = os.path.join(
                    path_CNN_cutouts, "candidate_%d.fits" % (cand["idx"])
                )
                outnames.append(outname)
            info_dict = {}
            info_dict["RA"] = cand["_RAJ2000"]
            info_dict["DEC"] = cand["_DEJ2000"]
            info_dict["XPOS"] = cand["Xpos"]
            info_dict["YPOS"] = cand["Ypos"]
            info_dict["FILE"] = cand["filenames"]
            info_dict["CANDID"] = cand["idx"]
            info_dict["MAG"] = cand["mag_calib"]
            info_dict["MAGERR"] = cand["mag_calib_err"]
            info_dict["FWHM"] = cand["FWHM"]
            info_dict["FWHMPSF"] = cand["FWHMPSF"]

            args_data.append(
                [
                    cand["filenames"],
                    coords,
                    "world",
                    [size_cnn, size_cnn],
                    -1,
                ]
            )

            info_dicts.append(info_dict)

    # Use a trained CNN model to filter candidates.
    if CNN_model is not None:
        print("Create fits cutouts for CNN")

        # Create sub-array
        args_data = np.array(args_data)
        make_sub_image(
            args_data[:, 0],
            outnames,
            args_data[:, 1],
            args_data[:, 2],
            args_data[:, 3],
            args_data[:, 4],
            info_dicts,
            [None] * len(outnames),
            [fmt] * len(outnames),
            nb_threads,
        )

        # The size of the cutout should be the same as the ones used
        # for the CNN training
        print("Use trained CNN model")
        infer(path_CNN_cutouts, CNN_model, 0.1)

        # Add the probability to the canidates table.
        infer_table = ascii.read(os.path.join(path_CNN_cutouts, "infer_results.dat"))
        sources = join(
            sources, infer_table["cand_ID", "label0", "label1"], join_type="left"
        )

        # keep only transients that are above the threshold
        mask_CNN = sources["label1"] >= CNN_thres
        mask_tot = mask_tot & mask_CNN

    # Write output file.
    candidates = sources[mask_tot]

    if fmt == "fits":
        info_dicts_filtered = np.array(info_dicts)[mask_tot]

    # if no sources skip the following
    if len(candidates) == 0:
        print("No candidates, no need to filter.")
        return 0

    # Create ID to start from 1.
    candidates["cand_ID"] = np.arange(len(candidates)) + 1
    # Rename colums
    if "label0" in candidates.colnames:
        candidates.rename_column("label0", "P_False")
        candidates.rename_column("label1", "P_True")
    candidates.write(
        os.path.join(path, fname2 + "_candidates.dat"),
        format="ascii.commented_header",
        overwrite=True,
    )

    # Update
    print("Make cutouts")
    # Extract small image centered on candidates passing the filters
    if makecutout:
        path_cutout = os.path.join(path, "cutouts")
        mkdir_p(path_cutout)
        args_data = []
        outnames = []
        titles = []
        if combined:
            args_combined = []
            path_cutout_combined = os.path.join(path_cutout, "combined")
            mkdir_p(path_cutout_combined)

        for cand in candidates:
            coords = [cand["_RAJ2000"], cand["_DEJ2000"]]
            outname = os.path.join(
                path_cutout, "candidate_%d.%s" % (cand["cand_ID"], fmt)
            )
            if combined:
                outname_combined = os.path.join(
                    path_cutout_combined,
                    "candidate_%d_comb.%s" % (cand["cand_ID"], "png"),
                )

            header = fits.getheader(cand["OriginalIma"])
            try:
                date = Time(header["DATE-OBS"], format="fits")
                # convert in GPS time
                # date_JD = date.jd
            except BaseException:
                date = Time(header["MJD-OBS"], format="mjd")
            date.format = "iso"
            if fmt != "fits" or combined:
                _coords = SkyCoord(
                    cand["_RAJ2000"],
                    cand["_DEJ2000"],
                    unit=(u.degree, u.degree),
                    frame="icrs",
                )
                coords_sexa = _coords.to_string(style="hmsdms")
                title = (
                    "RA Dec: %s \n" % coords_sexa
                    + "Time (UTC): %s \n" % (date.value)
                    + "Mag: %.2f +/- %.2f     "
                    % (cand["mag_calib"], cand["mag_calib_err"])
                    + "     FWHM_ratio: %.2f" % (cand["FWHM"] / cand["FWHMPSF"])
                )
                if CNN_model is not None:
                    title += "     CNN proba: %.2f " % cand["P_True"]

                titles.append(title)

            args_data.append(
                [
                    cand["filenames"],
                    coords,
                    "world",
                    [size, size],
                    -1,
                ]
            )
            outnames.append(outname)

            if combined:
                args_combined.append(
                    [
                        [cand["OriginalIma"], cand["RefIma"], cand["filenames"]],
                        coords,
                        "world",
                        outname_combined,
                        [size, size],
                        -1,
                        title,
                        [fmt] * len(outname_combined),
                    ]
                )

        if fmt != "fits":
            info_dicts_filtered = [None] * len(outnames)
        else:
            titles = [None] * len(outnames)

        # Create sub-array
        args_data = np.array(args_data)

        make_sub_image(
            args_data[:, 0],
            outnames,
            args_data[:, 1],
            args_data[:, 2],
            args_data[:, 3],
            args_data[:, 4],
            info_dicts_filtered,
            titles,
            [fmt] * len(outnames),
            nb_threads,
        )

        if combined:
            print("Make combined cutouts")
            p = mp.Pool(nb_threads)
            p.starmap(combine_cutouts, args_combined)
            p.close()
