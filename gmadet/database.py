#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Scripts to send data to a given database

"""

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.time import Time, TimeDelta
from shapely.geometry import Point, Polygon
import requests
import json
import voeventparse as vp
from gmadet.utils import make_sub_image, get_corner_coords


def send_data2DB(
    filename,
    candidates,
    Nb_cuts,
    owncloud_path,
    VOE_path,
    usrpwd_path,
    FoV=60,
    coords_type="wolrd",
    corner_cut=32,
    debug=False,
    fmt="png",
    subFiles=None,
):
    """Send candidates information to database"""

    # Load original image header to retrieve date of observation and airmass
    # These information might have been lost in the analyis images
    header = fits.getheader(filename)
    dateObs = header["DATE-OBS"]
    Tstart = Time(dateObs, format="fits", scale="utc")
    try:
        exposure = float(header["EXPOSURE"])
        Tend = Tstart + TimeDelta(exposure, format="sec")
    except BaseException:
        Tend = Time(dateObs, format="fits", scale="utc")
        exposure = Tend - Tstart
    # Try to get airmass rom header, else set it to -1
    try:
        Airmass = header["AIRMASS"]
    except BaseException:
        Airmass = -1

    #  Compute the image corner RA, Dec coordinates
    ra, dec = get_corner_coords(filename)
    pix_im_coord = np.array([ra, dec]).T
    im_poly = Polygon([tuple(co) for co in pix_im_coord])

    # Do not consider candidates found in the image edge
    # imsize = data.shape
    # print (imsize, header['NAXIS1'], header['NAXIS2'])
    # Get the physical pixels of the original size if image were split into
    # different quadrants.
    """
    for i, candidate in enumerate(candidates):
        quadrant_idx = candidate['quadrant']
        if quadrant_idx == 'None':
            quadrant = None
            index_i = 0
            index_j = 0
        else:
            quadrant, index_i, index_j = quadrant_idx.split('_')
            quadrant = quadrant[1:]

        candidates['Xpos'][i] = candidate['Xpos'] + (int(imsize[0]/Nb_cuts[0])
                                                     * int(index_j))
        candidates['Ypos'][i] = candidate['Ypos'] + (int(imsize[1]/Nb_cuts[1])
                                                     * int(index_i))
    #print (candidates)
    mask = (candidates['Xpos'] > corner_cut) & \
            (candidates['Ypos'] > corner_cut) & \
            (candidates['Xpos'] < imsize[1] - corner_cut) & \
            (candidates['Ypos'] < imsize[0] - corner_cut)

    candidates_cut = candidates[mask]
    """
    # Get information about the current alert from the xml file containing
    # observation plan
    with open(VOE_path, "rb") as f:
        obsplan = vp.load(f)

    dict_event = {}
    dict_event["event_type"] = obsplan.find(
        ".//Param[@name='Event_type']").attrib["value"]
    dict_event["event_name"] = obsplan.find(
        ".//Param[@name='Event_ID']").attrib["value"]
    dict_event["event_status"] = obsplan.find(
        ".//Param[@name='Event_status']").attrib["value"]
    dict_event["revision"] = obsplan.find(
        ".//Param[@name='Revision']").attrib["value"]
    dict_event["telescope"] = obsplan.find(
        ".//Param[@name='Name_tel']").attrib["value"]
    tiles_info = get_obsplan(obsplan)

    # Get user email adress and password to login in
    # https://grandma-fa-interface.lal.in2p3.fr
    with open(usrpwd_path) as f:
        usrpwd = json.load(f)

    # Set up the output repository path to store sub-images
    outputDir = (
        owncloud_path
        + "/"
        + dict_event["event_type"]
        + "/"
        + dict_event["event_name"]
        + "/"
        + dict_event["event_status"]
        + "_"
        + dict_event["revision"]
        + "/OTs/"
    )

    # Create a sub image centered on each candidate found, and gather
    # information
    # tile_id_list = [1] * len(candidates)
    # filter_list = candidates['filter_DB']
    # Tstart_list = [Tstart.fits] * len(candidates)
    # Tend_list = [Tend.fits] * len(candidates)
    # Airmass_list = [Airmass] * len(candidates)
    ImFits_path = []
    RefFits_path = []
    SubFits_path = []
    tile_id_list = []
    if subFiles:
        mask = candidates["FlagSub"] == "Y"
        candidates = candidates[mask]
    masktest = (
        (candidates["_RAJ2000"] < 244.01)
        & (candidates["_RAJ2000"] > 244.0)
        & (candidates["_DEJ2000"] < 22.27)
        & (candidates["_DEJ2000"] > 22.26)
    )
    candidates = candidates[masktest]
    for i, row in enumerate(candidates):
        name = (
            dict_event["telescope"]
            + "_"
            + str(round(float(row["_RAJ2000"]), 5))
            + "_"
            + str(round(float(row["_DEJ2000"]), 5))
            + "_"
            + dateObs
            + "."
            + fmt
        )

        OT_coords_wcs = [row["_RAJ2000"], row["_DEJ2000"]]
        OT_coords_pix = [row["Ypos"], row["Xpos"]]
        # Extract the region given wcs coordinates.
        # If substraction was performed, images are realigned and some
        # cut on the edges might performed, so the physical pixels position
        # do not match exactly the original image.
        # Astrometry is performed for the original image and sustracted image,
        # so no problem.
        make_sub_image(
            row["OriginalIma"],
            OT_coords_wcs,
            coords_type="world",
            output_name=outputDir + name,
            size=[128, 128],
            FoV=FoV,
            fmt=fmt,
        )
        ImFits_path.append(name)
        if subFiles:
            name_ref = (
                dict_event["telescope"]
                + "_"
                + str(round(float(row["_RAJ2000"]), 5))
                + "_"
                + str(round(float(row["_DEJ2000"]), 5))
                + "_"
                + dateObs
                + "_ref."
                + fmt
            )
            name_sub = (
                dict_event["telescope"]
                + "_"
                + str(round(float(row["_RAJ2000"]), 5))
                + "_"
                + str(round(float(row["_DEJ2000"]), 5))
                + "_"
                + dateObs
                + "_sub."
                + fmt
            )
            RefFits_path.append(name_ref)
            SubFits_path.append(name_sub)
            #  Here use physical pixels coordinates.
            #  Should work as well with wcs coordinates
            make_sub_image(
                row["RefIma"],
                OT_coords_wcs,
                coords_type="world",
                output_name=outputDir + name_ref,
                size=[128, 128],
                FoV=FoV,
                fmt=fmt,
            )
            make_sub_image(
                row["filenames"],
                OT_coords_wcs,
                coords_type="world",
                output_name=outputDir + name_sub,
                size=[128, 128],
                FoV=FoV,
                fmt=fmt,
            )

        else:
            RefFits_path.append("")
            SubFits_path.append("")
        #  Check in which tile the OT is located
        #  Check that the RA, Dec lies within the FoV rectangle
        #  Stop when finding one. Assuming that tiles are not overlapping.
        tile_id = 0  #  by default
        for tile in tiles_info:
            tile_center = Point(tiles_info["RA"], tiles_info["Dec"])
            if tile_center.intersects(im_poly):
                tile_id = tiles_info["Id"]
                break
        tile_id_list.append(tile_id)

    alias = ["new"] * len(candidates)
    new = [1] * len(candidates)
    RA_list = candidates["_RAJ2000"]
    Dec_list = candidates["_DEJ2000"]
    filter_list = candidates["filter_DB"]
    Tstart_list = [Tstart.fits] * len(candidates)
    Tend_list = [Tend.fits] * len(candidates)
    exp_list = [exposure] * len(candidates)
    Mag_list = candidates["mag_calib"]
    Mag_err_list = candidates["mag_calib_err"]
    Magsys_list = candidates["magsys"]
    Airmass_list = [Airmass] * len(candidates)

    candidates_2DB = Table(
        [
            alias,
            new,
            tile_id_list,
            RA_list,
            Dec_list,
            filter_list,
            Tstart_list,
            Tend_list,
            exp_list,
            Mag_list,
            Mag_err_list,
            Magsys_list,
            Airmass_list,
            ImFits_path,
            RefFits_path,
            SubFits_path,
        ],
        names=[
            "alias",
            "new",
            "tile_id",
            "RA",
            "DEC",
            "filter",
            "Tstart",
            "Tend",
            "Exposure",
            "Magnitude",
            "Magnitude_error",
            "Magsys",
            "Airmass",
            "im_fits_name",
            "ref_fits_name",
            "sub_fits_name",
        ],
    )

    # Set url to report tile or galaxy observations
    # url = "https://grandma-fa-interface.lal.in2p3.fr/obs_report_OT.php"
    # url = "http://localhost/grandma/obs_report_OT.php"
    #  Loop over the observations
    for i in range(len(candidates_2DB)):
        data2DB = {}
        for col in candidates_2DB.colnames:
            data2DB[col] = candidates_2DB[col][i]

        # Add obsplan info to data dictionary
        for key, value in dict_event.items():
            data2DB[key] = value

        # Add username and password to data dictionary
        for key, value in usrpwd.items():
            data2DB[key] = value

        # Add compulsory keys
        data2DB["method"] = "POST"
        data2DB["submit"] = "ok"

        response = requests.post(url, data=data2DB)

        if response.status_code == 200:
            print("Data sent succesfully to database.")
            forced_debug = False
        else:
            print("Data not sent to database. See information below.")
            forced_debug = True

        if debug or forced_debug:
            print("\nDEBUG:\n")
            print("Data sent to DB:")
            print(data2DB)
            print("\n\n")
            print("Request response text:")
            print(response.text)
            print("\n\n")
            print("Request response status code:")
            print(response.status_code)
            print("\n\n")
            print("Request response history:")
            print(response.history)

            forced_debug = False


def get_obsplan(v):
    """Extract the tiles id and RA, Dec from voevent"""
    ID = []
    Ra = []
    Dec = []
    Grade = []
    Header = []
    for element in v.What.iterchildren():
        tag_type = str(element.tag)
        if tag_type == "Table":
            for subelement in element.iterchildren():
                tag_type = str(subelement.tag)
                if tag_type == "Field":
                    Header.append(str(subelement.attrib["name"]))
                    if tag_type == "Data":
                        for lines in subelement.iterchildren():
                            ID.append(int(lines.TD[0]))
                            Ra.append(float(lines.TD[1]))
                            Dec.append(float(lines.TD[2]))
                            Grade.append(float(lines.TD[3]))

    tiles_info = Table([ID, Ra, Dec, Grade], names=[
                       "Id", "RA", "Dec", "Grade"])

    return tiles_info
