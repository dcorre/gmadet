#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import errno
import glob
import os
import shutil
import subprocess
import sys
import importlib
import time
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii, fits
from astropy.table import Table
from astropy import wcs
from astropy.wcs import WCS
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import (
    MinMaxInterval,
    SqrtStretch,
    LogStretch,
    SinhStretch,
    LinearStretch,
    ImageNormalize,
    ZScaleInterval,
)
from copy import deepcopy


def cp_p(src, dest):
    try:
        shutil.copy(src, dest)
    except BaseException:
        pass


def mv_p(src, dest):
    try:
        shutil.move(src, dest)
    except BaseException:
        pass


def rm_p(src):
    fileList = glob.glob(src, recursive=False)
    for filePath in fileList:
        try:
            os.remove(filePath)
        except BaseException:
            pass


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def getpath():
    """Get the path to gmadet module"""
    try:
        findspec = importlib.util.find_spec("gmadet")
        path = findspec.submodule_search_locations[0]
    except BaseException:
        print("path to gmadet can not be found.")

    return path


def getTel():
    """Get the list of all telescopes"""
    path_gmadet = getpath()
    telList = [
        name
        for name in os.listdir(os.path.join(path_gmadet, "config"))
        if os.path.isdir(os.path.join(path_gmadet, "config", name))
        and name != "conv_kernels"
    ]
    return telList


def is_subdir(path, basepath):
    """ Checks whether the path is inside basepath """
    path = os.path.abspath(path)
    basepath = os.path.abspath(basepath)

    return os.path.commonpath([path, basepath]) == basepath


def is_psf(filename, patterns=["_psf.fits"]):
    """
    Check whether a file is a PSF file with a standard name, typically
    '_psf.fits'. This is used for the simulation of sources during the CNN
    training for instance.
    """
    flag = False
    for ptn in patterns:
        if ptn in filename:
            flag = True
            break

    return flag


def list_files(
    paths,
    pattern=["*.fit", "*.fits", "*.fts"],
    recursive=True,
    get_subdirs=True,
    exclude=None,
):
    """(Recursively) list the files matching the pattern from the list of
    filenames or directories, omitting the ones containing file paths
    specified by 'exclude' option (either string or list of strings)"""

    filenames = []
    subdirs = []

    # Make sure paths is array-like
    paths = np.atleast_1d(paths)

    # Check if paths or files provided by user do exist.
    # filenames or folders can be a list
    for path in paths:
        if not os.path.exists(path):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

    # If there are than one pth provided keep the exact same structure
    # Otherwise skip the basepath
    for path in paths:
        # List all the files in the given path
        if os.path.isdir(path):
            # Recursively get all files matching given pattern(s)
            filenames_path = []

            for ptn in np.atleast_1d(pattern):
                filenames_path += glob.glob(path + "/**/" + ptn, recursive=recursive)
            # Sort alphanumerically
            filenames_path.sort()
            filenames += filenames_path
            # Get the relative path folder namesrelative to input paths
            reldirs = [
                os.path.dirname(os.path.relpath(_, path)) for _ in filenames_path
            ]
            # When more than one argument, add the base directory
            # being the last directory name provided by the user.
            if len(paths) > 1:
                basedir = os.path.basename(os.path.abspath(path))
                subdirs += [os.path.join(basedir, _) for _ in reldirs]
            else:
                subdirs += reldirs
        else:
            # if path is not a directory, assume it is a file
            filenames += [path]
            subdirs += [""]

    # Do not keep files in previous gmadet results folder
    # folder2skip = ["gmadet_results", "gmadet_astrometry",
    #                "gmadet_stacking", "gmadet_subBkg",
    #                "gmadet_psf", "gmadet_remove_cosmics",
    #                "candidates"]

    if isinstance(exclude, str):
        folder2skip = [exclude]
    elif exclude:
        folder2skip = exclude
    else:
        folder2skip = []

    idx = []  # Boolean mask whether to keep the filename or not

    for f in filenames:
        flag_2keep = True
        # If file is a standard PSF fits file
        if is_psf(f):
            flag_2keep = False
        else:
            for text in folder2skip:
                if is_subdir(f, text):
                    flag_2keep = False
                    break

        idx.append(flag_2keep)

    idx = np.array(idx)
    filenames = np.array(filenames)
    subdirs = np.array(subdirs)

    if get_subdirs:
        return filenames[idx], subdirs[idx]
    else:
        return filenames[idx]


def load_config(telescope, convFilter="default"):
    """Load the path to the configuration files required by the softs.
    They are telescope dependent.
    """
    path = getpath()
    path2tel = os.path.join(path, "config", telescope)
    config = {
        "telescope": telescope,
        "sextractor": {
            "conf": os.path.join(path2tel, "sourcesdet.sex"),
            "param": os.path.join(path2tel, "sourcesdet.param"),
            "convFilter": os.path.join(
                path, "config/conv_kernels/%s.conv" % convFilter
            ),
        },
        "scamp": {
            "sextractor": os.path.join(path2tel, "prepscamp.sex"),
            "param": os.path.join(path2tel, "prepscamp.param"),
            "conf": os.path.join(path2tel, "scamp.conf"),
        },
        "swarp": {},
        "psfex": {
            "sextractor": os.path.join(path2tel, "preppsfex.sex"),
            "param": os.path.join(path2tel, "preppsfex.param"),
            "conf": os.path.join(path2tel, "psfex.conf"),
        },
        "hotpants": {
            "conf": os.path.join(path2tel, "hotpants.hjson"),
            "conf2": os.path.join(path2tel, "hotpants_2.hjson"),
            "conf3": os.path.join(path2tel, "hotpants_3.hjson"),
        },
    }

    return config


def clean_folder(filelist, subFiles=None):
    """ Remove output files from previous iraf run. No need for sextractor  """

    types = ("*coo.*", "*mag.*", "*.magwcs", "*.magfiltered*")
    files2delete = []
    for filename in filelist:
        path = os.path.dirname(filename)

        for f in types:
            files2delete.extend(glob.glob(os.path.join(path, f)))
            files2delete.extend(glob.glob(f))

    files2delete = np.unique(files2delete)
    for f in files2delete:
        os.remove(f)


def make_results_dir(
    filename, outputDir="gmadet_results", keep=False, skip=False, copy=True
):
    """Copy original image to newly created output subdir inside given dir,
    optionally making backup of subfolder if it exists already"""

    # Filename without dir
    basename = os.path.basename(filename)

    # Subdir to store results for this file
    dirname = os.path.splitext(basename)[0]
    # TODO: optionally implement old behaviour?..
    dirname = os.path.join(outputDir, dirname)

    # Full path for the file
    newname = os.path.join(dirname, basename)

    # Make a backup of output dir if necessary
    if os.path.isdir(dirname):
        if skip:
            return None
        elif keep:
            mv_p(dirname, dirname + "_" + time.strftime("%Y%m%d-%H%M%S"))
        else:
            shutil.rmtree(dirname)

    mkdir_p(dirname)

    if copy:
        cp_p(filename, newname)

    return newname


def cut_image(filename, config, Nb_cuts=(2, 2), doAstrometry="scamp"):

    path, filename_ext = os.path.split(filename)
    filename2, extension = os.path.splitext(filename_ext)

    quadrant_list = []
    quadrant_ID = []

    if Nb_cuts == (1, 1):
        quadrant_list.append(filename)
        quadrant_ID.append("None")
    else:
        print("\nCutting %s into %d quadrants.\n" % (filename, np.sum(Nb_cuts)))
        if doAstrometry == "scamp":
            # Perform astrometry before cutting image into quadrants
            # It will ensure the astrometic calibration for each quadrant will
            # run smoothly
            from astrometry import scamp

            scamp(filename, config, accuracy=0.5, itermax=3, verbose="QUIET")

        data, header = fits.getdata(filename, header=True)
        w = wcs.WCS(header)
        Naxis1 = header["NAXIS1"]
        Naxis2 = header["NAXIS2"]
        Naxis11 = int(Naxis1 / Nb_cuts[0])
        Naxis22 = int(Naxis2 / Nb_cuts[1])

        index = 0
        for i in range(Nb_cuts[0]):
            for j in range(Nb_cuts[1]):
                index += 1
                if i == 0:
                    x1 = 1
                else:
                    x1 = Naxis11 * i + 1
                if j == 0:
                    y1 = 1
                else:
                    y1 = Naxis22 * j + 1
                x2 = Naxis11 * (i + 1)
                y2 = Naxis22 * (j + 1)

                filename_out = os.path.join(
                    path, filename2 + "_Q%d" % (index) + extension
                )

                # No need to update the header if astrometric calibration is
                # performed with scamp, this will be updated later.
                # datacut = data[x1-1:x2-1,y1-1:y2-1]
                datacut = data[y1 - 1 : y2 - 1, x1 - 1 : x2 - 1]
                newheader = deepcopy(header)
                """
                #Set center center of quadrant as CRPIX1,2
                #And compute the RA, Dec at this position for CRVAL1,2
                coeff_astro = ['PV', 'PC']
                #keys2delete=['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
                keys2delete=['CD1_2', 'CD2_1']
                for  key, value in newheader.items():
                   for coeff in coeff_astro:
                       _len = len(coeff)
                       if key[:_len] in [coeff]:
                           keys2delete.append(key)

                for keyword in keys2delete:
                    if keyword in newheader:
                        del newheader[keyword]
                crpix2 = int((y2-y1)/2)
                crpix1 = int((x2-x1)/2)
                ra, dec = w.wcs_pix2world(crpix1+Naxis11*i,
                                          crpix2+Naxis22*j,
                                          1)
                newheader['CRVAL1'] = float(ra)
                newheader['CRVAL2'] = float(dec)
                """
                # Use astrometric solution of the whole image by only
                # transforming the coordinate of CVAL1,2 for each quadrant
                newheader["CRPIX1"] = header["CRPIX1"] - Naxis11 * i
                newheader["CRPIX2"] = header["CRPIX2"] - Naxis22 * j

                fits.writeto(filename_out, datacut, newheader, overwrite=True)
                quadrant_list.append(filename_out)
                quadrant_ID.append("Q%d_%d_%d" % (index, i, j))

    image_table = Table([quadrant_list, quadrant_ID], names=["filenames", "quadrant"])

    return image_table


def make_sub_image(
    filenames,
    OT_coords,
    coords_type="world",
    sizes=[200, 200],
    FoVs=-1,
):
    """
    Extract sub-image around OT coordinates for the given size.

    Parameters
    ----------
    filename : path to image, string
        The file to read, with its extension. For ex: '/home/image.fits.gz'
    OT_coords : OT coordinates, list
        Coordinates of the OT, for instance [129.23, 45.27]. This coordinates
        are used as the center of the sub-image to create.
    coords_type : string, optional
        Either 'world' or 'pix'. 'world' means that coordinates are ra, dec
        expected in degrees format. 'pix' is the physical pixel coordinate
        on the detector, for instance [1248,2057].
        Default: 'world'
    size : list, optional
        define the size in pixels of the new sub-image.
        Default: [200,200]
    FoV: float, optional
        define the FoV in arcsec for the subimage. If -1 then the size is
        defined by `size`
    Returns
    -------
    subimage: array, float
              array of the data sub-image
    origin: string
            orientation of the figure

    """
    # Ensure all inputs are 1D
    filenames = np.atleast_1d(filenames)
    OT_coords = np.atleast_1d(OT_coords)
    coords_type = np.atleast_1d(coords_type)
    sizes = np.atleast_1d(sizes)
    FoVs = np.atleast_1d(FoVs)

    # Otherwise there is a problem using zip()
    if OT_coords.shape == (2,):
        OT_coords = [OT_coords]
    if sizes.shape == (2,):
        sizes = [sizes]

    subimages = []
    headers = []
    size_list = []
    pixref = []
    origin = []

    for fname, coords, _type, size, FoV in zip(
        filenames, OT_coords, coords_type, sizes, FoVs
    ):
        # Load file
        #data, header = fits.getdata(fname, header=True)
        hdul = fits.open(fname)
        data = hdul[0].data
        header = hdul[0].header
        hdul.close()
        headers.append(header)
        # Get physical coordinates of OT
        w = WCS(header)
        if _type == "world":
            # Get physical coordinates
            c = coord.SkyCoord(coords[0], coords[1], unit=(u.deg, u.deg), frame="icrs")
            world = np.array([[c.ra.deg, c.dec.deg]])
            # print (world)
            pix1, pix2 = w.all_world2pix(world, 1)[0]
            pix = [pix2, pix1]
            pixref.append(coords)
        elif _type == "pix":
            pix = coords
            # print (pix)
            # ra, dec = w.all_pix2world(np.array(pix), 0)
            ra, dec = w.all_pix2world(pix[0], pix[1], 0)
            pixref.append([float(ra), float(dec)])

        if FoV > 0:
            # Get pixel size in degrees
            try:
                pixSize = abs(float(header["CDELT1"]))
            except BaseException:
                pixSize = abs(float(header["CD1_1"]))
            # Compute number of pixels to reach desired FoV in arcseconds
            size = [int(FoV / (pixSize * 3600)), int(FoV / (pixSize * 3600))]

        size_list.append(size)
        # Extract subimage from image starting from reference pixel
        x1 = int(pix[0]) - int(size[0] / 2)
        if x1 < 0:
            x1 = 0
        x2 = int(pix[0]) + int(size[0] / 2)
        y1 = int(pix[1]) - int(size[1] / 2)
        if y1 < 0:
            y1 = 0
        y2 = int(pix[1]) + int(size[1] / 2)
        subimages.append(data[x1:x2, y1:y2])

        # Highest declination on top
        ra1, dec1 = w.all_pix2world(pix[0], y1, 0)
        ra2, dec2 = w.all_pix2world(pix[0], y2, 0)
        if dec1 > dec2:
            origin.append("upper")
        else:
            origin.append("lower")
    return [subimages, headers, size_list, pixref, origin]


def make_fits(data, output_name, header, size, pixref, info_dict=None):
    """Create a fits file from a data array"""
    hdu = fits.PrimaryHDU()
    hdu.data = data.astype(np.float32)
    # FIXME: need to adapt header here !!! Likely not correct
    header["CRPIX1"] = int(size[0] / 2)
    header["CRPIX2"] = int(size[1] / 2)
    header["CRVAL1"] = pixref[0]
    header["CRVAL2"] = pixref[1]
    if info_dict is not None:
        # Add information regarding the transients
        for key, value in info_dict.items():
            header[key] = value
        if [data.shape[0], data.shape[1]] != size:
            header["edge"] = "True"
        else:
            header["edge"] = "False"
    hdu.header = header
    hdu.writeto(output_name, overwrite=True)


def make_figure(data, output_name, origin, fmt, title=None):
    """Make a figure from a data array"""

    norm = ImageNormalize(
        # subimage - np.median(subimage),
        data,
        interval=ZScaleInterval(),
        stretch=LinearStretch(),
    )

    plt.figure()
    plt.imshow(data, cmap="gray", origin=origin, norm=norm)
    # plt.imshow(norm(subimage),cmap="gray", origin=origin)
    if title is not None:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(output_name, format=fmt)
    plt.close()

    # Much faster but without annotations
    # plt.imsave(output_name, norm(subimage),
    #            cmap="gray", origin=origin, format=fmt)
    # plt.close()


def combine_cutouts(
    filenames,
    OT_coords,
    coords_type="world",
    output_name="cutout_comb.png",
    size=[200, 200],
    FoV=-1,
    title=None,
):
    """Create a png file with the cutouts from science image,
    reference image and substarcted image."""

    # FIXME: need to optimise this function
    # May be provide the list of cutouts to create as inputs
    # to reuse the plt axes and avoid recreating a new pyplot
    # frame each time.
    data1, _, _, _, origin1 = make_sub_image(
        filenames[0],
        OT_coords,
        coords_type,
        size,
        FoV,
    )
    data2, _, _, _, origin2 = make_sub_image(
        filenames[1],
        OT_coords,
        coords_type,
        size,
        FoV,
    )
    data3, _, _, _, origin3 = make_sub_image(
        filenames[2],
        OT_coords,
        coords_type,
        size,
        FoV,
    )
    # Outputs of make_sub_image are lists
    data1 = data1[0]
    data2 = data2[0]
    data3 = data3[0]
    origin1 = origin1[0]
    origin2 = origin2[0]
    origin3 = origin3[0]

    norm1 = ImageNormalize(
        data1,  # - np.median(data1),
        interval=ZScaleInterval(),
        stretch=LinearStretch(),
    )
    norm2 = ImageNormalize(
        data2,  # - np.median(data2),
        interval=ZScaleInterval(),
        stretch=LinearStretch(),
    )
    norm3 = ImageNormalize(
        data3,  # - np.median(data3),
        interval=ZScaleInterval(),
        stretch=LinearStretch(),
    )

    # stretch=SinhStretch())
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    # Substracting by the median is a trick to highlight the source
    # when skybackground is important. It is not correct but just
    # used for illustration.
    # axs[1].imshow(data1 - np.median(data1),
    axs[1].imshow(data1, cmap="gray", origin=origin1, norm=norm1)
    axs[1].set_xlabel("Science", size=20)
    # axs[2].imshow(data2 - np.median(data2),
    axs[2].imshow(data2, cmap="gray", origin=origin2, norm=norm2)
    axs[2].set_xlabel("Reference", size=20)
    # axs[0].imshow(data3 #- np.median(data3),
    axs[0].imshow(data3, cmap="gray", origin=origin3, norm=norm3)
    axs[0].set_xlabel("Residuals", size=20)
    if title is not None:
        fig.suptitle(title, size=20)
    # Tight_layout() does not support suptitle so need to do it manually.
    fig.tight_layout(rect=[0, 0.03, 1, 0.80])
    fig.savefig(output_name)
    plt.close()


def get_corner_coords(filename):
    """Get the image coordinates of an image"""

    header = fits.getheader(filename)
    # Get physical coordinates
    Naxis1 = header["NAXIS1"]
    Naxis2 = header["NAXIS2"]

    pix_coords = [[0, 0, Naxis1, Naxis1], [0, Naxis2, Naxis2, 0]]

    # Get physical coordinates
    w = WCS(header)
    ra, dec = w.all_pix2world(pix_coords[0], pix_coords[1], 1)

    return [ra, dec]


def get_phot_cat(filename, telescope):
    """Get the name of the filter band from header and telescope name
    And associate the correct name from DB
    """
    header = fits.getheader(filename)

    # FILTER keyword required
    try:
        band = header["FILTER"]
    except Exception:
        print("No FILTER keyword found in header.")

    if band in ["C", "Clear", "NoFilter", "L"]:
        band_DB = "Clear"
        band_cat = "g+r"
    elif band in ["g", "gSDSS", "SDSS g"]:
        band_DB = "g/AB"
        band_cat = "g"
    elif band in ["r", "rSDSS", "rPATH", "SDSS r"]:
        band_DB = "r/AB"
        band_cat = "r"
    elif band in ["i", "iSDSS", "SDSS i"]:
        band_DB = "i/AB"
        band_cat = "i"
    elif band in ["z", "zSDSS", "SDSS z"]:
        band_DB = "z/AB"
        band_cat = "z"
    elif band in ["B"]:
        band_DB = "B/Johnson"
        band_cat = "B"
    elif band in ["V"]:
        band_DB = "V/Johnson"
        band_cat = "V"
    elif band in ["R"]:
        band_DB = "R/Johnson"
        band_cat = "R"
    elif band in ["I"]:
        band_DB = "I/Johnson"
        band_cat = "I"

    # Chose photometric catalog
    try:
        RA = header["CRVAL1"]
    except Exception:
        print("No CRVAL1 keyword found header. " "Astrometric calibration required.")
    try:
        DEC = header["CRVAL2"]
    except Exception:
        print("No CRVAL2 keyword found header. " "Astrometric calibration required.")

    # Use Pan-Starrs if Dec > -30 degrees
    if float(DEC) > -30.0:
        catalog = "II/349/ps1"
    # Else use SDSS if available.
    elif Vizier.query_region(
        field, width=rad_deg, height=rad_deg, catalog="V/147/sdss12"
    )[0]:
        catalog = "V/147/sdss12"
    # Else use Gaia, but no calibration available for z band.
    elif filter not in z_band:
        catalog = "I/345/gaia2"
    # Last choice: USNO-B1. All-sky but bad photometric calibration.
    else:
        catalog = "I/284/out"

    return band_DB, band_cat, catalog


def unpackbits(x, num_bits):
    """Unpack bits with any dimension ndarray.
    Can unpack however many bits"""
    xshape = list(x.shape)
    x = x.reshape([-1, 1])
    to_and = 2 ** np.arange(num_bits).reshape([1, num_bits])
    return (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])


def filter_catalog_data(data, catalogName):
    """Remove extended sources and bad measurements from reference catalogs
    before performing photometric calibration"""

    # Keep only point source objects and good measurements

    # Panstarrs flags
    if catalogName == "II/349/ps1":
        # First remove data using the general 'Qual' flag.
        # ------------------------------------------------
        # There are 8 bits
        # Bit 1: extended object in PS1
        # Bit 2: Extended in external data (e.g. 2MASS)
        # Bit 3: Good-quality measurement in PS1
        # Bit 4: Good-quality measurement in external data (eg, 2MASS)
        # Bit 5: Good-quality object in the stack (>1 good stack measurement)
        # Bit 6: The primary stack measurements are the best measurements
        # Bit 7: Suspect object in the stack (no more than 1 good measurement,
        #        2 or more suspect or good stack measurement)
        # Bit 8:  Poor-quality stack object (no more than 1 good or suspect
        # measurement)

        Quality_flags = np.array(data["Qual"])
        Qual_flag = unpackbits(Quality_flags, 8)

        qual_bits = [1, 2, 3, 7, 8]
        qual_values = [0, 0, 1, 0, 0]
        counter = 0
        for i, j in zip(qual_bits, qual_values):
            condition = Qual_flag[:, i - 1] == j
            if counter == 0:
                quality_mask = condition
            else:
                quality_mask = np.bitwise_and(quality_mask, condition)
            counter += 1
        # flag for individual bands
        # -------------------------
        # There are 25 bits. Use only the bit stating whether it is
        # an extended object in this band.
        # Bit 1: Used within relphot (SECF_STAR_FEW): skip star
        # Bit 2: Used within relphot (SECF_STAR_POOR): skip star
        # Bit 3: Synthetic photometry used in average measurement
        # Bit 4: Ubercal photometry used in average measurement
        # Bit 5: PS1 photometry used in average measurement
        # Bit 6: PS1 stack photometry exists
        # Bit 7: Tycho photometry used for synthetic magnitudes
        # Bit 8: Synthetic magnitudes repaired with zeropoint map
        # Bit 9: Average magnitude calculated in 0th pass
        # Bit 10: Average magnitude calculated in 1th pass
        # Bit 11: Average magnitude calculated in 2th pass
        # Bit 12: Average magnitude calculated in 3th pass
        # Bit 13: Average magnitude calculated in 4th pass
        # Bit 14: Extended in this band (PSPS only)
        # Bit 15: PS1 stack photometry comes from primary skycell
        # Bit 16: PS1 stack best measurement is a detection (not forced)
        # Bit 17: PS1 stack primary measurement is a detection (not forced)
        # Bit 18:
        # Bit 19:
        # Bit 20:
        # Bit 21: This photcode has SDSS photometry
        # Bit 22: This photcode has HSC photometry
        # Bit 23: This photcode has CFH photometry (mostly megacam)
        # Bit 24: This photcode has DES photometry
        # Bit 25: Extended in this band

        band_bits_and = [1, 2, 14, 25]
        band_values_and = [0, 0, 0, 0]
        band_bits_or = [9, 10, 11, 12]
        band_values_or = [1, 1, 1, 1]

        # bands = ['g', 'r', 'i', 'z', 'y']
        # No need to consider y band, and fainter sensitivity, so
        # might remove good reference stars
        bands = ["g", "r", "i", "z"]
        band_flags = []
        # Unpack bits from individual band flags
        for band in bands:
            _temp = np.array(data["%sFlags" % band])
            band_flags.append(unpackbits(_temp, 25))
        band_flags = np.array(band_flags)

        # Apply mask conditions
        for i in range(len(bands)):
            counter = 0
            for j1, k1 in zip(band_bits_and, band_values_and):
                condition = band_flags[i][:, j1 - 1] == k1
                quality_mask = np.bitwise_and(quality_mask, condition)
                counter += 1
            counter2 = 0
            # At least one Average magnitude calculated
            for j2, k2 in zip(band_bits_or, band_values_or):
                condition_or = band_flags[i][:, j2 - 1] == k2
                if counter2 == 0:
                    quality_mask_or = condition_or
                else:
                    quality_mask_or = np.bitwise_or(quality_mask_or, condition_or)
                counter2 += 1
            # Combine both masks
            quality_mask = np.bitwise_and(quality_mask, quality_mask_or)

    elif catalogName == "V/147/sdss12":
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    elif catalogName == "I/345/gaia2":
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    elif catalogName == "I/284/out":
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    return data[quality_mask]


def clean_outputs(filenames, outLevel):
    """Delete non required output files"""

    imagelist = np.atleast_1d(filenames)
    for ima in imagelist:
        print("\nCleaning up output files for %s" % ima)
        path = os.path.dirname(ima)

        rm_p(os.path.join(path, "*.head"))
        if outLevel == 0:
            rm_p("coadd.weight.fits")
            for _ in [
                "*.magwcs*",
                "*.oc*",
                "*.cat",
                "*.png",
                "*.fits",
                "*_ZP_*",
                "substraction/*.magwcs*",
                "substraction/*.cat",
                "substraction/*mask*",
                "substraction/*background",
                "substraction/*segmentation",
            ]:
                rm_p(os.path.join(path, _))

        elif outLevel == 1:
            rm_p("coadd.weight.fits")

            for _ in [
                "*.magwcs*",
                "*.oc*",
                "*.cat",
                "*.png",
                "*.fits",
                "*_ZP_*",
                "*_CRmask.fits",
                "*_CR_notcleaned.fits",
                "substraction/*.magwcs*",
                "substraction/*mask*",
                "substraction/*background",
                "substraction/*segmentation",
            ]:
                rm_p(os.path.join(path, _))

        elif outLevel == 2:
            pass
