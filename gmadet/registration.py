#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import shutil
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from copy import deepcopy
from gmadet.utils import rm_p, mkdir_p
from gmadet.astrometry import scamp


def registration(filelist, config, resultDir="", reference=None,
                 useweight=False,gain=1, normalise_exp=True,
                 verbose="NORMAL"):
    """Register images"""
    #  Initialise lists used for creating output astropy table
    inim_list = []
    refim_list = []
    mask_list = []
    XY_lim = []
    in_lo = []
    in_up = []
    ref_lo = []
    ref_up = []
    gain_in = []
    gain_ref = []

    for i in range(len(filelist)):
        inim = filelist[i][0]
        refim = filelist[i][1]
        refim_mask = filelist[i][2]

        # Create list of images to register
        files = [inim, refim]
        if refim_mask is not None:
            files = files + [refim_mask]
        # Save them in a file to give it as an argument to swarp
        np.savetxt("register.list", files, fmt="%s")

        # Get pixel scale from input image header
        header = fits.getheader(inim)
        pixScale = abs(float(header["CDELT1"])) * 3600

        imalists = ["@" + "register.list"]
        # File name to store the common header that will be shared by all
        # images in filelist
        point = "registration"
        # Delete if already exists
        rm_p(point + ".head")
        # First run swarp to create a .head file containing the shared header
        subprocess.call(
            [
                "swarp",
                "-HEADER_ONLY", "Y",
                "-IMAGEOUT_NAME", point + ".head",
                "-GAIN_DEFAULT", str(gain),
                # '-VERBOSE_TYPE', verbose] + imalists)
                "-VERBOSE_TYPE", verbose,
            ]
            #+ [inim]
            + [refim]
        )
        # Some keywords manipulation using sed
        subprocess.call(
            [
                "sed",
                "-i",
                "s/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /",
                point + ".head",
            ]
        )
        outFiles = []
        # Run swarp to perform the registration on each image in filelist
        for j, ima in enumerate(files):

            path, filename_ext = os.path.split(ima)
            epoch = resultDir + os.path.splitext(filename_ext)[0] + \
                    "_reg_%s" % i
            outFiles.append(epoch + ".fits")

            if "mask" in ima:
                subBackground = "N"
            else:
                subBackground = "Y"

            # use weight for PS1 image
            if reference == 'ps1' and 'rings_v3_skycell' in ima and \
                    'mask' not in ima:
                #weight_type = "MAP_WEIGHT"
                #weight_type = "MAP_VARIANCE"
                weight_type = "MAP_RMS"
                #weight_type = "NONE"

                weight_name = path + '/' + \
                        os.path.splitext(filename_ext)[0] + \
                        ".weight.fits"
            else:
                weight_type = "NONE"
                weight_name = path + '/' + \
                        os.path.splitext(filename_ext)[0] + \
                        ".weight.fits"

            # Copy the common header in the .head file
            # So that it is read by sawrp for each image
            shutil.copy(point + ".head", epoch + ".head")
            
            # Use bilinear to avoid artefact, but worst for
            # noise so would need to check in more details.
            subprocess.call(
                [
                    "swarp",
                    "-IMAGEOUT_NAME", epoch + ".fits",
                    "-WEIGHT_TYPE", weight_type,
                    "-WEIGHT_IMAGE", weight_name,
                    "-WEIGHTOUT_NAME", '.weight.fits',
                    # Arbitrary threshold. 
                    # Pixels at the edsge after resampling are 0 so
                    # it is enough here
                    "-WEIGHT_THRESH", "0.1",
                    "-RESCALE_WEIGHTS", "N",
                    # '-GAIN_DEFAULT', str(gain),
                    "-FSCALE_KEYWORD", "NONE",
                    "-FSCALE_DEFAULT", "1, 1",
                    "-SUBTRACT_BACK", subBackground,
                    "-COMBINE", "Y",
                    "-COMBINE_TYPE", "MEDIAN",
                    "-BACK_SIZE", "128",
                    "-BACK_FILTERSIZE", "3",
                    "-RESAMPLE", "Y",
                    "-RESAMPLE_DIR", resultDir,
                    "-RESAMPLE_SUFFIX", '_test.fits',
                    "-PIXELSCALE_TYPE", "MANUAL",
                    "-PIXEL_SCALE", str(pixScale),
                    # '-CENTER', '%s, %s' % (header['CRVAL1'],
                    #                         header['CRVAL2']),
                    # '-RESAMPLING_TYPE', 'LANCZOS3',
                    "-RESAMPLING_TYPE", "BILINEAR",
                    # '-RESAMPLING_TYPE', 'NEAREST',
                    "-OVERSAMPLING", "0",
                    "-VERBOSE_TYPE", verbose,
                    "-COPY_KEYWORDS", "FILTER",
                ]
                + [ima]
            )

            # replace borders with NaNs in ref image if there are
            # any that are == 0,
            # hdulist=fits.open(epoch + '.fits')
            # hdulist[0].data[hdulist[0].data==0]=np.nan
            # hdulist.writeto(epoch + '.fits',overwrite=True)
            
            rm_p(epoch + ".head")
        rm_p(point + ".head")
        rm_p("register.list")
        rm_p("coadd.weight.fits")
        rm_p('swarp.xml')

        inim_regist = outFiles[0]
        refim_regist = outFiles[1]
        maskim_regist = outFiles[2]

        # Rescale flux to 1s
        # In the future can try to rescale flux of ref image 
        # to match input image.
        if normalise_exp:
            rescale_flux(inim_regist)
            rescale_flux(refim_regist)

        # Set masked pixels to same value
        mask_pix = flag_bad_pixels(inim_regist,
                                   mask_ref=maskim_regist,
                                   value=1e-30)
        # Update mask map
        _ = flag_bad_pixels(maskim_regist, value=1e8, mask_map=mask_pix)
        # Apply mask on ref data
        _ = flag_bad_pixels(refim_regist, mask_ref=maskim_regist, value=1e-30)

        print ('Remove bad pixels on the edge.')
        # Take only part of image with data
        # This will decrease image size and speed up the substraction
        # Might also avoid probelm with masked values, depending on how 
        # good they are deal with in hotpants.
        limits = keep_useful_area(inim_regist, image_ref=refim_regist)    
        _ = keep_useful_area(maskim_regist, limits_force=limits)

        # Perform a second time, as the edge are not straight, we can still
        # remove some pixels after the first cut.
        limits = keep_useful_area(inim_regist, image_ref=refim_regist)
        _ = keep_useful_area(maskim_regist, limits_force=limits)

        # Get info to tune hotpants parameters
        filelist_regist = [inim_regist, refim_regist, maskim_regist]
        hotpants_info = get_hotpants_info(filelist_regist, config, verbose)

        inim_list.append(outFiles[0])
        refim_list.append(outFiles[1])
        mask_list.append(outFiles[2])
        XY_lim.append(hotpants_info[0])
        in_lo.append(hotpants_info[1][0])
        in_up.append(hotpants_info[1][1])
        ref_lo.append(hotpants_info[1][2])
        ref_up.append(hotpants_info[1][3])
        gain_in.append(hotpants_info[2][0])
        gain_ref.append(hotpants_info[2][1])

    info = Table(
        [
            inim_list,
            refim_list,
            mask_list,
            XY_lim,
            in_lo,
            in_up,
            ref_lo,
            ref_up,
            gain_in,
            gain_ref,
        ],
        names=[
            "inim",
            "refim",
            "mask",
            "XY_lim",
            "in_lo",
            "in_up",
            "ref_lo",
            "ref_up",
            "gain_in",
            "gain_ref",
        ],
    )
    return info

def flag_bad_pixels(image, mask_ref=None, value=1e-30, mask_map=None):
    """
    Set masked pixels to same value, to homogeneise.
    1e-30 is the value used by hotpants to ignore pixels.
    """
    hdulist1 = fits.open(image)

    if mask_map is not None:
        # Apply an additional mask map to data
        hdulist1[0].data[mask_map] = value

    if mask_ref is not None:
        hdulist2 = fits.open(mask_ref)
        mask = hdulist2[0].data != 0.0
        hdulist2.close()

    if value == 1e-30:
        # Used to flag bad pixels in science and ref images by hotpants
        # Can actually be around the 1e-30 if resampled have been
        # done iwith weight on science image.
        hdulist1[0].data[hdulist1[0].data == 0] = value
        # Flag bad pixels from the mask
        if mask_ref is not None:
            hdulist1[0].data[mask] = value
    else:
        # For the mask map, pixels with high values are ignored.
        # 0 means they are not masked, i.e. good.
        hdulist1[0].data[hdulist1[0].data != 0] = value

    output_mask = hdulist1[0].data == value
    hdulist1.writeto(image, overwrite=True)

    return output_mask


def rescale_flux(image):
    """
    Rescale flux scale to 1s.
    In the future can also compute the flux factor scaling to match science 
    and input image flux scale.
    """

    # Normalise pixel flux to 1s.
    hdulist = fits.open(image)
    hdulist[0].data = hdulist[0].data / hdulist[0].header["EXPTIME"]
    # Try to update headers that are affected by this change.
    # RN and DC standing for Read Noise and Dark Current are currently
    # not kept during sanitising of the header. Need to update it.
    keywords=['SATURATE', 'RN', 'DC']
    for key in keywords:
        try:
            hdulist[0].header[key] = (
                hdulist[0].header[key] /
                hdulist[0].header["EXPTIME"]
            )
        except BaseException:
            pass
    # Set Exposure time to 1s from now on.
    hdulist[0].header["EXPTIME"] = 1
    hdulist.writeto(image, overwrite=True)


def keep_useful_area(image, image_ref=None, limits_force=None):
    """Keep only part of image with non masked pixels"""
    hdulist1 = fits.open(image)
    if image_ref is not None:
        hdulist2 = fits.open(image_ref)

    if limits_force is None:
        ymin_im, xmin_im = np.min(np.where(hdulist1[0].data > 1e-30), axis=1)
        ymax_im, xmax_im = np.max(np.where(hdulist1[0].data > 1e-30), axis=1)

        if image_ref is None:
            # Do it on a single image
            xmin = xmin_im
            xmax = xmax_im
            ymin = ymin_im
            ymax = ymax_im

        else:
            ymin_ref, xmin_ref = np.min(np.where(hdulist2[0].data > 1e-30), axis=1)
            ymax_ref, xmax_ref = np.max(np.where(hdulist2[0].data > 1e-30), axis=1)
           
            xmin = np.max([xmin_ref, xmin_im])
            xmax = np.min([xmax_ref, xmax_im])
            ymin = np.max([ymin_ref, ymin_im])
            ymax = np.min([ymax_ref, ymax_im])

    else:
        xmin, xmax, ymin, ymax = limits_force
    hdulist1[0].data = hdulist1[0].data[ymin:ymax, xmin:xmax]

    # Need to update center position info in header
    # Not CRVAL ?? Need to check
    xcenter1, ycenter1 = (hdulist1[0].header["CRPIX1"],
                          hdulist1[0].header["CRPIX2"])

    xcenter, ycenter = int(xcenter1 - xmin) + 1, int(ycenter1 - ymin) + 1
    hdulist1[0].header["CRPIX1"] = xcenter
    hdulist1[0].header["CRPIX2"] = ycenter

    hdulist1.writeto(image, overwrite=True)

    if image_ref is not None:
        hdulist2[0].data = hdulist2[0].data[ymin:ymax, xmin:xmax]

        # Need to update center position info in header
        xcenter1, ycenter1 = (hdulist2[0].header["CRPIX1"],
                              hdulist2[0].header["CRPIX2"])

        xcenter, ycenter = int(xcenter1 - xmin) + 1, int(ycenter1 - ymin) + 1
        hdulist2[0].header["CRPIX1"] = xcenter
        hdulist2[0].header["CRPIX2"] = ycenter

        hdulist2.writeto(image_ref, overwrite=True)

    return [xmin, xmax, ymin, ymax]


def get_hotpants_info(filelist, config, verbose):
    """Get some information for tuning hotpants
    """
    # Get the min and max pix with non 0 values
    # to delimate the new frame
    imdata, imheader = fits.getdata(filelist[0], header=True)
    refdata, refheader = fits.getdata(filelist[1], header=True)

    inmin = -10  # np.nanmin(imdata)
    if "SATURATE" in imheader:
        inmax = 0.9 * imheader["SATURATE"]
    else:
        inmax = 0.9 * np.nanmax(imdata)

    refmin = -10  # np.nanmin(refdata)
    if "SATURATE" in refheader:
        refmax = 0.9 * refheader["SATURATE"]
    else:
        refmax = 0.9 * np.nanmax(refdata)

    imgain = imheader.get('GAIN', 1)
    refgain = refheader.get('GAIN', 1)

    if imgain == 0 or imgain > 10:
        imgain = 1.0

    if refgain == 0 or refgain > 10:
        refgain = 1.0

    # Not needed anymore since using the weight output of swarp.
    # That flags correctly the pixels at the edge after the resampling.
    xmin = None
    ymin = None
    xmax = None
    ymax = None
    
    return [[xmin, xmax, ymin, ymax], [
        inmin, inmax, refmin, refmax], [imgain, refgain]]


if __name__ == "__main__":

    im1 = "Test_Atlascow/gmadet_results/ATLAS18qqn-S001-R001-C001-SDSS_g.fits"
    im2 = "ps1_rescaled/ps1_mosaic.fits"
    registration(im1, im2, False)
