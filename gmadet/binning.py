#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

"""
Rebin an astronomical image by summing all the pixels
for the required binning.

Usage example for a 2x2 binning:
    python binning.py --filename /home/corre/codes/gmadet/gmadet/data/ --binning 2 2

If need to be corrected for RN and expressed in ADU:
    python binning.py --filename /home/corre/codes/gmadet/gmadet/data/ --binning 2 2 --RN 10 --gain 1.8

As a result a rebinned_image/ folder will be created contataining the binned images. 
"""


import errno, glob, shutil, os
import numpy as np
from astropy.io import fits
import argparse

def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise


def rebin(arr, binning):
    """
    Rebinned an array with the requested binning. 
    The binned pixel is the sum of the orignal pixel values.
    """
    new_shape = np.array([arr.shape[0] / binning[0],
                          arr.shape[1]/binning[1]],
                          dtype=int )
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).sum(-1).sum(1)


def rebin_images(filename, binning, RN, gain):
    """
    Rebin images. Assumes that data are in te first hdu.

    Parameters
    ----------
    filename : path to image, string
        The file to read, with its extension. For ex: '/home/image.fits'
        If a directory is given, loop through all the fits file it contains
    binning : list
        define the biining factor along the 2 axis.
    RN: float or None
        Readout noise (rms) in electrons. Used to correct RN in binned pixels.
    gain: float or None
        gain (e-/ADU). Used to convert e- to ADU when correcting RN in binned pixels.

    Returns
    -------
    No variable is returned. The rebinned images are saved in a rebinned/ folder

    """

    # List all the files in the given path
    if os.path.isdir(filename):
        # expected extensions: everything having '.fit' in the file name
        filenames = glob.glob(filename + '*.fit*')
    else:
        filenames = [filename]

    # Create new folder where rebinned are saved
    path, filename_ext = os.path.split(filename) 

    if path:
        folder = path + '/'
    else:
        folder = ''

    resultDir = folder + 'rebinned_images/'
    # Create results folder
    mkdir_p(resultDir)

    # Loop through all files
    for ima in filenames:
        _, filename2 = os.path.split(ima) 

        hdul = fits.open(ima)
        hdul[0].data = rebin(hdul[0].data, binning)

        # If Readout noise is provided, correct the binned pixels.
        # Useful if one try to rebin an image (RN was applied to each pixels)
        # to make some comparisons/calculations with other images that were
        # initially binned (without applying RN to each individual pixel,
        # but only to the binned pixel).
        # Example: 2x2 binning. The resulting binned pixel has 4 times the RN,
        # Need to substract 3 x RN and divides by gain if data are expressed
        # in ADU.
        if RN:
            tot_pixels = binning[0] * binning[1]
            if gain is None:
                print ('No gain was provided. Set to 1.')
                gain = 1
            hdul[0].data = hdul[0].data - RN * (tot_pixels-1) / gain

        # Modify wcs information
        # Indexes for axis 1 and 2 might be switched between the fits information
        # and astropy. If so simply switch binning[0] and [1] below.
        hdul[0].header['Xbinning'] = binning[0]
        hdul[0].header['Ybinning'] = binning[1]
        hdul[0].header['CRPIX1'] = hdul[0].header['CRPIX1'] / binning[0]
        hdul[0].header['CRPIX2'] = hdul[0].header['CRPIX2'] / binning[1]
        # The following keywords are not necessarily present so use try/except
        try:
            hdul[0].header['CDELT1'] = hdul[0].header['CDELT1'] * binning[0]
            hdul[0].header['CDELT2'] = hdul[0].header['CDELT2'] * binning[1]
        except:
            pass
        try:
            hdul[0].header['CD1_1'] = hdul[0].header['CD1_1'] * binning[0]
            hdul[0].header['CD1_2'] = hdul[0].header['CD1_2'] * binning[0]
            hdul[0].header['CD2_2'] = hdul[0].header['CD2_2'] * binning[1]
            hdul[0].header['CD2_1'] = hdul[0].header['CD2_1'] * binning[1]
        except:
            pass

        hdul.writeto(resultDir+filename2,overwrite=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Rebin astronomical images.')

    parser.add_argument('--filename',
                        dest='filename',
                        required=True,
                        type=str,
                        help='Path to files. Do not forget a `/` at the end of the path if it is a directory.')

    parser.add_argument('--binning',
                       dest='binning',
                       nargs='+',
                       type=float,
                       default=[2, 2],
                       help='The binning factor to apply to each axis of the input images. Default: [2, 2]')

    parser.add_argument('--RN',
                        dest='RN',
                        required=False,
                        type=float,
                        help='Readout noise (rms) in electrons. If not passed as an argument, binned data are unchanged.')

    parser.add_argument('--gain',
                        dest='gain',
                        required=False,
                        type=float,
                        help='Gain (e-/ADU). Used to convert RN to substract in ADU, if data are expressed in ADU. If not do not pass it as an argument.')

    args = parser.parse_args()

    if args.RN:
        RN = args.RN
    else:
        RN = None
    if args.gain:
        gain = args.gain
    else:
        gain = None

    rebin_images(args.filename, args.binning, RN, gain)

