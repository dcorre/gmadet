#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import shutil, subprocess
import numpy as np
from astropy.io import fits
from gmadet.astrometry import scamp
from gmadet.utils import rm_p

def create_mosaic(file_list, inputimage, outputDir, outName,
                  config, useweight=False, verbose='NORMAL',
                  astrometry=False):
    """Create a mosaic of images"""

    np.savetxt('mosaic.list', file_list, fmt='%s')

    mosaic_name = outputDir + outName 

    # Get pixel scale from input image header
    header = fits.getheader(inputimage)

    try:
        pixScale = abs(header['CDELT1'])
    except Exception:
        try:
            pixScale = abs(header['CD1_1'])
        except Exception:
            print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1 or CD1_1')
    pixScale = pixScale * 3600
    #print (inputimage, pixScale)
    crval1 = header['CRVAL1']
    crval2 = header['CRVAL2']
    #print (crval1, crval2)
    #print (header['CRPIX1'], header['CRPIX2'])
    imagesize = [header['NAXIS1'], header['NAXIS2']]

    # Force reference pixel to be in the center

    # File name to store the common header that will be shared by all
    # images in filelist
    point = 'registration'
    # Delete if already exists
    rm_p(point + '.head')
    # First run swarp to create a .head file containing the shared header
    subprocess.call(['swarp', '-HEADER_ONLY', 'Y', '-IMAGEOUT_NAME', \
                        point + '.head' , '-VERBOSE_TYPE', verbose] + [inputimage])
    # Some keywords manipulation using sed
    subprocess.call(['sed', '-i', \
                             's/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /', \
                     point + '.head'])

    imalists=['@' + 'mosaic.list']
    # Remove mosaic if already exists
    rm_p(mosaic_name + '.fits')
    if 'mask' in mosaic_name:
        subBackground = 'N'
    else:
        subBackground = 'Y'

    # Copy the common header in the .head file
    # So that it is read by sawrp for each image
    shutil.copy(point + '.head', mosaic_name + '.head')

    if useweight:
        subprocess.call(['swarp',
                         '-IMAGEOUT_NAME', mosaic_name + '.fits', \
                         '-WEIGHTOUT_NAME', mosaic_name + '.weight.fits', \
                         '-VERBOSE_TYPE', verbose] + imalists)
    else:
        subprocess.call(['swarp',
                         '-IMAGEOUT_NAME', mosaic_name + '.fits',\
                         '-SUBTRACT_BACK', subBackground, \
                         '-COMBINE', 'Y', \
                         '-BACK_SIZE', '128', \
                         '-BACK_FILTERSIZE', '3',\
                         #'-CENTER_TYPE', 'MANUAL', \
                         #'-CENTER', '%s, %s' % (crval1,crval2), \
                         '-RESAMPLE', 'Y',\
                         '-RESAMPLING_TYPE', 'LANCZOS3',\
                         #'-RESAMPLING_TYPE', 'BILINEAR',\
                         '-PIXELSCALE_TYPE', 'MANUAL', \
                         '-PIXEL_SCALE', str(pixScale), \
                         #'-IMAGE_SIZE', '%s, %s' % (imagesize[0], imagesize[1]), \
                         '-OVERSAMPLING', '0',\
                         '-COMBINE_TYPE', 'MEDIAN', \
                         '-COPY_KEYWORDS', ' PIXEL_SCALE', \
                         '-VERBOSE_TYPE', verbose] + imalists)
    rm_p(mosaic_name+'.head')
    # Perform astrometric calibration of the mosaic with scamp
    if astrometry:
        scamp(mosaic_name + '.fits', config, useweight=False, CheckPlot=False, verbose=verbose)

    rm_p('mosaic.list')
    rm_p('swarp.xml')
    rm_p(point+'.head')

    return True


