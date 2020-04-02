#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, subprocess, shutil
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from utils import rm_p, mkdir_p 
from copy import deepcopy
from astrometry import scamp

def registration(inim, refim, config, refim_mask=None, useweight=False, gain=1, debug=False):
    """Register images"""

    # Create folder with substraction results
    path, filename_ext = os.path.split(inim)
    # Get rid of the extension to keep only the name
    filename2 = filename_ext.split('.')[0]

    if path:
        folder = path + '/'
    else:
        folder = ''

    resultDir = folder + 'substraction/'
    # Create results folder
    mkdir_p(resultDir)

    # Read header from input image
    hdulist = fits.open(inim)
    hdulist[0].data = hdulist[0].data / hdulist[0].header['EXPTIME']
    try:
        hdulist.header[0]['SATURATE'] = hdulist[0].header['SATURATE'] / hdulist[0].header['EXPTIME']
    except:
        pass
    hdulist[0].header['EXPTIME'] = 1
    hdulist.writeto(inim, overwrite=True)
    #header = fits.getheader(inim)

    # Create list of images to register
    filelist=[inim, refim]
    if refim_mask is not None:
        filelist = filelist + [refim_mask]
    # Save them in a file to give it as an argument to swarp
    np.savetxt('register.list', filelist, fmt='%s')

    imalists=['@' + 'register.list']
    # File name to store the common header that will be shared by all
    # images in filelist
    point = 'registration'
    # Delete if already exists
    rm_p(point + '.head')
    # First run swarp to create a .head file containing the shared header
    subprocess.call(['swarp', '-HEADER_ONLY', 'Y', '-IMAGEOUT_NAME', \
                        point + '.head', '-GAIN_DEFAULT', str(gain)] + imalists)
    # Some keywords manipulation using sed
    subprocess.call(['sed', '-i', \
                             's/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /', \
                     point + '.head'])
    # Run swarp to perform the registration on each image in filelist
    for i, ima in enumerate(filelist):
        if 'ps1' in ima:
            subBackground = 'N'
        else:
            subBackground = 'Y'
        path, filename_ext = os.path.split(ima)
        epoch = resultDir + filename_ext.split('.')[0] + '_regist'
        # Copy the common header in the .head file 
        # So that it is read by sawrp for each image
        shutil.copy(point + '.head',epoch + '.head')
        if useweight:
            subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', epoch + '.fits', \
                             '-WEIGHTOUT_NAME', epoch + '.weight.fits', \
                             '-GAIN_DEFAULT', str(gain)] + [ima])
        else:
            subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', epoch + '.fits',\
                             #'-GAIN_DEFAULT', str(gain),\
                             '-FSCALE_KEYWORD', 'NONE', \
                             '-FSCALE_DEFAULT', '1, 1', \
                             '-SUBTRACT_BACK', subBackground, \
                             '-COMBINE', 'Y', \
                             '-BACK_SIZE', '128', \
                             '-BACK_FILTERSIZE', '3',\
                             #'-CENTER', '%s, %s' % (header['CRVAL1'],header['CRVAL2']), \
                             #'-RESAMPLING_TYPE', 'LANCZOS3',\
                             '-RESAMPLING_TYPE', 'BILINEAR',\
                             #'-RESAMPLING_TYPE', 'NEAREST', \
                             '-OVERSAMPLING', '0',\
                             '-COMBINE_TYPE', 'MEDIAN', \
                             '-COPY_KEYWORDS', 'FILTER'] + [ima])

        # replace borders with NaNs in ref image if there are any that are == 0,
        #hdulist=fits.open(epoch + '.fits')
        #hdulist[0].data[hdulist[0].data==0]=np.nan
        #hdulist.writeto(epoch + '.fits',overwrite=True)

        rm_p(epoch + '.head')
    rm_p('*.head')
    rm_p('register.list')
    rm_p('coadd.weight.fits')

    path, filename_ext = os.path.split(inim)
    inim_regist = resultDir + filename_ext.split('.')[0] + '_regist.fits'
    path, filename_ext = os.path.split(refim)
    refim_regist = resultDir + filename_ext.split('.')[0] + '_regist.fits'
    path, filename_ext = os.path.split(refim_mask)
    maskim_regist = resultDir + filename_ext.split('.')[0] + '_regist.fits'

    # Set specific value for mask image
    hdulist = fits.open(inim_regist)
    # 1e-30 is the default value of bad pixels for hotpants
    hdulist[0].data[hdulist[0].data == 0] = 1e-30
    #hdulist[0].data[hdulist[0].data != 1e-30] /= 20
    hdulist.writeto(inim_regist, overwrite=True)

    hdulist = fits.open(refim_regist)
    # 1e-30 is the default value of bad pixels for hotpants
    hdulist[0].data[hdulist[0].data == 0] = 1e-30
    #hdulist[0].data[hdulist[0].data != 1e-30] /= 20
    hdulist.writeto(refim_regist, overwrite=True)

    hdulist = fits.open(maskim_regist)
    hdulist[0].data[hdulist[0].data != 0] = 1e8
    hdulist.writeto(maskim_regist, overwrite=True)

    filelist_regist = [inim_regist, refim_regist, maskim_regist]
    overlap_info = get_exact_overlap(filelist_regist, resultDir, config)

    return overlap_info

def get_exact_overlap(filelist, resultDiri, config):
    """Cut the registered images to get rid of the edges filled by 0 or nan"""

    # Get the min and max pix with non 0 values 
    # to delimate the new frame
    imdata, imheader = fits.getdata(filelist[0], header=True)
    refdata, refheader = fits.getdata(filelist[1], header=True)
    ymin1, xmin1 = np.min(np.where(imdata>1e-30),axis=1)
    ymax1, xmax1 = np.max(np.where(imdata>1e-30),axis=1)
    ymin2, xmin2 = np.min(np.where(refdata>1e-30),axis=1)
    ymax2, xmax2 = np.max(np.where(refdata>1e-30),axis=1)

    xmin=np.max([xmin1,xmin2])
    xmax=np.min([xmax1,xmax2])
    ymin=np.max([ymin1,ymin2])
    ymax=np.min([ymax1,ymax2])
    
    hdulist = fits.open(filelist[2])
    # 1e-30 is the default value of bad pixels for hotpants
    good_region = deepcopy(hdulist[0].data[xmin:xmax,ymin:ymax])
    hdulist[0].data[:] = 1e8
    hdulist[0].data[xmin:xmax,ymin:ymax] = good_region
    #hdulist[0].data[hdulist[0].data != 1e-30] /= 20
    hdulist.writeto(filelist[2], overwrite=True)
    
    inmin = -10#np.nanmin(imdata)
    try:
        inmax = imheader['SATURATE']
    except:
        inmax = 0.8 * np.nanmax(imdata)
    refmin = -10#np.nanmin(refdata)
    try:
        refmax = refheader['SATURATE']
    except:
        refmax = 0.8 * np.nanmax(refdata)
    #refmax = np.nanmax(refdata)

    try:
        imgain = imheader['GAIN']
    except:
        imgain = 1

    try:
        refgain = refheader['GAIN']
    except:
        refgain = 1
    if imgain == 0 or imgain > 10:
        imgain = 1.0
        
    if refgain == 0 or refgain > 10:
       refgain = 1.0

    
    # replace borders with NaNs in ref image if there are any that are == 0,
    hdulist=fits.open(filelist[0])
    ycenter1, xcenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
    hdulist[0].data=hdulist[0].data[ymin:ymax, xmin:xmax]

    ycenter, xcenter = int(ycenter1 - ymin)+1,  int(xcenter1- xmin)+1
    hdulist[0].header['CRPIX1'] = xcenter
    hdulist[0].header['CRPIX2'] = ycenter
    #print (np.min(hdulist[0].data))
    inimmed=np.median(hdulist[0].data)
    inimmax=np.max(hdulist[0].data)
    inimmin=np.min(hdulist[0].data)
    hdulist.writeto(filelist[0],overwrite=True)

    # perform astrometric calibration after cutting image on input image only
    scamp(filelist[0], config, useweight=False, CheckPlot=False, verbose='NORMAL')
    

    hdulist=fits.open(filelist[1])
    ycenter1, xcenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
    newimg = hdulist[0].data[ymin:ymax, xmin:xmax]
    #set max pixel value to 65000
    factor = np.max(newimg)
    hdulist[0].data=newimg#/200#/factor*10000
    #hdulist[0].data[hdulist[0].data == 0] = 1e-30

    ycenter, xcenter = int(ycenter1 - ymin)+1,  int(xcenter1- xmin)+1
    hdulist[0].header['CRPIX1'] = xcenter
    hdulist[0].header['CRPIX2'] = ycenter

    #hdulist[0].data=hdulist[0].data *1000
    #hdulist[0].header['GAIN']=5.0
    #print (np.min(hdulist[0].data))
    refmed=np.median(hdulist[0].data)
    refmax=np.max(hdulist[0].data)
    refmin=np.min(hdulist[0].data)
    hdulist.writeto(filelist[1],overwrite=True)

    if len(filelist) == 3:
    
        hdulist=fits.open(filelist[2])
        ycenter1, xcenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
        newdata=hdulist[0].data[ymin:ymax, xmin:xmax]
        newdata[newdata!=0]=1e8
        hdulist[0].data=newdata

        ycenter, xcenter = int(ycenter1 - ymin)+1,  int(xcenter1- xmin)+1
        hdulist[0].header['CRPIX1'] = xcenter
        hdulist[0].header['CRPIX2'] = ycenter
        hdulist.writeto(filelist[2],overwrite=True)
    
    return [[xmin, xmax, ymin, ymax], [inmin, inmax, refmin, refmax], [imgain, refgain]]

if __name__ == '__main__':

    im1 = 'Test_Atlascow/gmadet_results/ATLAS18qqn-S001-R001-C001-SDSS_g.fits'
    im2 = 'ps1_resample/ps1_mosaic.fits'
    registration(im1, im2, False)
