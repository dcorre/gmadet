#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, subprocess, shutil
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from utils import rm_p, mkdir_p 
from copy import deepcopy
from astrometry import scamp

def registration(filelist, config, resultDir='', useweight=False, gain=1, normalise_exp=True, verbose='NORMAL'):
    """Register images"""

    # Initialise lists used for creating output astropy table
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
        np.savetxt('register.list', files, fmt='%s')

        # Get pixel scale from input image header
        header = fits.getheader(inim)
        pixScale = abs(float(header['CDELT1'])) * 3600

        imalists=['@' + 'register.list']
        # File name to store the common header that will be shared by all
        # images in filelist
        point = 'registration'
        # Delete if already exists
        rm_p(point + '.head')
        # First run swarp to create a .head file containing the shared header
        subprocess.call(['swarp', '-HEADER_ONLY', 'Y', '-IMAGEOUT_NAME', \
                         point + '.head', '-GAIN_DEFAULT', str(gain), \
                         #'-VERBOSE_TYPE', verbose] + imalists)
                         '-VERBOSE_TYPE', verbose] + [inim])
        # Some keywords manipulation using sed
        subprocess.call(['sed', '-i', \
                         's/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /', \
                         point + '.head'])
        outFiles=[]
        # Run swarp to perform the registration on each image in filelist
        for j, ima in enumerate(files):
            if 'mask' in ima:
                subBackground = 'N'
            else:
                subBackground = 'Y'
            path, filename_ext = os.path.split(ima)
            epoch = resultDir + filename_ext.split('.')[0] + '_reg_%s' % i
            outFiles.append(epoch+'.fits')
            # Copy the common header in the .head file 
            # So that it is read by sawrp for each image
            shutil.copy(point + '.head',epoch + '.head')
            if useweight:
                subprocess.call(['swarp',
                                 '-IMAGEOUT_NAME', epoch + '.fits', \
                                 '-WEIGHTOUT_NAME', epoch + '.weight.fits', \
                                 '-VERBOSE_TYPE', verbose, \
                                 '-GAIN_DEFAULT', str(gain)] + [ima])
            else:
                # Use bilinear to avoid artefact, but worst for noise so need to check
                # in more details
                subprocess.call(['swarp',
                                 '-IMAGEOUT_NAME', epoch + '.fits',\
                                 #'-GAIN_DEFAULT', str(gain),\
                                 '-FSCALE_KEYWORD', 'NONE', \
                                 '-FSCALE_DEFAULT', '1, 1', \
                                 '-SUBTRACT_BACK', subBackground, \
                                 '-COMBINE', 'Y', \
                                 '-BACK_SIZE', '128', \
                                 '-BACK_FILTERSIZE', '3',\
                                 '-RESAMPLE', 'Y',\
                                 '-PIXELSCALE_TYPE', 'MANUAL', \
                                 '-PIXEL_SCALE', str(pixScale), \
                                 #'-CENTER', '%s, %s' % (header['CRVAL1'],header['CRVAL2']), \
                                 #'-RESAMPLING_TYPE', 'LANCZOS3',\
                                 '-RESAMPLING_TYPE', 'BILINEAR',\
                                 #'-RESAMPLING_TYPE', 'NEAREST', \
                                 '-OVERSAMPLING', '0',\
                                 '-COMBINE_TYPE', 'MEDIAN', \
                                 '-VERBOSE_TYPE', verbose, \
                                '-COPY_KEYWORDS', 'FILTER'] + [ima])

            # replace borders with NaNs in ref image if there are any that are == 0,
            #hdulist=fits.open(epoch + '.fits')
            #hdulist[0].data[hdulist[0].data==0]=np.nan
            #hdulist.writeto(epoch + '.fits',overwrite=True)

            rm_p(epoch + '.head')
        rm_p(point+'.head')
        rm_p('register.list')
        rm_p('coadd.weight.fits')

        inim_regist = outFiles[0] 
        refim_regist = outFiles[1]
        maskim_regist = outFiles[2]

        # Set specific value for mask image
        hdulist = fits.open(inim_regist)
        if normalise_exp:
            hdulist[0].data = hdulist[0].data / hdulist[0].header['EXPTIME']
            try:
                hdulist[0].header['SATURATE'] = hdulist[0].header['SATURATE'] / hdulist[0].header['EXPTIME']
            except:
                pass
            hdulist[0].header['EXPTIME'] = 1
        # 1e-30 is the default value of bad pixels for hotpants
        hdulist[0].data[hdulist[0].data == 0] = 1e-30
        #hdulist[0].data[hdulist[0].data != 1e-30] /= 20
        hdulist.writeto(inim_regist, overwrite=True)

        hdulist = fits.open(refim_regist)
        if normalise_exp:
            hdulist[0].data = hdulist[0].data / hdulist[0].header['EXPTIME']
            try:
                hdulist[0].header['SATURATE'] = hdulist[0].header['SATURATE'] / hdulist[0].header['EXPTIME']
            except:
                pass
            hdulist[0].header['EXPTIME'] = 1

        # 1e-30 is the default value of bad pixels for hotpants
        hdulist[0].data[hdulist[0].data == 0] = 1e-30
        #hdulist[0].data[hdulist[0].data != 1e-30] /= 20
        hdulist.writeto(refim_regist, overwrite=True)

        hdulist = fits.open(maskim_regist)
        hdulist[0].data[hdulist[0].data != 0] = 1e8
        hdulist.writeto(maskim_regist, overwrite=True)

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

    info = Table([inim_list, refim_list, mask_list, XY_lim, in_lo, in_up, ref_lo, ref_up, gain_in, gain_ref], names=['inim', 'refim', 'mask', 'XY_lim', 'in_lo', 'in_up', 'ref_lo', 'ref_up', 'gain_in', 'gain_ref'])
    return info

def get_hotpants_info(filelist, config, verbose):
    """Get some information for tuning hotpants
    """

    # Get the min and max pix with non 0 values 
    # to delimate the new frame
    imdata, imheader = fits.getdata(filelist[0], header=True)
    refdata, refheader = fits.getdata(filelist[1], header=True)
    ymin1, xmin1 = np.min(np.where(imdata>1e-30),axis=1)
    ymax1, xmax1 = np.max(np.where(imdata>1e-30),axis=1)
    ymin2, xmin2 = np.min(np.where(refdata>1e-30),axis=1)
    ymax2, xmax2 = np.max(np.where(refdata>1e-30),axis=1)

    # Take only non 1e-30 values
    # Even if it removes some part of the image on the edges.
    # Hotpants works better like this
    xmin=np.max([xmin1,xmin2])
    xmax=np.min([xmax1,xmax2])
    ymin=np.max([ymin1,ymin2])
    ymax=np.min([ymax1,ymax2])
    
    inmin = -10 #np.nanmin(imdata)
    try:
        inmax = 0.9 * imheader['SATURATE']
    except:
        inmax = 0.9 * np.nanmax(imdata)
    refmin = -10 #np.nanmin(refdata)
    try:
        refmax = 0.9 * refheader['SATURATE']
    except:
        refmax = 0.9 * np.nanmax(refdata)
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
    xcenter1, ycenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
    hdulist[0].data=hdulist[0].data[ymin:ymax, xmin:xmax]

    xcenter, ycenter = int(xcenter1 - xmin)+1,  int(ycenter1- ymin)+1
    hdulist[0].header['CRPIX1'] = xcenter
    hdulist[0].header['CRPIX2'] = ycenter
    #print (np.min(hdulist[0].data))
    inimmed=np.median(hdulist[0].data)
    inimmax=np.max(hdulist[0].data)
    inimmin=np.min(hdulist[0].data)
    hdulist.writeto(filelist[0],overwrite=True)
    # perform astrometric calibration after cutting image on input image only
    # This astrometric solution will be shared with the reference file and 
    # resulting substracted image where sources are analysed.
    #scamp(filelist[0], config, useweight=False, CheckPlot=False, verbose=verbose)

    hdulist=fits.open(filelist[1])
    xcenter1, ycenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
    newimg = hdulist[0].data[ymin:ymax, xmin:xmax]
    #set max pixel value to 65000
    factor = np.max(newimg)
    hdulist[0].data=newimg#/200#/factor*10000
    #hdulist[0].data[hdulist[0].data == 0] = 1e-30

    xcenter, ycenter = int(xcenter1 - xmin)+1,  int(ycenter1- ymin)+1
    hdulist[0].header['CRPIX1'] = xcenter
    hdulist[0].header['CRPIX2'] = ycenter

    #hdulist[0].data=hdulist[0].data *1000
    #hdulist[0].header['GAIN']=5.0
    #print (np.min(hdulist[0].data))
    #refmed=np.median(hdulist[0].data)
    #refmax=np.max(hdulist[0].data)
    #refmin=np.min(hdulist[0].data)
    hdulist.writeto(filelist[1],overwrite=True)

    if len(filelist) == 3:
    
        hdulist=fits.open(filelist[2])
        xcenter1, ycenter1 = hdulist[0].header['CRPIX1'],  hdulist[0].header['CRPIX2']
        newdata=hdulist[0].data[ymin:ymax, xmin:xmax]
        newdata[newdata!=0]=1e8
        hdulist[0].data=newdata

        xcenter, ycenter = int(xcenter1 - xmin)+1,  int(ycenter1- ymin)+1
        hdulist[0].header['CRPIX1'] = xcenter
        hdulist[0].header['CRPIX2'] = ycenter
        hdulist.writeto(filelist[2],overwrite=True)
    
    return [[xmin, xmax, ymin, ymax], [inmin, inmax, refmin, refmax], [imgain, refgain]]

if __name__ == '__main__':

    im1 = 'Test_Atlascow/gmadet_results/ATLAS18qqn-S001-R001-C001-SDSS_g.fits'
    im2 = 'ps1_rescaled/ps1_mosaic.fits'
    registration(im1, im2, False)
