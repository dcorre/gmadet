#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, subprocess
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from registration import registration
from ps1_survey import ps1_grid, download_ps1_cells
from utils import get_phot_cat
import hjson

def get_corner_coords(filename):
    """ Compute the RA, Dec coordinates at the corner of one image"""

    header = fits.getheader(filename)
    Naxis1 = header['NAXIS1']
    Naxis2 = header['NAXIS2']

    pix_coords = [[0,0,Naxis1,Naxis1], [0,Naxis2,Naxis2,0]]

    # Get physical coordinates of OT
    w = WCS(header)
    ra, dec = w.all_pix2world(pix_coords[0],pix_coords[1], 1)

    return [ra, dec]

def substraction(filenames, reference, config, method='hotpants'):
    """Substract a reference image to the input image"""

    imagelist = np.atleast_1d(filenames)
    subFiles = []
    for ima in imagelist:
        
        path, filename = os.path.split(ima)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get coordinates of input image 
        im_coords = get_corner_coords(ima)

        # Define the reference image
        if reference == 'ps1':
            _, band, _ = get_phot_cat(ima, None)
            if band == 'B':
                band = 'g'
            elif band == 'V':
                band = 'g'
            elif band == 'R':
                band = 'r'
            elif band == 'I':
                band = 'i'
            elif band == 'g+r':
                band = 'r'
            #band = 'g'
            cell_table = ps1_grid(im_coords)
            download_ps1_cells(cell_table, band, ima)
            refim = folder +  filename.split('.')[0] + '_ps1_mosaic.fits' 
            refim_mask = folder +  filename.split('.')[0] + '_ps1_mosaic_mask.fits'

        sub_info = registration(ima, refim, config, refim_mask, False)
        
        if method == 'hotpants':
            ima_regist = folder + 'substraction/' + filename.split('.')[0] + '_regist.fits'
            refim_regist = folder + 'substraction/' + filename.split('.')[0] + '_ps1_mosaic_regist.fits'
            refim_regist_mask = folder + 'substraction/' + filename.split('.')[0] + '_ps1_mosaic_mask_regist.fits'
            subfiles = hotpants(ima_regist, refim_regist, config, sub_info, refim_mask=refim_regist_mask)

        subFiles.append(subfiles)

    return subFiles

def hotpants(inim, refim, config, sub_info, refim_mask=None):
    """Image substraction using hotpants"""

    path, _ = os.path.split(inim)
    resfile = inim.split('.')[0] + '_sub.fits'
    maskfile = inim.split('.')[0] + '_sub_mask.fits'

    with open(config['hotpants']['conf']) as json_file: 
        hotpants_conf = hjson.load(json_file)

    # Set min and max acceptable values for input and template images
    # Too simple, need to adapt it in the future
    il = str(sub_info[1][0])
    iu = str(sub_info[1][1])
    tl = str(sub_info[1][2])
    tu = str(sub_info[1][3])
   
    overlap = '%s, %s, %s, %s' % (sub_info[0][0],
                                  sub_info[0][1],
                                  sub_info[0][2],
                                  sub_info[0][3],)

    hotpants_cmd = 'hotpants -inim %s -tmplim %s -outim %s -omi %s ' % (inim, refim, resfile, maskfile)
    #hotpants_cmd += '-il %s -iu %s -tl %s -tu %s -gd %s ' % (il, iu, tl, tu, overlap)
    hotpants_cmd += '-il %s -iu %s -tl %s -tu %s ' % (il, iu, tl, tu)
    hotpants_cmd += '-tuk %s -iuk %s ' % (tu, iu)
    hotpants_cmd += '-ig %s -tg %s ' % (sub_info[2][0], sub_info[2][1])

    if refim_mask:
        hotpants_cmd += '-tmi %s ' % refim_mask

    # Add params from the hjson conf file
    for key, value in hotpants_conf.items():
        hotpants_cmd += '-%s %s ' % (key, value)

    print (hotpants_cmd)

    
    os.system(hotpants_cmd)
    #subprocess.call([hotpants_cmd])
    
    # Set bad pixel values to 0 and others to 1 for sextractor
    hdulist = fits.open(maskfile)
    hdulist[0].data[hdulist[0].data == 0] = 1
    hdulist[0].data[hdulist[0].data != 1] = 0
    hdulist.writeto(maskfile, overwrite=True)

    return [inim, refim, resfile, maskfile]
