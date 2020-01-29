#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, errno, subprocess, glob, math, shutil, os
from astropy.io import ascii, fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from astropy.table import Table
from shapely.geometry import Polygon

from utils import rm_p, mkdir_p


def get_crpix(proj_crpix1,proj_crpix2, Xcell, Ycell, x, y):
    """Compute CRPIX1 and CRPIX2 for cell based on the CRPIX values of the projcell """
    x_center, y_center = 5, 5
    cprix1 = proj_crpix1 + (x_center - x) * (Xcell - 480)
    crpix2 = proj_crpix2 + (y_center - y) * (Ycell - 480)
    
    return cprix1, crpix2

def get_RADEC_coord(proj_crpix1, proj_crpix2, Xcell, Ycell, x, y, RA, Dec):

    pixscale=0.25/3600
    crpix1, crpix2 = get_crpix(proj_crpix1,proj_crpix2, Xcell, Ycell, x, y)

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)

    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [float(crpix1), float(crpix2)]
    w.wcs.cdelt = np.array([-pixscale, pixscale])
    w.wcs.crval = [RA, Dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    #w.wcs.set_pv([(2, 1, 45.0)])

    # Three pixel coordinates of interest.
    # The pixel coordinates are pairs of [X, Y].
    # The "origin" argument indicates whether the input coordinates
    # are 0-based (as in Numpy arrays) or
    # 1-based (as in the FITS convention, for example coordinates
    # coming from DS9).
    pixcrd = np.array([[0,0],[Xcell, 0], [Xcell, Ycell], [0, Ycell]], dtype=np.float64)

    # Convert pixel coordinates to world coordinates.
    # The second argument is "origin" -- in this case we're declaring we
    # have 0-based (Numpy-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    world = np.array(world).T

    return world, w

def ps1_cell_coord(im_corner_coords, projcell_id, Xcell, Ycell, projcell_ra_center, projcell_dec_center, proj_crpix1, proj_crpix2):
    """
    Template for computing the 10x10 cells composing a PS1 projcell
    """
    # cell size, 0yx
    ny = 10   # bottom to top of projell (ascending dec)
    nx = 10   # left to right of projcell (ascending ra)

    id_list = []
    RA_min = []
    RA_max = []
    dec_min = []
    dec_max = []
    projcell_id_list = []
    for y in range(ny):
        for x in range(nx):

            corner_coords, w = get_RADEC_coord(proj_crpix1, proj_crpix2, Xcell, Ycell, x, y, projcell_ra_center, projcell_dec_center)
            #print (projcell_id,y,x, corner_coords)
            # Check whether at least one cell corner is contained in the input image
            # using its coordinates
            pix_im_coord = np.array([im_corner_coords[0], im_corner_coords[1]]).T
            pix_cell_coord = np.array([corner_coords[0], corner_coords[1]]).T

            shapely_poly = Polygon([tuple(co) for co in pix_im_coord])
            shapely_poly2 = Polygon([tuple(co) for co in pix_cell_coord])

            #print (polygon_im.contains(cell_corner_coords, w))
            if shapely_poly.intersects(shapely_poly2):
                projcell_id_list.append(projcell_id)
                id_list.append('0%d%d' % (y,x))
                RA_min.append(np.min(corner_coords[0], axis=0))
                RA_max.append(np.max(corner_coords[0], axis=0))
                dec_min.append(np.min(corner_coords[1], axis=0))
                dec_max.append(np.max(corner_coords[1], axis=0))

    overlap_cells = Table([projcell_id_list, id_list, RA_min, RA_max, dec_min, dec_max],
                         names=('projcell_id', 'cell_id', 'RA_min', 'RA_max', 'dec_min', 'dec_max'))
    return overlap_cells


def ps1_grid(im_corner_coords):
    """
    Return the ps1 cell IDs for a given image dimension
    Skycell images have names like skycell.nnnn.0yx where nnnn 
    is the projection cell number (which ranges from 635 to 2643)
    and 0yx gives the skycell location in the image, with y and x
    ranging from 0 to 9 indicating the respective y and x section
    of the projection cell. The 000 skycell is in the bottom left
    corner of the projection cell, 010 is just above it, and 099 
    is in the upper right corner.
    RA = n *360 deg / M
    """
    # RA, dec min and max of input image
    ra_min = np.min(im_corner_coords[0], axis=0)
    ra_max = np.max(im_corner_coords[0], axis=0)
    dec_min = np.min(im_corner_coords[1], axis=0)
    dec_max = np.max(im_corner_coords[1], axis=0)
    
    ps1grid = Table.read('ps1_survey/ps1grid.fits', hdu=1)
 
    # Get the min and max declination zones
    mask = ((ps1grid['DEC_MIN'] < dec_min) & (ps1grid['DEC_MAX'] > dec_min)) | ((ps1grid['DEC_MIN'] < dec_max) & (ps1grid['DEC_MAX'] > dec_max))

    # Get all declinations zones
    all_zones_id = np.arange(ps1grid[mask]['ZONE'][0],ps1grid[mask]['ZONE'][-1]+1)
    
    #print (all_zones_id)
    projcell_idx_list = []
    all_cells = []

    # Loop over the different zones
    for zone in all_zones_id:
        mask = ps1grid['ZONE'] == zone
        idx_bkp = -1
        for ra in [ra_min, ra_max]:
            
            # Get the cells covering the input ra
            closet_projcell_idx = float(ps1grid[mask]['PROJCELL']) + ra * float(ps1grid[mask]['NBAND']) / 360
            projcell_idx = int(np.rint(closet_projcell_idx))

            if projcell_idx != idx_bkp:
                projcell_idx_list.append(projcell_idx)
            idx_bkp = projcell_idx
    
        total_proj_cell_idx = np.arange(projcell_idx_list[0],projcell_idx_list[-1]+1)
        for cell_id in total_proj_cell_idx:
            diff_projcell_idx = cell_id - float(ps1grid[mask]['PROJCELL'])
            ra_center_projcell = diff_projcell_idx * 360 / float(ps1grid[mask]['NBAND'])
            all_cells.append(ps1_cell_coord(im_corner_coords,
                                                cell_id,
                                                ps1grid[mask]['XCELL'],
                                                ps1grid[mask]['YCELL'],
                                                ra_center_projcell,
                                                ps1grid[mask]['DEC'],
                                                ps1grid[mask]['CRPIX1'],
                                                ps1grid[mask]['CRPIX2']
                                                ))
    
    if len(all_cells) == 1:
        all_cells = all_cells[0]
    else:
        all_cells = vstack([tab for tab in all_cells])
        
    return all_cells

def download_ps1_cells(cell_table, band):
    """Download the required cell from PS1 DR1"""
    
    ps1Dir = '/home/corre/codes/Test_hotpants/ps1Dir/'
    if not os.path.isdir(ps1Dir):
        os.makedirs(ps1Dir)

    ps1ResampleDir = '/home/corre/codes/Test_hotpants/ps1_resample/'
    if not os.path.isdir(ps1ResampleDir):
        os.makedirs(ps1ResampleDir)

    #if os.path.isfile(imagefile):
    #    return
    
    file_list = []
    
    BaseURL = "http://ps1images.stsci.edu/"
    for cell in cell_table:
        cell_url_path = '/rings.v3.skycell/%s/%s/' % (cell['projcell_id'],
                                                      cell['cell_id'])
        cell_file = 'rings.v3.skycell.%s.%s.stk.%s.unconv.fits' % (cell['projcell_id'],
                                                                   cell['cell_id'],
                                                                   band)
        Link = BaseURL + cell_url_path + cell_file
        local_cell_file = cell_file.replace(".","_").replace("_fits",".fits")
        FileNameFitsPath = ps1Dir + local_cell_file
        if os.path.isfile(FileNameFitsPath):
            print ('File %s already downloaded' % FileNameFitsPath)
        else:
            wget_command = "wget %s -O %s"%(Link,FileNameFitsPath)
            os.system(wget_command)

            # do not why, size reduction?
            funpack_command = "fpack %s; rm %s; funpack %s.fz"%(FileNameFitsPath,FileNameFitsPath,FileNameFitsPath)
            os.system(funpack_command)

            rm_command = "rm funpack %s.fz"%(FileNameFitsPath)
            os.system(rm_command)

            linear_rescale_ps1(local_cell_file, ps1Dir, ps1ResampleDir)

        file_list.append(ps1ResampleDir+local_cell_file)
            
        
    create_ps1_mosaic(file_list, ps1ResampleDir)

    return True

def linear_rescale_ps1(file, inputDir, outputDir):
    """resample the PS1 DR1 fits file"""
    
    # Transform into linear flux scale
    hdulist=fits.open(inputDir+file)
    
    boffset = hdulist[0].header['BOFFSET']
    bsoften = hdulist[0].header['BSOFTEN']
    a = 2.5/np.log(10)
    hdulist[0].data = boffset + 2 *bsoften* np.sinh(hdulist[0].data/a)
    
    hdulist.writeto(outputDir+file,overwrite=True)
    hdulist.close()

    return True

def create_ps1_mosaic(file_list, inputDIR):
    """Create a single mosaic of PS1 image using swarp"""
    useweight=False
    gain = 0
    
    np.savetxt('register.list', file_list, fmt='%s')
    
    imalists=['@' + 'register.list']
    if useweight:
        subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', epoch + '.fits', \
                             '-WEIGHTOUT_NAME', epoch + '.weight.fits', \
                             '-GAIN_DEFAULT', str(gain)] + [imalist])
    else:
        subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', inputDIR+'ps1_mosaic.fits',\
                             '-GAIN_DEFAULT', str(gain),\
                             '-SUBTRACT_BACK', 'N', \
                             '-COMBINE', 'Y', \
                             '-BACK_SIZE', '128', \
                             '-BACK_FILTERSIZE', '3',\
                             #'-CENTER', '%s, %s' % (ra,dec), \
                             '-RESAMPLING_TYPE', 'LANCZOS3',\
                             #'-RESAMPLING_TYPE', 'NEAREST', \
                             '-OVERSAMPLING', '0',\
                             '-COMBINE_TYPE', 'MEDIAN'] + imalists)

        # replace borders with NaNs in ref image if there are any that are == 0,
        #hdulist=fits.open(epoch + '.fits')
        #hdulist[0].data[hdulist[0].data==0]=np.nan
        #hdulist.writeto(epoch + '.fits',overwrite=True)

