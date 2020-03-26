#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, errno, subprocess, glob, math, shutil, os
from astropy.io import ascii, fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from astropy.table import Table, vstack
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

    # Set up projection
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

    # Need to append a 0 when < 1000 otherwise download crashes
    if projcell_id < 1000:
        projcell_id = '0' + str(projcell_id)
    else:
        projcell_id = str(projcell_id)

    id_list = []
    RA_min = []
    RA_max = []
    dec_min = []
    dec_max = []
    projcell_id_list = []
    for y in range(ny):
        for x in range(nx):
            # Estimate corner coordinates for each cell
            corner_coords, w = get_RADEC_coord(proj_crpix1, proj_crpix2, Xcell, Ycell, x, y, projcell_ra_center, projcell_dec_center)
            #print (projcell_id,y,x, corner_coords)
            # Check whether one cell is contained in the input image
            pix_im_coord = np.array([im_corner_coords[0], im_corner_coords[1]]).T
            pix_cell_coord = np.array([corner_coords[0], corner_coords[1]]).T

            im_poly = Polygon([tuple(co) for co in pix_im_coord])
            cell_poly = Polygon([tuple(co) for co in pix_cell_coord])

            #print (polygon_im.contains(cell_corner_coords, w))
            if im_poly.intersects(cell_poly):
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
            cell_query = ps1_cell_coord(im_corner_coords,
                                        cell_id,
                                        ps1grid[mask]['XCELL'],
                                        ps1grid[mask]['YCELL'],
                                        ra_center_projcell,
                                        ps1grid[mask]['DEC'],
                                        ps1grid[mask]['CRPIX1'],
                                        ps1grid[mask]['CRPIX2']
                                        )
            if len(cell_query) > 0:
                all_cells.append(cell_query)

    if len(all_cells) == 1:
        all_cells = all_cells[0]
    else:
        all_cells = vstack([tab for tab in all_cells])
        
    return all_cells

def download_ps1_cells(cell_table, band, inputimage):
    """Download the required cell from PS1 DR1"""
   
    path, filenameInput = os.path.split(inputimage)
    if path:
       folder = path + '/'
    else:
       folder = ''

    ps1Dir = 'ps1Dir/'
    if not os.path.isdir(ps1Dir):
        os.makedirs(ps1Dir)

    ps1ResampleDir = 'ps1_resample/'
    if not os.path.isdir(ps1ResampleDir):
        os.makedirs(ps1ResampleDir)

    #if os.path.isfile(imagefile):
    #    return
    
    file_list = []
    # extension to get auxiliary images
    # See https://outerspace.stsci.edu/display/PANSTARRS/PS1+Stack+images
    #auxiliaryFiles = ['.mask', '.wt', '.num', '.exp', '.expwt', '']
    auxiliaryFiles = ['.mask', '']

    BaseURL = "http://ps1images.stsci.edu/"
    for cell in cell_table:
        cell_url_path = 'rings.v3.skycell/%s/%s/' % (cell['projcell_id'],
                                                      cell['cell_id'])


        for aux in auxiliaryFiles:
            cell_file = 'rings.v3.skycell.%s.%s.stk.%s.unconv%s.fits' % (cell['projcell_id'],
                                                                         cell['cell_id'],
                                                                         band,
                                                                         aux)
            Link = BaseURL + cell_url_path + cell_file
            local_cell_file = cell_file.replace(".","_").replace("_fits",".fits")
            FileNameFitsPath = ps1Dir + local_cell_file
            if os.path.isfile(FileNameFitsPath):
                print ('File %s already downloaded' % FileNameFitsPath)
            else:
                wget_command = "wget %s -O %s"%(Link,FileNameFitsPath)
                os.system(wget_command)
                if os.path.isfile(FileNameFitsPath):
                    # do not really understand what this is doing
                    funpack_command = "fpack %s; rm %s; funpack %s.fz"%(FileNameFitsPath,FileNameFitsPath,FileNameFitsPath)
                    os.system(funpack_command)

                    rm_command = "rm funpack %s.fz"%(FileNameFitsPath)
                    os.system(rm_command)
                else:
                    print ("File %s was not downloaded or found on the server." % Link)

            if aux == '':
                # Check if targeted file was downloaded to continue
                if os.path.isfile(FileNameFitsPath):
                    if os.path.isfile(ps1ResampleDir + local_cell_file):
                        pass
                    else:
                        linear_rescale_ps1(local_cell_file, ps1Dir, ps1ResampleDir)
                
                    file_list.append(ps1ResampleDir+local_cell_file)

    # Create mosaic file if it does not exist
    mosaicfile = folder + filenameInput.split('.')[0] + '_ps1_mosaic.fits'
    if os.path.isfile(mosaicfile):
        print ('PS1 mosaic image already exists in this location: %s. If you want to recompute it, delete it.' % mosaicfile)
    else:
        create_ps1_mosaic(file_list, inputimage, folder)

    return True

def linear_rescale_ps1(filename, inputDir, outputDir, normalise=True):
    """resample the PS1 DR1 fits file"""
    
    # Transform into linear flux scale
    hdulist=fits.open(inputDir+filename)
    #hdulist_exp=fits.open(inputDir+filename.split('.')[0]+'_exp.fits')
    #hdulist_expwt=fits.open(inputDir+filename.split('.')[0]+'_expwt.fits')
    
    boffset = hdulist[0].header['BOFFSET']
    bsoften = hdulist[0].header['BSOFTEN']
    a = 2.5/np.log(10)
    hdulist[0].data = boffset + 2 *bsoften* np.sinh(hdulist[0].data/a)

    if normalise:
        
        hdr = hdulist[0].header
        SCLlist = hdr['SCL_*']
        scale_flux = []
        for SCL in SCLlist:
            if float(hdr[SCL]) > 0:
                scale_flux.append(1)
            else:
                scale_flux.append(0)
        explist = hdr['EXP_*']
        exposure = 0
        for i, exp in enumerate(explist):
            exposure += float(scale_flux[i]) * float(hdr[exp])
        #hdulist[0].data /= hdulist[0].header['EXPTIME']
        hdulist[0].data /= exposure

        try:
            #hdulist[0].header['SATURATE'] /= hdulist[0].header['EXPTIME']
            hdulist[0].header['SATURATE'] /= exposure
        except:
            pass 
        
        # Normalise by exacte exposure time in each pixels
        #hdulist[0].data = hdulist[0].data / hdulist_exp[0].data
        hdulist[0].header['EXPTIME'] = 1 
    hdulist.writeto(outputDir+filename,overwrite=True)
    hdulist.close()

    # replace pixels == 0 with NaNs. Mostly the border, saturated pixels
    hdulist=fits.open(outputDir+filename)
    hdulist[0].data[hdulist[0].data==0]=np.nan
    hdulist.writeto(outputDir+filename,overwrite=True)

    # Create a mask to propagate the nan pixels 
    hdulist=fits.open(outputDir+filename)
    hdulist[0].data[np.isfinite(hdulist[0].data)]=0
    hdulist[0].data[np.isnan(hdulist[0].data)]=1
    hdulist.writeto(outputDir+filename.split('.')[0] + '_mask.fits',overwrite=True)

    return True

def create_ps1_mosaic(file_list, inputimage, outputDir, useweight=False):
    """Create a single mosaic of PS1 image using swarp"""
    _, filenameInput = os.path.split(inputimage)

    # Create list of mask fits
    mask_list = [ima.split('.')[0] + '_mask.fits' for ima in file_list]

    np.savetxt('mosaic.list', file_list, fmt='%s')
    np.savetxt('mask.list', mask_list, fmt='%s')
    
    imagefiles = [outputDir + filenameInput.split('.')[0] + '_ps1_mosaic',
                  outputDir + filenameInput.split('.')[0] + '_ps1_mosaic_mask']

    # Get pixel scale from input image header
    header = fits.getheader(inputimage)

    try:
        pixScale = abs(header['CDELT1'])
    except Exception:
        try:
            pixScale = abs(header['_DELT1'])
        except Exception:
            try:
                pixScale = abs(header['CD1_1'])
            except Exception:
                print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1')
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
                        point + '.head' ] + [inputimage])
    # Some keywords manipulation using sed
    subprocess.call(['sed', '-i', \
                             's/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /', \
                     point + '.head'])



    imalists=[['@' + 'mosaic.list'], ['@' + 'mask.list']]
    for i, imagefile in enumerate(imagefiles):
        # Remove mosaic if already exists
        rm_p(imagefile + '.fits')
        if 'mask' in imagefile:
            subBackground = 'N'
        else:
            subBackground = 'Y'

        # Copy the common header in the .head file
        # So that it is read by sawrp for each image
        shutil.copy(point + '.head', imagefile + '.head')

        if useweight:
            subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', imagefile + '.fits', \
                             '-WEIGHTOUT_NAME', imagefile + '.weight.fits'] + imalists[i])
        else:
            subprocess.call(['swarp',
                             '-IMAGEOUT_NAME', imagefile + '.fits',\
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
                             '-COPY_KEYWORDS', ' PIXEL_SCALE'] + imalists[i])

    # replace pixels == 0 with NaNs. Mostly the border, saturated pixels
    hdulist=fits.open(imagefiles[0] + '.fits')
    hdulist[0].data[hdulist[0].data==0]=np.nan
    hdulist.writeto(imagefiles[0] + '.fits',overwrite=True)

    # Create a mask to propagate the nan pixels 
    hdulist=fits.open(imagefiles[1] + '.fits')
    hdulist[0].data[hdulist[0].data > 0] = 1
    hdulist[0].data[np.isnan(hdulist[0].data)]=1
    hdulist.writeto(imagefiles[1] + '.fits',overwrite=True)

    #for ima in file_list:
    #    rm_p(ima)
    #for ima in mask_list:
    #    rm_p(ima)
    rm_p('mosaic.list')
    rm_p('mask.list')
    rm_p('swarp.xml')
    #rm_p('coadd.weight.fits')

if __name__ == "__main__":

    im_coords = [[121.18, 121.18, 121.05, 121.05],
                 [ 38.3,  38.46,  38.46,  38.3]]


    inputimage = 'Test_Atlascow/gmadet_results/ATLAS18qqn-S001-R001-C001-SDSS_g.fits'
    
    header = fits.getheader(inputimage)
    Naxis1 = header['NAXIS1']
    Naxis2 = header['NAXIS2']
    print (fits.getdata(inputimage).data.shape)
    print (Naxis1, Naxis2)

    pix_coords = [[0,0,Naxis1,Naxis1], [0,Naxis2,Naxis2,0]]

    from astropy import wcs
    from astropy.wcs import WCS

    from astropy.coordinates import SkyCoord
    from astropy import units as u
 
    # Get physical coordinates of OT
    w = WCS(header)
    ra, dec = w.all_pix2world(pix_coords[0],pix_coords[1], 1) 
    print (ra, dec)
    #stop
    im_coords = [ra, dec]
    cell_table = ps1_grid(im_coords)
    print (cell_table)
    download_ps1_cells(cell_table, 'g', inputimage)
