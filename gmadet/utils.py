#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import errno, glob, os, shutil, subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import vstack, Table

from astropy import wcs
from astropy.wcs import WCS
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.visualization import (MinMaxInterval, SqrtStretch, LogStretch,SinhStretch,
                                       ImageNormalize, ZScaleInterval)
import voeventparse as vp
import json
import requests
from copy import deepcopy
from shapely.geometry import Point, Polygon

def cp_p(src, dest):
  try:
    shutil.copy(src, dest)
  except:
    pass

def mv_p(src, dest):
  try:
    shutil.move(src, dest)
  except:
    pass

def rm_p(src):
  try:
    os.remove(src)
  except:
    pass

def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise

def load_config(telescope):
    """Load the path to the configuration files required by the softs.
       They are telescope dependent.
    """

    path2tel = 'config/' + telescope + '/'
    config = {
            'sextractor': {
                'conf': path2tel + 'sourcesdet.sex',
                'param': path2tel + 'sourcesdet.param',
                },
            'scamp': {
                'sextractor': path2tel + 'prepscamp.sex',
                'param': path2tel + 'prepscamp.param',
                'conf': path2tel + 'scamp.conf',
                },
            'swarp': {

                },
            'psfex': {
                'sextractor': path2tel + 'preppsfex.sex',
                'param': path2tel + 'preppsfex.param',
                'conf': path2tel + 'psfex.conf',
                },
            'hotpants': {
                'conf': path2tel + 'hotpants.hjson'
                }
            }

    return config


def clean_folder(filelist):
    """ Remove output files from previous iraf run. No need for sextractor  """

    types = ('*coo.*', '*mag.*', '*.magwcs', '*.magfiltered*')
    files2delete = []
    for filename in  filelist:
        path = os.path.split(filename)
        if path[0]:
            folder = path[0] + '/'
        else:
            folder = ''

        for f in types:
            files2delete.extend(glob.glob(folder+f))
            files2delete.extend(glob.glob(f))

    files2delete = np.unique(files2delete)
    for f in files2delete:
        os.remove(f)


def cut_image(filename, resultDir, Nb_cuts=(2,2)):

    path, filename_ext = os.path.split(filename)
    if path:
        folder = path + '/'
    else:
        folder = ''

    filename2 = filename_ext.split('.')[0]
    extension = ''
    for ext in filename_ext.split('.')[1:]:
        extension = extension + '.' + ext

    if Nb_cuts == (1,1):
        filename_out = resultDir + filename2 + extension
        cp_p(filename,filename_out)
    else:
        filename_out = resultDir + filename2 + extension
        cp_p(filename,filename_out)

        data, header = fits.getdata(filename, header=True)
        # if keyword with scamp PV keywords
        # perform a quick and dirty fix
        try:
            pvlist = header['PV*']
            for pv in pvlist:
                tpv = 'T'+pv
                header.rename_keyword(pv, tpv, force=False)
        except:
            pass
        w = wcs.WCS(header)
        Naxis1 = header['NAXIS1']
        Naxis2 = header['NAXIS2']
        Naxis11 = int(Naxis1/Nb_cuts[0])
        Naxis22 = int(Naxis2/Nb_cuts[1])

        index=0
        for i in range(Nb_cuts[0]):
            for j in range(Nb_cuts[1]):
                index+=1
                if i == 0 :
                    x1 = 1
                else:
                    x1 = Naxis11 * i + 1
                if j == 0 :
                    y1 = 1
                else:
                    y1 = Naxis22 * j + 1
                x2 = Naxis11 * (i+1)
                y2 = Naxis22 * (j+1)

                filename_out = resultDir + filename2 + "_Q%d" % (index) + extension
 
                # No need to update the header if astrometric calibration is performed with scamp
                # this will be updated later
                #datacut = data[x1-1:x2-1,y1-1:y2-1]
                datacut = data[y1-1:y2-1,x1-1:x2-1]
                newheader = deepcopy(header)
                crpix2 = int((y2-y1)/2)
                crpix1 = int((x2-x1)/2)
                ra, dec = w.wcs_pix2world(crpix1+Naxis11*i, crpix2+Naxis22*j, 1)
                newheader['CRPIX1'] = crpix1
                newheader['CRPIX2'] = crpix2
                newheader['CRVAL1'] = float(ra)
                newheader['CRVAL2'] = float(dec)
                fits.writeto(filename_out, datacut, newheader,overwrite=True)


def make_sub_image(filename, OT_coords, coords_type='world',
                   output_name='subimage.fits.gz', size=[200,200],
                   FoV=-1, fmt='png'):
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
    output_name : string, optional
        path, including the name, where to write the new image to be created.
        Without extension as the extension is automatically set to .fits.gz.
        Default: 'subimage'
    size : list, optional
        define the size in pixels of the new sub-image.
        Default: [200,200]
    FoV: float, optional
        define the FoV in arcsec for the subimage. If -1 then size is defined by size
    fmt: string, optional
        define the format of the subimage, 'png' or 'fits'
    Returns
    -------
    No variable is returned.
    A '.fits.gz' file is created using the path defined through 'output_name'.

    """
    # Load file
    data, header = fits.getdata(filename, header=True)

    # Get physical coordinates of OT
    w = WCS(header)
    if coords_type == 'world':
        # Get physical coordinates
        c = coord.SkyCoord(OT_coords[0], OT_coords[1], unit=(u.deg, u.deg),frame='icrs')
        world = np.array([[c.ra.deg, c.dec.deg]])
        #print (world)
        pix = w.all_world2pix(world,1)[0]
        pix_ref = OT_coords
    elif coords_type == 'pix':
        pix = OT_coords
        #print (pix)
        #ra, dec = w.all_pix2world(np.array(pix), 0)
        ra, dec = w.all_pix2world(pix[0],pix[1], 0)
        pix_ref = [float(ra),float(dec)]

    if FoV > 0:    
        # Get pixel size in degrees
        try:
            pixSize = abs(float(header['CDELT1']))
        except:
            pixSize = abs(float(header['CD1_1']))
        # Compute number of pixels to reach desired FoV in arcseconds
        size = [int(FoV/(pixSize*3600)), int(FoV/(pixSize*3600))]
    
    # Extract subimage from image starting from reference pixel
    subimage = data[int(pix[0]) - int(size[0]/2) : int(pix[0]) + int(size[0]/2),
                    int(pix[1]) - int(size[1]/2) : int(pix[1]) + int(size[1]/2)]

    if fmt == 'fits':
        # write new sub-image
        hdu = fits.PrimaryHDU()
        hdu.data = subimage.astype(np.uint16)
        # Need to adapt header here !!!
        header['CRPIX1'] = int(size[0]/2)
        header['CRPIX2'] = int(size[1]/2)
        header['CRVAL1'] = pix_ref[0]
        header['CRVAL2'] = pix_ref[1]

        hdu.header = header
        hdu.writeto(output_name,overwrite=True)

    elif fmt == 'png':
       norm = ImageNormalize(subimage, interval=ZScaleInterval(),
                      stretch=SinhStretch())
       plt.imshow(subimage, cmap='gray',origin='upper',norm=norm)
       plt.savefig(output_name)

def get_corner_coords(filename):
    """Get the image coordinates of an image"""

    header = fits.getheader(filename)
    # Get physical coordinates
    Naxis1 = header['NAXIS1']
    Naxis2 = header['NAXIS2']

    pix_coords = [[0,0,Naxis1,Naxis1], [0,Naxis2,Naxis2,0]]

    # Get physical coordinates
    w = WCS(header)
    ra, dec = w.all_pix2world(pix_coords[0],pix_coords[1], 1)

    return ra, dec

def send_data2DB(filename, candidates, Nb_cuts, owncloud_path, VOE_path, usrpwd_path, FoV=60, coords_type='pix', corner_cut=32,debug=False,fmt='png'):
    """Send candidates information to database"""

    # Load data and header
    data, header = fits.getdata(filename, header = True)

    dateObs = header['DATE-OBS']
    Tstart = Time(dateObs, format='fits', scale='utc')
    try:
        exposure = float(header['EXPOSURE'])
        Tend = Tstart + TimeDelta(exposure, format='sec')
    except:
        exposure = -1
        Tend = Time(dateObs, format='fits', scale='utc')
    
    # Try to get airmass rom header, else set it to -1
    try:
        Airmass = header['AIRMASS']
    except:
        Airmass = -1
    
    # Compute the image corner RA, Dec coordinates
    ra, dec = get_corner_coords(filename)
    pix_im_coord = np.array([ra, dec]).T
    im_poly = Polygon([tuple(co) for co in pix_im_coord])
    
    # Do not consider candidates found in the image edge
    imsize = data.shape
    #print (imsize, header['NAXIS1'], header['NAXIS2'])
    # Get the physical pixels of the original size if image were split into different quadrants.
    for i, candidate in enumerate(candidates):
        quadrant_idx = candidate['quadrant']
        if quadrant_idx == 'None':
            quadrant = None
            index_i = 0
            index_j = 0
        else:
            quadrant, index_i, index_j = quadrant_idx.split('_')
            quadrant = quadrant[1:]
        
        candidates['Xpos'][i] = candidate['Xpos'] + int(imsize[0]/Nb_cuts[0]) * int(index_j)
        candidates['Ypos'][i] = candidate['Ypos'] + int(imsize[1]/Nb_cuts[1]) * int(index_i)
    #print (candidates)
    mask = (candidates['Xpos'] > corner_cut) & \
            (candidates['Ypos'] > corner_cut) & \
            (candidates['Xpos'] < imsize[1] - corner_cut) & \
            (candidates['Ypos'] < imsize[0] - corner_cut)
    candidates_cut = candidates[mask]
    #print (candidates_cut)
    # Get information about the current alert from the xml file containing observation plan
    with open(VOE_path,'rb') as f:
        obsplan = vp.load(f)

    dict_event={}
    dict_event['event_type'] = obsplan.find(".//Param[@name='Event_type']").attrib['value']
    dict_event['event_name'] = obsplan.find(".//Param[@name='Event_ID']").attrib['value']
    dict_event['event_status'] = obsplan.find(".//Param[@name='Event_status']").attrib['value']
    dict_event['revision'] = obsplan.find(".//Param[@name='Revision']").attrib['value']
    dict_event['telescope'] = obsplan.find(".//Param[@name='Name_tel']").attrib['value']
    tiles_info = get_obsplan(obsplan)

    # Get user email adress and password to login in https://grandma-fa-interface.lal.in2p3.fr
    with open(usrpwd_path) as f:
        usrpwd = json.load(f)

    # Set up the output repository path to store sub-images
    outputDir = owncloud_path + '/' + dict_event['event_type'] + '/' + \
            dict_event['event_name'] + '/' + dict_event['event_status'] + \
            '_' + dict_event['revision'] + '/OTs/'

    # Create a sub image centered on each candidate found, and gather information
    #tile_id_list = [1] * len(candidates_cut)
    #filter_list = candidates_cut['filter_DB']
    #Tstart_list = [Tstart.fits] * len(candidates_cut)
    #Tend_list = [Tend.fits] * len(candidates_cut)
    #Airmass_list = [Airmass] * len(candidates_cut)
    Fits_path = []
    tile_id_list = []
    for i, row in enumerate(candidates_cut):
        name = dict_event['telescope'] + '_' + \
                str(round(float(row['_RAJ2000']),5)) + '_' + \
                str(round(float(row['_DEJ2000']),5)) + '_' + \
                dateObs + '.' + fmt
        #name = 'test' + str(i) + '.' + fmt
        Fits_path.append(name)
        if coords_type == 'world':
            OT_coords = [row['_RAJ2000'], row['_DEJ2000']]
        elif coords_type == 'pix':
            OT_coords = [row['Ypos'], row['Xpos']]
        make_sub_image(filename, OT_coords, coords_type=coords_type,
                       output_name=outputDir+name, size=[128,128], FoV=FoV, fmt=fmt)

        # Check in which tile the OT is located
        # Check that the RA, Dec lies within the FoV rectangle
        # Stop when finding one. Assuming that tiles are not overlapping.
        tile_id = 0 # by default
        for tile in tiles_info:
            tile_center = Point(tiles_info['RA'], tiles_info['Dec'])
            if tile_center.intersects(im_poly):
                tile_id = tiles_info['Id']
                break
        tile_id_list.append(tile_id)

    alias = ['new'] * len(candidates_cut)
    new = [1] * len(candidates_cut)
    RA_list = candidates_cut['_RAJ2000']
    Dec_list = candidates_cut['_DEJ2000']
    filter_list = candidates_cut['filter_DB']
    Tstart_list = [Tstart.fits] * len(candidates_cut)
    Tend_list = [Tend.fits] * len(candidates_cut)
    Mag_list = candidates_cut['mag_calib']
    Mag_err_list = candidates_cut['mag_calib_err']
    Magsys_list = candidates_cut['magsys']
    Airmass_list = [Airmass] * len(candidates_cut)

    candidates_2DB = Table([alias,new,tile_id_list,RA_list,Dec_list,filter_list,Tstart_list,Tend_list,Mag_list,Mag_err_list,Magsys_list,Airmass_list,Fits_path], names=['alias','new','tile_id','RA','DEC','filter','Tstart','Tend','Magnitude','Magnitude_error','Magsys','Airmass','fits_name'])

    
    # Set url to report tile or galaxy observations
    url = "https://grandma-fa-interface.lal.in2p3.fr/obs_report_OT.php"
    #url = "http://localhost/grandma/obs_report_OT.php"

    # Loop over the observations
    for i in range(len(candidates_2DB)):
        data2DB={}
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
            print ('Data sent succesfully to database.')
            forced_debug = False
        else:
            print ('Data not sent to database. See information below.')
            forced_debug = True

        if debug or forced_debug:
            print ('\nDEBUG:\n')
            print ('Data sent to DB:')
            print (data2DB)
            print ('\n\n')
            print ('Request response text:')
            print (response.text)
            print ('\n\n')
            print ('Request response status code:')
            print (response.status_code)
            print ('\n\n')
            print ('Request response history:')
            print (response.history)

            forced_debug = False

def get_obsplan(v):
    """Extract the tiles id and RA, Dec from voevent"""
    ID=[]
    Ra=[]
    Dec=[]
    Grade=[]
    Header=[]
    for element in v.What.iterchildren():
        tag_type=str(element.tag)
        if tag_type=="Table":
            for subelement in element.iterchildren():
                tag_type=str(subelement.tag)
                if tag_type=="Field":
                    Header.append(str(subelement.attrib['name']))
                    if tag_type=="Data":
                        for lines in subelement.iterchildren():
                            ID.append(int(lines.TD[0]))
                            Ra.append(float(lines.TD[1]))
                            Dec.append(float(lines.TD[2]))
                            Grade.append(float(lines.TD[3]))

    tiles_info = Table([ID,Ra,Dec,Grade], names=['Id', 'RA', 'Dec', 'Grade'])

    return tiles_info

def get_phot_cat(filename, telescope):
    """ Get the name of the filter band from header and telescope name
        And associate the correct name from DB
    """
    header = fits.getheader(filename)

    # FILTER keyword required
    try: 
        band = header['FILTER']
    except Exception:
        print('No FILTER keyword found in header.')

    if band in  ['C', 'Clear', 'NoFilter']:
        band_DB = 'Clear'
        band_cat = 'g+r'
    elif band in ['g', 'gSDSS']:
        band_DB = 'g/AB'
        band_cat = 'g'
    elif band in ['r', 'rSDSS', 'rPATH']:
        band_DB = 'r/AB'
        band_cat = 'r'
    elif band in ['i', 'iSDSS']:
        band_DB = 'i/AB'
        band_cat = 'i'
    elif band in ['z', 'zSDSS']:
        band_DB = 'z/AB'
        band_cat = 'z'
    elif band in ['B']:
        band_DB = 'B/Johnson'
        band_cat = 'B'
    elif band in ['V']:
        band_DB = 'V/Johnson'
        band_cat = 'V'
    elif band in ['R']:
        band_DB = 'R/Johnson'
        band_cat = 'R'
    elif band in ['I']:
        band_DB = 'I/Johnson'
        band_cat = 'I'

    # Chose photometric catalog
    try:
        RA = header['CRVAL1']
    except Exception:
        print('No CRVAL1 keyword found header. Astrometric calibration required.')
    try:
        DEC = header['CRVAL2']
    except Exception:
        print('No CRVAL2 keyword found header. Astrometric calibration required.')

    # Use Pan-Starrs if Dec > -30 degrees
    if float(DEC) > -30.:
        catalog = 'II/349/ps1' 
    # Else use SDSS if available.
    elif Vizier.query_region(field, width=rad_deg, height=rad_deg,
                             catalog="V/147/sdss12")[0]:
        catalog = "V/147/sdss12" 
    # Else use Gaia, but no calibration available for z band.
    elif filter not in z_band:
        catalog = "I/345/gaia2"
    # Last choice: USNO-B1. All-sky but bad photometric calibration.
    else:
        catalog = "I/284/out"

    return band_DB, band_cat, catalog


def unpackbits(x,num_bits):
    """ Unpack bits with any dimension ndarray.
        Can unpack however many bits"""
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    return (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])

def filter_catalog_data(data, catalogName):
    """Remove extended sources and bad measurements from reference catalogs
       before performing photometric calibration"""

    # Keep only point source objects and good measurements

    # Panstarrs flags
    if catalogName == 'II/349/ps1':
        # First remove data using the general 'Qual' flag.
        # ------------------------------------------------
        # There are 8 bits
        # Bit 1: extended object in PS1
        # Bit 2: Extended in external data (e.g. 2MASS)
        # Bit 3: Good-quality measurement in PS1
        # Bit 4: Good-quality measurement in external data (eg, 2MASS)
        # Bit 5: Good-quality object in the stack (>1 good stack measurement)
        # Bit 6: The primary stack measurements are the best measurements
        # Bit 7: Suspect object in the stack (no more than 1 good measurement, 2 or more suspect or good stack measurement)
        # Bit 8:  Poor-quality stack object (no more than 1 good or suspect measurement)

        Quality_flags = np.array(data['Qual'])
        Qual_flag = unpackbits(Quality_flags, 8)

        qual_bits = [1, 2, 3, 7, 8]
        qual_values = [0, 0, 1, 0, 0]
        counter=0
        for i, j in zip(qual_bits, qual_values):
            condition = Qual_flag[:,i-1] == j
            if counter == 0:
                quality_mask = condition
            else:
                quality_mask = np.bitwise_and(quality_mask, condition)
            counter+=1
        # flag for individual bands
        # -------------------------
        # There 25 bits. Use only the bit stating whether it is an extended
        # object in this band
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
        bands = ['g', 'r', 'i', 'z']
        band_flags = []
        # Unpack bits from individual band flags
        for band in bands:
            _temp = np.array(data['%sFlags' % band])
            band_flags.append(unpackbits(_temp, 25))
        band_flags = np.array(band_flags)

        # Apply mask conditions
        for i in range(len(bands)):
            counter=0
            for j1, k1 in zip(band_bits_and, band_values_and):
                condition = band_flags[i][:,j1-1] == k1
                quality_mask = np.bitwise_and(quality_mask, condition)
                counter+=1
            counter2=0
            # At least one Average magnitude calculated
            for j2, k2 in zip(band_bits_or, band_values_or):
                condition_or = band_flags[i][:,j2-1] == k2
                if counter2 == 0: 
                    quality_mask_or = condition_or
                else:
                    quality_mask_or = np.bitwise_or(quality_mask_or, condition_or)
                counter2+=1
            # Combine both masks
            quality_mask = np.bitwise_and(quality_mask, quality_mask_or)

    elif catalogName == 'V/147/sdss12':
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    elif catalogName == 'I/345/gaia2':
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    elif catalogName == 'I/284/out':
        # No mask yet
        quality_mask = np.ones(len(data), dtype=bool)

    return data[quality_mask]
