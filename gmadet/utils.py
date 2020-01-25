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
                'conf': path2tel+'sourcesdet.sex',
                'param': path2tel+'sourcesdet.param',
                'default_conv': path2tel+'default.conv'
                },
            'scamp': {
                'sextractor': path2tel+'prepscamp.sex',
                'param': path2tel+'prepscamp.param',
                'conf': path2tel+'scamp.conf',
                },
            'swarp': {

                },
            'psfex': {

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
                crpix1 = int((y2-y1)/2)
                crpix2 = int((x2-x1)/2)
                newheader['CRPIX1'] = crpix1
                newheader['CRPIX2'] = crpix2
                w = wcs.WCS(header)
                ra, dec = w.wcs_pix2world(crpix1+Naxis11*i, crpix2+Naxis22*j, 1)
                newheader['CRVAL1'] = str(ra)
                newheader['CRVAL2'] = str(dec)
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


def send_data2DB(filename, candidates, Nb_cuts, owncloud_path, VOE_path, usrpwd_path, FoV=60, coords_type='pix', corner_cut=32,debug=False,fmt='png'):
    """Send candidates information to database"""

    # Load data and header
    data, header = fits.getdata(filename, header = True)
    dateObs = header['DATE-OBS']
    Tstart = Time(dateObs, format='fits', scale='utc')
    Tend = Tstart + TimeDelta(float(header['EXPOSURE']), format='sec')
    
    # Try to get airmass rom header, else set it to -1
    try:
        Airmass = header['AIRMASS']
    except:
        Airmass = -1
    
    # Do not consider candidates found in the image edge
    imsize = data.shape
    print (imsize, header['NAXIS1'], header['NAXIS2'])
    # Get the physical pixels of the original size if image were split into different quadrants.
    for i, candidate in enumerate(candidates):
        quadrant_idx = candidate['quadrant']
        quadrant, index_i, index_j = quadrant_idx.split('_')
        quadrant = quadrant[1:]
        
        candidates['Xpos'][i] = candidate['Xpos'] + int(imsize[0]/Nb_cuts[0]) * int(index_j)
        candidates['Ypos'][i] = candidate['Ypos'] + int(imsize[1]/Nb_cuts[1]) * int(index_i)
    print (candidates)
    candidates.write('test_all.oc', format='ascii.commented_header', overwrite=True)
    mask = (candidates['Xpos'] > corner_cut) & \
            (candidates['Ypos'] > corner_cut) & \
            (candidates['Xpos'] < imsize[1] - corner_cut) & \
            (candidates['Ypos'] < imsize[0] - corner_cut)
    candidates_cut = candidates[mask]
    print (candidates_cut)
    # Get information about the current alert from the xml file containing observation plan
    with open(VOE_path,'rb') as f:
        obsplan = vp.load(f)

    dict_event={}
    dict_event['event_type'] = obsplan.find(".//Param[@name='Event_type']").attrib['value']
    dict_event['event_name'] = obsplan.find(".//Param[@name='Event_ID']").attrib['value']
    dict_event['event_status'] = obsplan.find(".//Param[@name='Event_status']").attrib['value']
    dict_event['revision'] = obsplan.find(".//Param[@name='Revision']").attrib['value']
    dict_event['telescope'] = obsplan.find(".//Param[@name='Name_tel']").attrib['value']

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
    for i, row in enumerate(candidates_cut):
        name = dict_event['telescope'] + '_' + \
                str(round(float(row['_RAJ2000']),5)) + '_' + \
                str(round(float(row['_DEJ2000']),5)) + '_' + \
                dateObs + '.' + fmt
        name = 'test' + str(i) + '.' + fmt
        Fits_path.append(name)
        if coords_type == 'world':
            OT_coords = [row['_RAJ2000'], row['_DEJ2000']]
        elif coords_type == 'pix':
            OT_coords = [row['Ypos'], row['Xpos']]
        make_sub_image(filename, OT_coords, coords_type=coords_type,
                       output_name=outputDir+name, size=[128,128], FoV=FoV, fmt=fmt)

    alias = ['new'] * len(candidates_cut)
    new = [1] * len(candidates_cut)
    tile_id_list = [3] * len(candidates_cut)
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
    #url = "http://localhost/test2/obs_report_OT.php"

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

        if debug:
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




def get_filter(filename, telescope):
    """ Get the name of the filter band from header and telescope name
        And associate the correct name from DB
    """
    
    header = fits.getheader(filename)

    band = header['FILTER']

    if telescope == 'TRE':
        if band in  ['C', 'Clear', 'NoFilter']:
            band_DB = 'Clear'
            band_cat = 'g+r'
    elif telescope == 'Lisniky-AZT8':
        if band in ['R']:
            band_DB = 'R/Johnson'
            band_cat = 'r'
    elif telescope == 'UBAI-T60S':
        if band in ['R']:
            band_DB = 'R/Johnson'
            band_cat = 'r'
    elif telescope == 'UBAI-T60N':
        if band in ['R']:
            band_DB = 'R/Johnson'
            band_cat = 'r'
    # Need to add all the filters from telescopes


    return band_DB, band_cat
    
       


