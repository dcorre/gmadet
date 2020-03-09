#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

"""
Group astronomical images by fields and epochs
"""

import errno, glob, os, subprocess, shutil
from astropy.io import fits
from astropy.table import Table
import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import time, wcs
import numpy as np

def rm_p(src):
   try:
      #shutil.rmtree(src, ignore_errors=True)
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

def table_obs(path_data, radius, deltaT):
   """ Create astropy table to group epochs and fields """

   
   # List of all raw files
   filenames = glob.glob(path_data+'*.fit*')

   names = []
   RA = []
   Dec = []
   Time = []
   telescopes = []
   instruments = []
   filters = []

   for ima in filenames:
        #print("processing " + ima + " ...\x1b[2K", end='\r', flush=True),
        hdr = fits.open(ima, memmap=False)[0].header
        # Get time of observation in hours
        try:
            date = time.Time(hdr['DATE-OBS'], format='fits')
            # convert Julian day in hours
            hr = date.jd * 24
        except:
            try:
                hr = float(hdr["JD"]) * 24.0
            except:
                print ('No keyword is found for the date of observation.\nExpected: `DATE-OBS` or `JD`')
        w = wcs.WCS(hdr)

        names.append(ima)
        RA.append(w.wcs.crval[0])
        Dec.append(w.wcs.crval[1])
        Time.append(hr)
        telescopes.append(hdr['TELESCOP'])
        instruments.append(hdr['INSTRUME'])
        filters.append(hdr['FILTER'])

   # Add unique index identifier per image
   idx = np.arange(len(names))
   # id to identify same field of view within given radius
   field_id = np.zeros(len(names), dtype=int)
   # id to identify epoch of same field within given time
   epoch_id = np.zeros(len(names), dtype=int)
   # Column to indicate the name of the stacked image 
   stack_name = [None] * len(names)
   # RA and Dec took as reference for one field
   ra_ref = [None] * len(names)
   dec_ref = [None] * len(names)

   obs_table = Table([idx, names, telescopes, instruments, filters, RA, Dec, Time, field_id, epoch_id, ra_ref, dec_ref, stack_name], names=['idx', 'filename', 'Telescope', 'Instrument', 'Filter', 'RA', 'Dec', 'JD', 'fieldID', 'epochID', 'RA_ref', 'Dec_ref', 'stack_filename'])

   # Sort by obs-time
   obs_table.sort('JD')

   field_id = 0
   for tel, inst, filt in obs_table.group_by(['Telescope', 'Instrument', 'Filter']).groups.keys:
       mask = (obs_table['Telescope'] == tel) & (obs_table['Instrument'] == inst) & (obs_table['Filter'] == filt)
       # Group by field of view
       # initialise with first image data
       ccrval_ref = SkyCoord(obs_table[mask]['RA'][0], obs_table[mask]['Dec'][0], unit=(u.deg, u.deg), frame='icrs')
       field_id = 1
       mask_idx = obs_table['idx'] == obs_table[mask]['idx'][0]
       obs_table['fieldID'][mask_idx] = field_id
       obs_table['RA_ref'][mask_idx] = obs_table[mask]['RA'][0]
       obs_table['Dec_ref'][mask_idx] = obs_table[mask]['Dec'][0]

       for data in obs_table[mask]:
           if data['fieldID'] == 0:
               # If image has not been associated to a field yet
               # Check for the closest field 
               # otherwise create new field ID
               ccrval = SkyCoord(data['RA'], data['Dec'], unit=(u.deg, u.deg), frame='icrs')
               mask2 = (obs_table['fieldID'] != 0) & mask
               sep_min = 100 # in degrees
               field_ref = -1
               for j, key in enumerate(obs_table[mask2].group_by('fieldID').groups.keys):
                   # Assume that ra and dec of one field is defined by first image for that field
                   mask3 = (obs_table['fieldID'] == key[0]) & mask2
                   ra_ref = np.atleast_1d(obs_table[mask3]['RA'])[0]
                   dec_ref = np.atleast_1d(obs_table[mask3]['Dec'])[0]
                   ccrval_ref = SkyCoord(ra_ref, dec_ref, unit=(u.deg, u.deg), frame='icrs')
                   sep =  ccrval.separation(ccrval_ref).degree
                   if (sep < radius) & (sep < sep_min):
                       sep_min = sep
                       field_ref = key[0]

               if field_ref != -1:
                   mask_idx = obs_table['idx'] == data['idx'] 
                   obs_table['fieldID'][mask_idx] = field_ref
                   obs_table['RA_ref'][mask_idx] = ra_ref
                   obs_table['Dec_ref'][mask_idx] = dec_ref
               else:
                   field_id += 1
                   mask_idx = obs_table['idx'] == data['idx'] 
                   obs_table['fieldID'][mask_idx] = field_id
                   obs_table['RA_ref'][mask_idx] = data['RA']
                   obs_table['Dec_ref'][mask_idx] = data['Dec']

   # Group fields by epochs
   for tel, inst, filt in obs_table.group_by(['Telescope', 'Instrument', 'Filter']).groups.keys:
       mask = (obs_table['Telescope'] == tel) & (obs_table['Instrument'] == inst) & (obs_table['Filter'] == filt)

       for field_id in obs_table[mask].group_by('fieldID').groups.keys:
           mask_field = (obs_table['fieldID'] == field_id[0]) & mask
           JD_ref = obs_table[mask_field]['JD'][0]
           epoch_id = 1
           for data in obs_table[mask_field]:
               if data['JD'] <= JD_ref + deltaT:
                   mask_idx = obs_table['idx'] == data['idx']
                   obs_table['epochID'][mask_idx] = epoch_id
               else:
                   epoch_id += 1
                   JD_ref = data['JD']
                   mask_idx = obs_table['idx'] == data['idx']
                   obs_table['epochID'][mask_idx] = epoch_id
   obs_table.show_in_browser()
   return obs_table

def makelists(path_data, radius, deltaT):
    """
    Group images by fields and epochs
    
    Parameters
    ----------
    path_data : path to images, string
        directory path to loop through all the fits file it contains
    path_lists : path to folder containing list of grouped images, string
        
    radius : radius in arcmin, float
        radius in arcmin used to group fields based on CRVAL values

    deltaT : time in hours, float
        maximum time interval for one epoch, i.e. from one image taken at
        time t, all images of the same field taken before t + deltaT
        are stacked


    Returns
    -------
    No variable is returned. Files containing the images to stack are
    created in stacklists/ folder


    """

    path_lists = path_data + 'fieldlists/'

    #Create folder for lists, delete existing files
    rm_p(path_lists)
    mkdir_p(path_lists)

    # Convert radius in degrees
    radius = radius / 60
    # Create observation table with images grouped by field and epoch
    fields = table_obs(path_data, radius, deltaT)

    # Create ascii files containing images to stack.
    # These files are the input of SWARP
    l = open(path_lists + "fields.slist", "w")
    for tel, inst, filt in fields.group_by(['Telescope', 'Instrument', 'Filter']).groups.keys:
        mask = (fields['Telescope'] == tel) & (fields['Instrument'] == inst) & (fields['Filter'] == filt)

        for field_id, epoch_id in fields[mask].group_by(['fieldID', 'epochID']).groups.keys:
            mask_field = (fields['fieldID'] == field_id) & (fields['epochID'] == epoch_id) & mask
            tel = str(fields['Telescope'][mask_field][0]).replace(' ','')
            band = str(fields['Filter'][mask_field][0]).replace(' ', '')
            ra = str(np.round(fields['RA_ref'][mask_field][0],3)).replace('.','')
            dec = str(np.round(fields['Dec_ref'][mask_field][0],3)).replace('.','')
            filename = tel + '_' + band + '_' + ra + '_' + dec + "_field_%03d_%03d" % (field_id, epoch_id) 
            #filename = prefix + "_%03d_%03d" % (field_id, epoch_id) 
            f = open(path_lists + filename + ".list", "w")
            for data in fields[mask_field]:
                f.write(data['filename'] + "\n")
                mask_idx = fields['idx'] == data['idx']
                fields['stack_filename'][mask_idx] = filename
            f.close()
            l.write(filename + " ")
    l.close()
    #fields.show_in_browser()

def stacking(path_data, radius, deltaT,useweight=False, subBack=True, gain=1):
    """Stack images"""

    # Add '/' at the end of the paths if they are missing
    if path_data[-1] != '/':
        path_data = path_data + '/'

    path_lists = path_data + 'fieldlists/'
    path_stacks = path_data + 'stacks/'

    useweight = bool(useweight)

    # Whether to substrack background
    if subBack:
        subBack = 'Y'
    else:
        subBack = 'N'

    # Make list of images to stack
    makelists(path_data, radius, deltaT)

    mkdir_p(path_stacks)

    # Get all the prefixes corresponding to one field
    filenames = glob.glob(path_lists + '*.list')
    print (filenames)
    prefixes = []
    for filename in filenames:
        splitfilename = os.path.splitext(filename)[0].split('/')[-1].split('_')
        prefi = ''
        for i in range(len(splitfilename)-1):
            prefi += splitfilename[i] + '_'
        prefixes.append(prefi)
    # Discard duplicates
    prefixes = np.unique(prefixes)

    # Loop over fields
    for pref in prefixes:
        imalists = []
        epochs = []
        # Loop over epochs
        for imalist in glob.glob(path_lists + pref + '???.list'):
            # Check that there are at least 2 images to stack
            # Otherwise skip it
            file = np.genfromtxt(imalist, dtype=str)
            if len(np.atleast_1d(file)) < 2:
                continue

            epochs += [path_stacks +  os.path.splitext(imalist)[0].split('/')[-1]]
            imalists += ['@' + imalist]

        point = path_stacks + pref
        subprocess.call(['swarp', '-HEADER_ONLY', 'Y', '-IMAGEOUT_NAME', \
                            point + '.head', '-GAIN_DEFAULT', str(gain)] + imalists)
        subprocess.call(['sed', '-i', \
                             's/MJD-OBS/COMMENT/; s/EXPTIME/COMMENT/; s/GAIN   /COMMENT/; s/SATURATE/COMMENT /', \
                         point + '.head'])

        for i, imalist in enumerate(imalists):
            epoch = epochs[i]
            shutil.copy(point + '.head', epoch + '.head')
            if useweight:
                subprocess.call(['swarp',
                                 '-IMAGEOUT_NAME', epoch + '.fits', \
                                 '-SUBTRACT_BACK', subBack,\
                                 '-BACK_SIZE', '128', \
                                 '-BACK_FILTERSIZE', '3',\
                                 '-WEIGHTOUT_NAME', epoch + '.weight.fits', \
                                 '-RESAMPLING_TYPE', 'LANCZOS3',\
                                 '-OVERSAMPLING', '0',\
                                 '-COMBINE_TYPE', 'MEDIAN',\
                                 '-GAIN_DEFAULT', str(gain)] + [imalist])
            else:
                subprocess.call(['swarp',
                                 '-IMAGEOUT_NAME', epoch + '.fits',\
                                 '-GAIN_DEFAULT', str(gain),\
                                 '-SUBTRACT_BACK', subBack,\
                                 '-BACK_SIZE', '128', \
                                 '-BACK_FILTERSIZE', '3',\
                                 '-RESAMPLING_TYPE', 'LANCZOS3',\
                                 '-OVERSAMPLING', '0',\
                                 '-COMBINE_TYPE', 'MEDIAN'] + [imalist])
            rm_p(epoch + '.head')

        rm_p(point + '.head')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Rebin astronomical images.')
   
    parser.add_argument('--path_data',
                        dest='path_data',
                        required=True,
                        type=str,
                        help='Path where the files to be stacked are.')

    parser.add_argument('--radius',
                       dest='radius',
                       required=False,
                       type=float,
                       default=5,
                       help='Radius in arcmin to group fields. Default: 5 arcmin')
  
    parser.add_argument('--deltaT',
                       dest='deltaT',
                       required=False,
                       type=float,
                       default=1,
                       help='Time interval in hours to group fields into same epoch. Default: 1h')


    parser.add_argument('--no_BackSub',
                       dest='no_BackSub',
                       action='store_false',
                       help='If provided as argument, no background substraction is performed on each image prior to stacking. Default: Substraction is performed')


    args = parser.parse_args()
    #path_data = '/home/corre/codes/gmadet/gmadet/data/Atlascow/Iris10obs/rebinned_images'
    stacking(args.path_data,args.radius,args.deltaT,subBack=args.no_BackSub)



