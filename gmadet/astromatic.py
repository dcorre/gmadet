#! /usr/bin/env python
# -*- coding: utf-8 -*-
  
"""Detection of sources."""

import errno, glob, os, shutil, subprocess, sys
import numpy as np
from astropy.io import fits


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


def scamp(filename, useweight=False):
    """Compute PSF in astronomical images"""
    
    #imagelist=glob.glob(path+'/*.fits')
    path = os.path.dirname(filename)
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
         root = os.path.splitext(ima)[0]
         print ('Sextractor')
         #print("Processing " + ima + " ...", end='\r', flush=True),
         subprocess.call(['sex', '-c', 'prepscamp.sex', ima])

         cat = 'prepscamp.cat'
         print ('scamp')
         subprocess.call(['scamp', cat, '-c', 'scamp.conf'])
         mv_p('snap_prepsfex.fits', root + '.psf.fits')

    # Adding header
    hdulist=fits.open(filename)
    
    # First remove old ones
    keywords_to_remove = ["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CD1_1","CD1_2","CD2_1","CD2_2","CDELT1", "CDELT2","PIXSCALX","PIXSCALY","CUNIT1","CUNIT2","WCSAXES","WCSNAME","RADESYS","WCSVERS","CTYPE1","CTYPE2", "EQUINOX"]

    for keyword in keywords_to_remove:
        if keyword in hdulist[0].header:
                del hdulist[0].header[keyword]
    s=''
    with open('prepscamp.head') as f: 
        for line in f: 
            s = s + line + '\n' 
    newheader = fits.Header.fromstring(s, sep='\n')
    for key, value in newheader.items():
        print (key, value)
        #truncate long keys
        if len(key) > 8:
            key = key[:7]
        try:
            hdulist[0].header.set(key.upper(), value)
        except:
            try:
                hdulist[0].header.set(key.upper(), str(value))
            except:
                pass
    
    hdulist.writeto(path+'/test_scamp.fits',overwrite=True)



def psfex(filename, useweight=False):
    """Compute PSF in astronomical images"""
    
    #imagelist=glob.glob(path+'/*.fits')
    imagelist = np.atleast_1d(filename)
    print (imagelist)
    for ima in imagelist:
         root = os.path.splitext(ima)[0]
         if useweight:
             weight = root + '.weight.fits'
             #print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima, '-WEIGHT_IMAGE', weight])
         else:
             print ('Sextractor')
             #print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima])

         cat = 'prepsfex.cat'
         print ('psfex')
         #subprocess.call(['psfex', cat, '-c', 'config.psfex'])
         subprocess.call(['psfex', cat])
         mv_p('snap_prepsfex.fits', root + '.psf.fits')

def sextractor(path, outdir='sources_cat/', useweight=False):
    """Detect sources in astronomical images"""
    
    imagelist=glob.glob(path+'/*.fits')
    mkdir_p(outdir)

    for ima in imagelist:
        if '.psf' not in ima:
            root = os.path.splitext(ima)[0]
            if useweight:
                weight = root + '.weight.fits'
                #print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima, '-WEIGHT_IMAGE', weight])
            else:
                print (ima.split('/')[-1].split('.')[0])
                #print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima,
                                 '-CATALOG_NAME', outdir+ima.split('/')[-1].split('.')[0]+'.cat'])

