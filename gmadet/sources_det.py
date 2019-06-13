#! /usr/bin/env python
# -*- coding: utf-8 -*-
  
"""Detection of sources."""

import errno, glob, os, shutil, subprocess, sys

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



def psf(path, soft="sextractor", useweight=False):
    """Compute PSF in astronomical images"""
    
    imagelist=glob.glob(path+'/*.fits')
    print (imagelist)
    for ima in imagelist:
         root = os.path.splitext(ima)[0]
         if useweight:
             weight = root + '.weight.fits'
             print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima, '-WEIGHT_IMAGE', weight])
         else:
             print("Processing " + ima + " ...", end='\r', flush=True),
             subprocess.call(['sex', '-c', 'prepsfex.sex', ima])

         cat = 'prepsfex.cat'
         subprocess.call(['psfex', cat])
         mv_p('snap_prepsfex.fits', root + '.psf.fits')

def sources_det(path, outdir='sources_cat/', soft="sextractor", useweight=False):
    """Detect sources in astronomical images"""
    
    imagelist=glob.glob(path+'/*.fits')
    mkdir_p(outdir)

    for ima in imagelist:
        if '.psf' not in ima:
            root = os.path.splitext(ima)[0]
            if useweight:
                weight = root + '.weight.fits'
                print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima, '-WEIGHT_IMAGE', weight])
            else:
                print (ima.split('/')[-1].split('.')[0])
                print("Processing " + ima + " ...", end='\r', flush=True),
                subprocess.call(['sex', '-c', 'sourcesdet.sex', ima,
                                 '-CATALOG_NAME', outdir+ima.split('/')[-1].split('.')[0]+'.cat'])

