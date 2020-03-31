# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:36:04 2019

@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, math, os, shutil, subprocess, sys
import argparse

def mv_p(src, dest):
  try:
    shutil.move(src, dest)
  except:
    pass

def compute_psf(telescope, path_stacks, useweight):
    """ Compute psf cube for each stacked images """

    #path to data
    dir = path_stacks + telescope + '/'  
   
    imagelist = sorted(glob.glob(dir + '*.fits'))

    useweight = bool(useweight)

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
         subprocess.call(['psfex', cat, '-NTHREADS 4'])
         mv_p('snap_prepsfex.fits', root + '.psf.fits')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Formatting data before stacking.')

    parser.add_argument('--path_stacks',
                        dest='path_stacks',
                        required=True,
                        type=str,
                        help='Path to the stacked images')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--useweight',
                        dest='useweight',
                        required=True,
                        type=int,
                        help='Use weight map. 0:no / 1:yes')


    args = parser.parse_args()

    compute_psf(args.telescope, args.path_stacks, useweight=args.useweight)

