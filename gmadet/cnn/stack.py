#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, os, subprocess, shutil
from astropy.io import ascii, fits
import argparse
import numpy as np

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


def stacking(telescope, path_lists, path_stacks, useweight, nb_stacks, gain):
    """Stack images"""
    
    useweight = bool(useweight)

    dir_lists = path_lists + telescope + '/'
    dir_stacks = path_stacks + telescope + '/'

    mkdir_p(dir_stacks)

    # Get all the prefixes corresponding to one field
    filenames = glob.glob(dir_lists + '*.list')
    prefixes = []
    for filename in filenames:
        splitfilename = os.path.splitext(filename)[0].split('/')[-1].split('_')
        prefixes.append(splitfilename[0] + '_' + splitfilename[1])
    # Discard duplicates
    prefixes = np.unique(prefixes)
    
    # Loop over fields
    for prefix in prefixes:
        imalists = []
        epochs = []
        # Loop over epochs
        for imalist in glob.glob(dir_lists + prefix + '_??.list'):
            # Check that there are the right number of images to stack
            # Otherwise skip it
            file = np.genfromtxt(imalist, dtype=str)
            if len(np.atleast_1d(file)) < nb_stacks:
                continue

            epochs += [dir_stacks +  os.path.splitext(imalist)[0].split('/')[-1]]
            imalists += ['@' + imalist]

        point = dir_stacks + prefix
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
                                 '-WEIGHTOUT_NAME', epoch + '.weight.fits', \
                                 '-GAIN_DEFAULT', str(gain)] + [imalist])
            else:
                subprocess.call(['swarp',
	                             '-IMAGEOUT_NAME', epoch + '.fits',\
                                 '-GAIN_DEFAULT', str(gain),\
                                 '-BACK_SIZE', '128', \
                                 '-BACK_FILTERSIZE', '3',\
                                 '-RESAMPLING_TYPE', 'LANCZOS3',\
                                 '-OVERSAMPLING', '0',\
                                 '-COMBINE_TYPE', 'MEDIAN'] + [imalist])
            rm_p(epoch + '.head')
            
        rm_p(point + '.head')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Stacking images.')

    parser.add_argument('--path_lists',
                        dest='path_lists',
                        required=True,
                        type=str,
                        help='Path to the raw data')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--path_stacks',
                        dest='path_stacks',
                        required=True,
                        type=str,
                        help='Path where to save stacked images')

    parser.add_argument('--useweight',
                        dest='useweight',
                        required=True,
                        type=int,
                        help='Use weight map. 0:no / 1:yes')

    parser.add_argument('--nb_stacks',
                        dest='nb_stacks',
                        required=True,
                        type=int,
                        help='Number of images to stack')

    parser.add_argument('--gain',
                        dest='gain',
                        required=True,
                        type=float,
                        help='Gain in e-/ADU')

    args = parser.parse_args()
    
    stacking(args.telescope, args.path_lists, args.path_stacks, args.useweight, args.nb_stacks, args.gain)

