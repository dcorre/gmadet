#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @file		convert.py
# @brief	Convert cutouts to datacube
# @date		12/11/2018
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	P9 search scripts
#
#	Copyright:		(C) 2018 IAP/CNRS/SorbonneU
#
#	Author:			Emmanuel Bertin (IAP)
#
#	License:		GNU General Public License
#
#	Bertinoscopic is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	Bertinoscopic is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with Bertinoscopic. If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Original scrip modified by: David Corre (IJCLab/CNRS)

import sys

#if len(sys.argv) < 3:
#  print('syntax: ' + sys.argv[0] + ' <output_npz> <stack_prefix_1> [<stack_prefix2> ...]')
#  sys.exit(0)

import glob, re, os, errno
import numpy as np
from astropy.io import fits
import argparse
from gmadet.utils import getpath, rm_p, mkdir_p

def convert(path, telescope, cubename):
    """convert simulated data before starting training"""
    
    path_gmadet = getpath()

    outdir = path_gmadet + '/cnn/datacube/' + telescope + '/'
    mkdir_p(outdir)

    # Get all the prefixes corresponding to one field
    truelist = glob.glob(path + '/true/*.fits')
    falselist = glob.glob(path + '/false/*.fits')
    # output cube name
    npz_name = "%s.npz" % cubename

    Ncand = len(truelist)+len(falselist)
    cube = [] #np.zeros((Ncand, 64, 64))
    labels = []
    mags = []
    errmags = []
    filters = []
    for cand in truelist:
        hdus = fits.open(cand, memmap=False)
        head = hdus[0].header
        # Exclude cases too close to the edge
        # Meaning they are located at less than the defined size
        # of the small images
        if head['EDGE'] == 'False':  
            labels += [1]
            mags += [head['MAG']]
            errmags += [head['MAGERR']]
            filters += [head['FILTER']]
            cube.append(hdus[0].data)
        hdus.close()

    for cand in falselist:
        hdus = fits.open(cand, memmap=False)
        head = hdus[0].header
        #if hdus[0].data.shape != (64, 64):
        #    print ('skip %s as its shape is not (64,64): (%d,%d)' % (cand, hdus[0].data.shape[0], hdus[0].data.shape[1]))
        #    continue
        # Exclude cases too close to the edge
        # Meaning they are located at less than the defined size
        # of the small images
        if head['EDGE'] == 'False'
            labels += [0]
            mags += [head['MAG']]
            errmags += [head['MAGERR']]
            filters += [head['FILTER']]
            cube.append(hdus[0].data)
        hdus.close()


    print("Converting and reshaping arrays ...")
    # Convert lists to B.I.P. NumPy arrays
    print (len(cube))
    # Check whether all candidates has 64x64 pixels
    # If not, delete them
    # This can happen at the edge of images
    #for i in range(len(cube)):
    #    if np.array(cube[i]).shape != (64, 64):
    #        print (i, np.array(cube[i]).shape)
    #        del cube[i]
    print (len(cube))
    cube = np.asarray(cube, dtype=np.float32)
    print (cube.shape)
    if cube.ndim < 4:
        cube = np.reshape(cube, [cube.shape[0], cube.shape[1], cube.shape[2], 1])
    else:
        cube = np.moveaxis(cube, 1,-1)

    # Report dimensions of the data cube
    print("Saving %d %d×%d×%d image datacube ..." %cube.shape, end='\r', flush=True)
    np.savez(outdir+npz_name, cube=cube, labels=labels, mags=mags, errmags=errmags, filters=filters)

    print("Saved to " + outdir + npz_name)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Convert sim data to a .npz cube.')

    parser.add_argument('--path',
                        dest='path',
                        required=True,
                        type=str,
                        help='Path to simulated data')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--cubename',
                        dest='cubename',
                        required=True,
                        type=str,
                        help='Name of the datacube')

    args = parser.parse_args()

    convert(args.path,args.telescope,args.cubename)

