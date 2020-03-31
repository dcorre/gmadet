#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Group by fields and epochs

@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, os, subprocess, shutil
from astropy.io import fits
from astropy.table import Table
import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs

def rm_p(src):
  try:
    shutil.rmtree(src, ignore_errors=True)
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


def makelists(telescope, path_data, path_lists, prefix, radius):
    """ Group images by fields and epochs"""

    #path to data
    dir = path_data + telescope + '/'

    #Create folder for lists, delete existing files
    rm_p(path_lists + telescope + '/')
    mkdir_p(path_lists + telescope + '/')

    # Convert radius in degrees
    radius = radius / 60
    
    # List of all raw files
    expos = sorted(glob.glob(dir+'*.fits'))

    objname = []
    n = {}
    crvals = {}
    objo = ""
    hro = 0
    l = open(path_lists + telescope + '/' + prefix + ".slist", "w")
    first=True
    for expo in expos:
        print("processing " + expo + " ...\x1b[2K", end='\r', flush=True),
        hdr = fits.open(expo, memmap=False)[0].header

        if radius <= 0:
            obj = hdr["OBJECT"].replace("-", "_")
        else:
            w = wcs.WCS(hdr)
            ccrval = SkyCoord(w.wcs.crval[0], w.wcs.crval[1], unit=(u.deg, u.deg), frame='icrs')
            ind = 1
            for crval in crvals:
                if ccrval.separation(crvals[crval]).degree < radius:
                    break
                ind += 1
            obj = prefix + '_%03d' %ind
            crvals[obj] = ccrval
        hr = float(hdr["MJD-OBS"]) * 24.0
        # New epoch if exposure taken more than 30 minutes after first one
        # or coordinates not matching with any previous pointing
        if obj != objo or hr > hro + 0.5:
            if obj in n:
                n[obj] += 1
            else:
                n[obj] = 1
                if first:
                    l.write(obj)
                else:
                    l.write(' ' + obj)
                first=False
            if "f" in locals():
                f.close()
            f = open(path_lists + telescope + '/' + obj + "_%02d" %(n[obj]) + ".list", "w")
        f.write(expo + "\n")
        objo = obj
        hro = hr
    f.close()
    l.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Group images by field and epochs.')

    parser.add_argument('--path_data',
                        dest='path_data',
                        required=True,
                        type=str,
                        help='Path to the raw data')

    parser.add_argument('--path_lists',
                        dest='path_lists',
                        required=True,
                        type=str,
                        help='Path where to save list of fields and epochs')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--prefix',
                        dest='prefix',
                        required=True,
                        type=str,
                        help='Prefix for files listing fields and epochs')

    parser.add_argument('--radius',
                        dest='radius',
                        required=True,
                        type=float,
                        help='Crossmatch radius in arcminutes to group field of different epochs')

    args = parser.parse_args()
    
    makelists(args.telescope, args.path_data, args.path_lists, args.prefix, args.radius)

