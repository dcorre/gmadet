#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# @author: David Corre (IJCLab/CNRS)

import sys

import errno, glob, os
import numpy as np
from astropy.io import fits
import argparse
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astroquery.xmatch import XMatch
from astroquery.imcce import Skybot
from astroML.crossmatch import crossmatch_angular
from astropy.table import Table

def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise

def get_candidates_wcs(telescope, path_events):
    """ Get wcs coordinates of detected candidates  """

    dir_data = path_events + telescope + '/'

    # Get all the events images filenames
    filenames = glob.glob(dir_data + '*.fits')
    ra_list = []
    dec_list = []
    index = []
    filename_list = []

    for i, filename in enumerate(filenames):
        hdu = fits.open(filename, memmap=False)
        header = hdu[0].header
        data1 = hdu[0].data[0].astype(np.float32)
        data3 = hdu[0].data[2].astype(np.float32)
        xpos, ypos = np.unravel_index(data3.argmax(), data3.shape)
        # required to mimic that there is only one array when converting pix2world
        header['NAXIS'] = 1

        w = WCS(header)
        ra, dec = w.all_pix2world(xpos, ypos, 0)
        ra_list.append(ra)
        dec_list.append(dec)
        index.append(i)
        filename_list.append(filename)

    candidates = Table([index, ra_list, dec_list, filename_list], names=['id', '_RAJ2000', '_DEJ2000', 'filename'])

    return candidates

def xmatch(coordinates, catalog, radius):
    """
    Perform cross-match with a catalog using the CDS XMatch 
    parameters: coordinates, catalog, radius:
                coordinates: astropy table with RA, DEC of all detected sources
                catalog: Vizier identifier of the catalog
                radius in arcsecond
    returns: astropy.table object

    Vizier catalog identifiers:
    Gaia DR2: I/345/gaia2
    SDSS DR12: V/147/sdss12
    2MASS: II/246/out
    USNO B1: I/284/out
    USNO A2: I/252/out
    GLADE 2: VII/281/glade2
    Panstarrs DR1: II/349/ps1
    VSX: B/vsx/vsx
    """

    matched_stars = XMatch.query(coordinates,
                     cat2='vizier:%s' % catalog,
                     max_distance=radius * u.arcsec, colRA1='_RAJ2000',
                     colDec1='_DEJ2000')

    return matched_stars




def crossmatch(telescope, path_events, radius, catalog='B/vsx/vsx'):
    """Radius in arcseconds """

    dir_data = path_events + telescope + '/'

    candidates = get_candidates_wcs(telescope, path_events)

    # Save candidates coordinates in ascii file
    candidates.write(dir_data+'candidates.dat', format='ascii.commented_header', overwrite=True)
   
    candidates['_RAJ2000'] *= u.deg
    candidates['_DEJ2000'] *= u.deg

    # perform crossmatch with variable star catalog
    crossmatch = xmatch(candidates, catalog, radius)

    print (crossmatch)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Apply a previously trained CNN.')

    parser.add_argument('--path_events',
                        dest='path_events',
                        required=True,
                        type=str,
                        help='Path to the detected events')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')


    parser.add_argument('--radius',
                        dest='radius',
                        required=True,
                        type=float,
                        help='Crossmatch radius with catalogues')

    args = parser.parse_args()

    crossmatch(args.telescope, args.path_events, args.radius)

