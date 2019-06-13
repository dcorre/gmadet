#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord


def make_sub_image(filename, OT_coords, coords_type='world',
                   output_name='subimage.fits.gz', size=[200,200]):
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
        ra, dec = w.all_pix2world(np.array(pix), 0)
        pix_ref = [ra,dec]

    # Extract subimage from image starting from reference pixel
    subimage = data[int(pix[0]) - int(size[0]/2) : int(pix[0]) + int(size[0]/2),
                    int(pix[1]) - int(size[1]/2) : int(pix[1]) + int(size[1]/2)]
    
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

if __name__ == '__main__':

    # Set the path to your image, including the extension 
    filename = 'yourimage.fits'
    # Set coordinates to use as center of new sub-image
    OT_coords = [190.885035, 12.742086]
    # Type of coordinates: 'world' or 'pix' if wcs or physical resp.
    coords_type='world'
    # Set size of new sub-image (in pixels)
    size = [200,200]
    # Path to new subimage, including extension (you can let .fits.gz by default)
    output_name = 'subimage.fits.gz'

    make_sub_image(filename, OT_coords, coords_type=coords_type, size=size, output_name=output_name)


