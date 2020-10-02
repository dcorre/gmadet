#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test usage of hips2fits"""

from astropy.io import fits
from urllib.parse import quote

class hips2fits:

    def __init__(self, fov, ra, dec, pixel_scale, catalog, band):

        width = int(fov * 3600. / pixel_scale)
        height = int(fov * 3600. / pixel_scale)

        if catalog == 'Pan-STARRS':
            # for Pan-STARRS, available bands: g, r,i, z, y, color-i-r-g, color-z-zg-g
            hips = 'PanSTARRS/DR1/' + band

        elif catalog == 'DECaLS':
            # for DecaLS, only g, r and color are available
            hips = 'DECaLS/DR5/' + band

        self.fov = fov
        self.ra = ra
        self.dec = dec
        self.pixel_scale = pixel_scale
        self.band = band
        self.hips = hips
        self.width = width
        self.height = height

    def download(self, path):
        url = 'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={}' \
              '&width={}&height={}&fov={}&projection=TAN' \
              '&coordsys=icrs&ra={}&dec={}'.format(quote(self.hips),
                                                   self.width,
                                                   self.height,
                                                   self.fov,
                                                   self.ra,
                                                   self.dec)

        hdu = fits.open(url)
        file_name = path + '/{}-{}-{}-{}-{}-{}.fits'.format(self.hips.replace('/', '_'),
                                                            self.fov,
                                                            self.ra,
                                                            self.dec,
                                                            self.pixel_scale,
                                                            self.band)
        hdu.writeto(file_name, overwrite=True)


    def sources_detection():
        pass

    def crossmatch():
        pass

    def photometric_calibration():
        """
        Photometrically calibrate the hips2fits image.
        Compute the zeropoint with respect to the catalog
        and apply the correction on the image.
        """
        pass


if __name__ == "__main__":

    catalog = 'Pan-STARRS'
    catalog = 'DECaLS'

    fov = 0.2
    pixel_scale = 0.8
    ra = 244.00
    dec = 22.28
    band='g'

    test = hips2fits(fov, ra, dec, pixel_scale, catalog, band)
    test.download('.')
