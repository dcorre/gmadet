#! /usr/bin/env python
# -*- coding: utf-8 -*-


import subprocess, sys, os
import numpy as np
from astropy.io import fits
from lacosmic import lacosmic
from utils import cp_p

def run_lacosmic(filename, FWHM, flim=2, sigma=5):
    """Run lacosmic to remove cosmic rays from the input image"""

    imagelist = np.atleast_1d(filename)

    for i, ima in enumerate(imagelist):
        path, filename_ext = os.path.split(ima)
        if path:
            folder = path + '/'
        else:
            folder = ''

        filename2 = filename_ext.split('.')[0]

        # Make copy of original image
        cp_p(ima, folder+filename2+'_CR_notcleaned.fits')

        hdulist = fits.open(ima)
        hdr = hdulist[0].header
        try:
            gain = hdr['GAIN']
        except:
            gain = 1
        try:
            RN = hdr['RN']
        except:
            RN = 10

        if FWHM[i] > 2:
            flim = 2
        else:
            flim = 5
        data = np.asarray(hdulist[0].data, dtype=float)
        lacosmic_res = lacosmic(data,flim,sigma,sigma,effective_gain=gain,readnoise=RN)

        # Create image cleaned from cosmic rays
        hdulist[0].data = lacosmic_res[0]
        hdulist.writeto(ima, overwrite=True)

        # Create mask of cosmic rays
        hdulist[0].data = np.asarray(lacosmic_res[1],dtype=int)
        hdulist.writeto(folder+filename2+'_CRmask.fits', overwrite=True)


def update_headers_scamp(filename, scamphead, pixelscale):
    """Modify the header after running scamp"""

    hdulist = fits.open(filename)
    # Verify and try to fix issue with fits standard
    hdulist.verify('fix')
    hdr = hdulist[0].header
    # First remove old keywords related to astrometric calibration
    keywords_to_remove = ["CRPIX1","CRPIX2","CRVAL1","CRVAL2","CD1_1","CD1_2","CD2_1","CD2_2","CDELT1", "CDELT2","PIXSCALX","PIXSCALY","CUNIT1","CUNIT2","WCSAXES","WCSNAME","RADESYS","WCSVERS","CTYPE1","CTYPE2", "EQUINOX", "COORDSYS", "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER", "IMAGEW", "IMAGEH", "LONPOLE", "LATPOLE", "CTYPE1T","CTYPE2T", "CRPIX1T","CRPIX2T","CRVAL1T","CRVAL2T", "CDELT1T", "CDELT2T", "CROTA1", "CROTA2", "CROTA1T", "CROTA2T", "AIRMASS"]

    # List of the first letters for standard astrometric coefficients.
    # Such as TR, PV, SIA
    coeff_astro = ['TR', 'SIA', 'A_', 'B_', 'AP_', 'BP_', 'LT', 'PV']

    for  key, value in hdr.items():
        for coeff in coeff_astro:
            _len = len(coeff)
            if key[:_len] in [coeff]:
                keywords_to_remove.append(key)

    for keyword in keywords_to_remove:
        if keyword in hdr:
                del hdr[keyword]

    s=''
    with open(scamphead) as f:
        for line in f:
            s = s + line + '\n'
    newheader = fits.Header.fromstring(s, sep='\n')
    for key, value in newheader.items():
        #print (key, value)
        #truncate long keys
        if len(key) > 8:
            key = key[:7]
        try:
            hdr.set(key.upper(), value)
        except:
            try:
                hdr.set(key.upper(), str(value))
            except:
                pass
    hdr.set('CDELT1', pixelscale[0])
    hdr.set('CDELT2', pixelscale[1])

    if (hdr['CTYPE1'] != 'RA---TPV') and (hdr['CTYPE2'] != 'DEC--TPV'):
        print ('\nWARNING: scamp did not set CTYPE1, CTYPE2 to RA---TPV and DEC--TPV.')
        print ('Set them to to these values by hand, otherwise it can not be read by astropy.wcs')
        print ('Likely due to some SIP distortions parameters already present in headers.')
        print ('One might check the astrometry to be safe.\n')
        hdr['CTYPE1'] = 'RA---TPV'
        hdr['CTYPE2'] = 'DEC--TPV'
    hdulist.writeto(filename,overwrite=True)


def astrometrynet(filename, radius=4, scaleLow=3.5, scaleHigh=4, scaleUnits='arcsecperpix'):
    """Run astrometry.net on the input image"""

    # Get the RA and DEC before removing the keywords as required by astrometry.net. 
    # Is it really required? need to check
    hdul = fits.open(filename)
    hdr = hdul[0].header

    ra = str(hdr['CRVAL1'])
    dec = str(hdr['CRVAL2'])

    del hdr['CRVAL1']
    del hdr['CRVAL2']
    del hdr['CRPIX1']
    del hdr['CRPIX2']
    hdul[0].header = hdr
    hdul.writeto(filename,overwrite=True)
    hdul.close()

    # Run astrometry.net
    subprocess.call(['solve-field', filename, '--ra', ra, '--dec', dec, '--radius', str(radius), '--scale-units', str(scaleUnits), '--scale-low', str(scaleLow), '--scale-high', str(scaleHigh), '--crpix-center', '--tweak-order','5','--uniformize','0', '--no-plots','--overwrite'])
    
    # Delete temporary files
    clean_tmp_files(filename, soft='astrometrynet')


def scamp(filename, config, useweight=False, CheckPlot=False, verbose='NORMAL'):
    """Compute astrometric solution of astronomical image using scamp"""

    path = os.path.dirname(filename)
    imagelist = np.atleast_1d(filename)
    for ima in imagelist:
         # Make sure to use only the Primary hdu.
         # sextractor, scamp seems to crash otherwise.
         hdul = fits.open(filename)
         #print (hdul.info())
         # Delete empty keywords
         for key, val in hdul[0].header.items():
             if key == '':
                  del hdul[0].header[key]
         newhdu = fits.PrimaryHDU()
         newhdu.data = hdul[0].data
         newhdu.header = hdul[0].header
         newhdulist = fits.HDUList([newhdu])
         newhdulist.writeto(filename,overwrite=True)
         hdul.close()

         root = os.path.splitext(ima)[0]
         _name = root.split('/')[-1]

         #print ('Create FITS-LDAC file from SExtractor')
         subprocess.call(['sex', '-c', config['scamp']['sextractor'], \
                 '-PARAMETERS_NAME', config['scamp']['param'], \
                 #'-FILTER_NAME', config['sextractor']['default_conv'], \
                 ima])

         if CheckPlot:
             plotnames = '%s_fgroups,%s_distort,%s_astr_interror2d,%s_astr_interror1d,%s_astr_referror2d,%s_astr_referror1d,%s_astr_chi2,%s_psphot_error' % (root,root,root,root,root,root,root,root)
             subprocess.call(['scamp', 'prepscamp.cat', '-c', config['scamp']['conf'], \
                 '-CHECKPLOT_DEV', 'PNG', \
                 '-CHECKPLOT_NAME', plotnames, \
                 '-VERBOSE_TYPE', verbose])
         else:
             subprocess.call(['scamp', 'prepscamp.cat', '-c', config['scamp']['conf'], \
                 '-CHECKPLOT_DEV', 'NULL', \
                 '-VERBOSE_TYPE', verbose])

         # Check astrometry offset
         with open('scamp.xml') as fd:
             doc = xmltodict.parse(fd.read())
         offset = doc['VOTABLE']['RESOURCE']['RESOURCE']['TABLE'][0]['DATA']['TABLEDATA']['TR']['TD'][34]
         daxis = offset.split(' ')
         daxis1 = float(daxis[0])
         daxis2 = float(daxis[1])
         daxis_mean = np.mean([daxis1,daxis2])

         pixelscale = doc['VOTABLE']['RESOURCE']['RESOURCE']['TABLE'][0]['DATA']['TABLEDATA']['TR']['TD'][18].split('  ')
         
         pixelscale = [float(pixelscale[0])/3600, float(pixelscale[1])/3600]
         # Update header of input fits file
         update_headers_scamp(ima,'prepscamp.head', pixelscale)

         # Delete temporary files
         clean_tmp_files(ima, soft='scamp')

         return daxis_mean
