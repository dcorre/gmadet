#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Python module for cutting TAROT images and performing astrometry.net
# Martin Blazek, IAA, Granada, Spain, alf@iaa.es
# July 2019, v1.2
#
# Example:
# python tarot_astro.py tarot.fits 
#
# Image is cut into 4 pieces with prefixes C1_, C2_, C3 and C4_
# In each subimage the basic WCS keywords are erased
# Then astrometry.net is performed
# You need to install locally astrometry.net and to download following indeces files: index-4204-*.fits
# This version uses Heasoft ftcopy command

import subprocess, sys, os
from astropy.io import fits
import shutil


def cut_image(filename, Nb_cuts=(2,2)):
    
    path, filename2 = os.path.split(filename)
    if path:
        folder = path + '/'
    else:
        folder = ''


    header = fits.getheader(filename)
    Naxis1 = header['NAXIS1']
    Naxis2 = header['NAXIS2']
    Naxis11 = int(Naxis1/Nb_cuts[0])
    Naxis22 = int(Naxis2/Nb_cuts[1])

    index=0
    for i in range(Nb_cuts[0]):
        for j in range(Nb_cuts[1]):
            index+=1
            if i == 0 :
                x1 = 1
            else:
                x1 = Naxis11 * i + 1
            if j == 0 :
                y1 = 1
            else:
                y1 = Naxis22 * j + 1
            x2 = Naxis11 * (i+1)
            y2 = Naxis22 * (j+1)
            
            extension = filename + "[%d:%d,%d:%d]" % (x1,x2, y1,y2)
            filename_out = folder + "D%d_" % (index) + filename2

            shutil.copyfile(filename, filename_out)
            hdul = fits.open(filename_out)
            fulltable = hdul[0].data
            subtable = fulltable[x1-1:x2-1,y1-1:y2-1]
            hdul[0].data = subtable
            hdul.writeto(filename_out,overwrite=True)
            hdul.close()

    #subprocess.call(['ftcopy',extension1,fileout1,'clobber=yes'])
    #subprocess.call(['ftcopy',extension2,fileout2,'clobber=yes'])
    #subprocess.call(['ftcopy',extension3,fileout3,'clobber=yes'])
    #subprocess.call(['ftcopy',extension4,fileout4,'clobber=yes'])


def clean_astrometry_temp_files(filename):
    fileroot = filename.split('.fits')
    os.remove(fileroot[0] + "-indx.xyls")
    os.remove(fileroot[0] + ".axy")
    os.remove(fileroot[0] + ".corr")
    os.remove(fileroot[0] + ".match")
    os.remove(fileroot[0] + ".rdls")
    os.remove(fileroot[0] + ".solved")
    os.remove(fileroot[0] + ".wcs")
    os.remove(fileroot[0] + ".fits")
    os.rename(fileroot[0] + ".new",fileroot[0] + ".fits") 

def erase_astrometry_header(filename):
    hdul = fits.open(filename)
    hdr = hdul[0].header
    del hdr['CRVAL1']
    del hdr['CRVAL2']
    del hdr['CRPIX1']
    del hdr['CRPIX2']
    hdul[0].header = hdr
    hdul.writeto(filename,overwrite=True)
    hdul.close()
    
    #subprocess.call(['fthedit', filename, 'CRVAL1', 'delete'])
    #subprocess.call(['fthedit', filename, 'CRVAL2', 'delete'])
    #subprocess.call(['fthedit', filename, 'CRPIX1', 'delete'])
    #subprocess.call(['fthedit', filename, 'CRPIX2', 'delete'])
    

def perform_astrometry(filename, radius=4, scaleLow=3.5, scaleHigh=4, scaleUnits='arcsecperpix'):
    print (filename)
    header = fits.getheader(filename)
    ra = str(header['CRVAL1'])
    dec = str(header['CRVAL2'])

    erase_astrometry_header(filename)

    subprocess.call(['solve-field', filename, '--ra', ra, '--dec', dec, '--radius', str(radius), '--scale-units', str(scaleUnits), '--scale-low', str(scaleLow), '--scale-high', str(scaleHigh), '--no-plots','--overwrite'])
    clean_astrometry_temp_files(filename)
    

if __name__ == "__main__":
    Nb_cuts = (4,4) 
    filename_argument = sys.argv[1]
    cut_image(filename_argument, Nb_cuts = Nb_cuts)

    index = 0
    for i in range(Nb_cuts[0]):
        for j in range(Nb_cuts[1]):
                index+=1
                filein = "D%d_" % (index) + filename_argument
    
                perform_astrometry(filein)
    

