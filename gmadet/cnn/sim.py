#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:15:55 2019

@author: David Corre (IJCLab/CNRS)
"""

import errno
import glob
import os
import random
import shutil
import subprocess
import sys
import numpy as np
import cv2
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
import matplotlib.pyplot as plt
import argparse
from skimage.feature import register_translation
from gmadet.utils import rm_p, mkdir_p


def sim(datapath, filenames, Ntrans=50, size=48,
        magrange=[14, 22], gain=None, magzp=30,
        radec= None, magnitude=None):
    """Insert point sources in real images """

    filenames = np.atleast_1d(filenames)

    #simdir = os.path.join(datapath, "simulation")
    simdir = datapath
    mkdir_p(simdir)

    cutsize = np.array([size, size], dtype=np.int32)

    hcutsize = cutsize // 2

    #  List to store position of simulated transients
    trans_pix = []
    trans_wcs = []
    filelist = []
    maglist = []
    filterlist = []
    trans_count = []
    cpsf1 = np.zeros((cutsize[1], cutsize[0]))

    counter = 0
    # To limit the number of images in which
    # to simulate stars, psf are also included in filenames
    # filenames = filenames[:4]
    for filename in filenames:
        if "psf" not in filename and "weight" not in filename:
            name = os.path.basename(filename)
            # print("\x1b[2K", end='\r', flush=True),
            # print("Loading " + epoch1 + " image data ...", end='\r',
            # flush=True),
            hdusi1 = fits.open(filename, memmap=False)
            headi1 = hdusi1[0].header
            band = str(headi1['FILTER'])
            ima1 = hdusi1[0].data.astype(np.float32)
            hdusp1 = fits.open(filename.split('.')[0] + '_psf.fits',
                               memmap=False)
            headp1 = hdusp1[0].header
            step1 = headp1['PSF_SAMP']
            nb_psf_snaps = int(headp1['PSF_NB'])
            psfs1 = hdusp1[0].data.astype(np.float32)
            imsize = ima1.shape
            posfac = np.array([nb_psf_snaps, nb_psf_snaps]) / imsize
            w = wcs.WCS(headi1)

            # try to use GAIN from header
            if gain is None:
                try:
                    gain = headp1["GAIN"]
                except BaseException:
                    print("GAIN keyword not found in header, set to 1.0.")
                    gain = 1.0

            # Add the transients to image
            pos = np.zeros((Ntrans, 2), dtype=float)
            for j in range(Ntrans):
                # newfile = os.path.join(simdir, os.path.splitext(name)[
                #                       0] + "_" + str(counter) + ".fits")
                # Keep same name actually, if works then remove line above.
                newfile = filename
                filelist.append(os.path.abspath(newfile))

                filterlist.append(band)
                
                if radec is None:
                    pos[j] = np.random.random_sample(2) * (imsize - cutsize) + cutsize / 2.0
                    # store positions
                    ra, dec = w.wcs_pix2world(pos[j][1], pos[j][0], 1)
                    trans_pix.append(pos[j])
                    trans_wcs.append([ra, dec])
     
                    # get pixels indexes in the image
                    ipos = pos[j].astype(int)
                    # extract subimage centered on position of the object and store in cima1
                    # same for weight maps
                    iposrange = np.s_[ipos[1] - hcutsize[1] : ipos[1] + hcutsize[1], \
                                      ipos[0] - hcutsize[0] : ipos[0] + hcutsize[0]]

                else:
                    ra, dec = radec[j][0], radec[j][1]
                    pos[j] = w.all_world2pix(ra, dec, 0)
                    trans_pix.append(pos[j])
                    trans_wcs.append([ra, dec])
                    # get pixels indexes in the image
                    ipos = pos[j].astype(int)
                    # print(ipos)
                    # extract subimage centered on position of the object and store in cima1
                    # same for weight maps
                    iposrange = np.s_[ipos[1] - hcutsize[1] : ipos[1] + hcutsize[1], \
                                      ipos[0] - hcutsize[0] : ipos[0] + hcutsize[0]]
                    # print(iposrange)
                # Select the PSF corresponding to the image area, in case there
                # are more than one PSF estimated per axis
                if nb_psf_snaps == 1:
                    psf1 = psfs1
                else:
                    # get position with respect to number of PSF snapshots
                    ppos = (pos[j] * posfac).astype(int)
                    psf1 = psfs1[ppos[0], ppos[1]]
                # step1 is psf_samp parameter from psfex, used in mat1 and 2 to rescale psf to image
                # psf1 is from psfex, and has different sampling than image1.
                # New object put in the center of this matrix
                mat1 = np.array([[step1, 0.0, hcutsize[0] - psf1.shape[0]*step1 / 2.0], \
                                 [0.0, step1, hcutsize[1] - psf1.shape[0]*step1 / 2.0]])

                # transformation of the PSF, resampling. 
                cpsf1 = cv2.warpAffine(psf1, mat1, cpsf1.shape, flags=cv2.INTER_LANCZOS4)


                if magnitude is None:
                # define the object magnitude using this random number and predefined ranges
                    mag = np.random.uniform(low=magrange[0], high=magrange[1], size=(1,))
                # convert the magnitude in ADU using the zeropoint magnitude.
                # Note that the zeropoint magnitude is define as 30, so did not care of the exact value
                # simply needed to draw random magnitudes. We could estimate the proper one for our telescopes
                    maglist.append(mag[0])
                    # trans_count.append(np.exp(0.921034 * (magzp - mag[0])))
                else:
                    mag = np.array([magnitude[j]])
                    maglist.append(mag[0])
                    # trans_count.append(np.exp(0.921034 * (magzp - mag[0])))
                amp1 = np.exp(0.921034 * (magzp - mag))
                # Apply Poisson Noise to simulated object
                noisy_object1 = cpsf1 * amp1#np.random.poisson(cpsf1 * amp1)
                noisy_object1[noisy_object1<0] = 0
                noisy_object1 = np.random.poisson(noisy_object1)
                noisy_object1 = noisy_object1 / gain
                
                ima1[iposrange] += noisy_object1
                count = np.sum(noisy_object1) 
            
                trans_count.append(count)
                # print('amp1=', amp1)
                # print('count=',count)
                # print('ima1 = ', np.sum(ima1[iposrange]))
            hdusi1[0].data = ima1
            # Write new fits file
            hdusi1.writeto(newfile, overwrite=True)

            hdusi1.close()
            hdusp1.close()
            counter+=1
    xypos = np.array(trans_pix)
    wcspos = np.array(trans_wcs)
    idx = np.arange(len(xypos))

    table = Table(
        [
            idx,
            filelist,
            xypos[:, 1],
            xypos[:, 0],
            wcspos[:, 0],
            wcspos[:, 1],
            maglist,
            filterlist,
            trans_count
        ],
        names=[
            "idx",
            "filename",
            "Xpos",
            "Ypos",
            "RA",
            "Dec",
            "mag",
            "filter",
            "count_ADU"],
    )
    table.write(
        os.path.join(simdir, "simulated_objects.list"),
        format="ascii.commented_header",
        overwrite=True,
    )
    return table

