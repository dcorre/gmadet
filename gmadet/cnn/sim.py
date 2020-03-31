# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:15:55 2019

@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, os, random, shutil, subprocess, sys
import numpy as np
import cv2
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import argparse
from skimage.feature import register_translation

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


def sim(path_stacks, telescope, useweight):
    """Insert point sources with varying mag in real images """

    dir = path_stacks + telescope + '/'
    
    cutdir = 'cutouts/' + telescope + '/'
    eventdir = 'events/' + telescope + '/'
    mkdir_p(cutdir)
    mkdir_p(eventdir)

    useweight = bool(useweight)


#parametres a modifier##################################
    type = 'variable' #'variable' # 'moving' 'variable'
    data_compression=False
    cleandata=False
    same_image = False
    
########################################################

    npair = 24000
    size = 64
    eventfrac = 0.2
    minshift = 15.0
    maxshift = 30.0
    magzp = 30.0
    if type == 'moving':
        magmax = 19
        dmagmax = 0.05
        magrange = 4
    elif type == 'variable':
        magmax = 20
        dmagmax = 1
        magrange = 4
        dmag_min = 0.2
    elif type == 'appears':
        magmax = 22
        dmagmax = 0.001
        magrange = 3 #you can modify it
    
    fwhm = 2.0
    cutsize = np.array([size, size], dtype=np.int32)

    hcutsize = cutsize // 2

    ngoodpix = int(0.5*(cutsize[0]*cutsize[1]))

    cube = np.zeros((2,cutsize[1], cutsize[0]))
    cpsf1 = np.zeros((cutsize[1], cutsize[0]))
    cpsf2 = np.zeros((cutsize[1], cutsize[0]))

    
    # Get all the prefixes corresponding to one field
    filenames = glob.glob(dir + '*.fits')
    
    prefixes = []
    for filename in filenames:
        if 'psf' not in filename:
            splitfilename = os.path.splitext(filename)[0].split('/')[-1].split('_')
            prefixes.append(splitfilename[0] + '_' + splitfilename[1])
    # Discard duplicates
    prefixes = np.unique(sorted(prefixes))

    # take only one event
    prefixes=[prefixes[1], prefixes[2]]
    print (prefixes)
    
    for prefix in prefixes:
        #print (prefix)
        imalists = sorted(glob.glob(dir + prefix + '_??.fits'))
        #print (imalists)
        epochs = []
        for imalist in imalists:
            epochs += [os.path.splitext(imalist)[0]]
        
        e1 = 0
        # take only 2 images
        epochs=epochs[:5]
        print (epochs)
        for epoch1 in epochs:
            e1 += 1
          
            print("\x1b[2K", end='\r', flush=True),
            print("Loading " + epoch1 + " image data ...", end='\r', flush=True),
            hdusi1 = fits.open(epoch1 + '.fits', memmap=False)
            headi1 = hdusi1[0].header
            ima1 = hdusi1[0].data.astype(np.float32)
            if useweight:
                hdusw1 = fits.open(epoch1 + '.weight.fits', memmap=False)
                weight1 = hdusw1[0].data.astype(np.float32)
            hdusp1 = fits.open(epoch1 + '.psf.fits', memmap=False)
            headp1 = hdusp1[0].header
            step1 = headp1['PSF_SAMP']
            psfs1 = hdusp1[0].data.astype(np.float32)
            if psfs1.shape[0] == 1 and psfs1.shape[1] == 1:
                continue

            imsize = ima1.shape
            posfac = np.array([10.0, 10.0]) / imsize
            e2 = e1
            if same_image == True:
                Epochs = [epochs[e1-1]]
            else:
                Epochs = epochs[e1:]
            print('Epochs : ',Epochs)
            for epoch2 in Epochs:
                
                e2 += 1
                print("\x1b[2K", end='\r', flush=True),
                print("Loading " + epoch2 + " image data ...", end='\r', flush=True),
                hdusi2 = fits.open(epoch2 + '.fits', memmap=False)
                headi2 = hdusi2[0].header
                ima2 = hdusi2[0].data.astype(np.float32)
                if useweight:
                    hdusw2 = fits.open(epoch2 + '.weight.fits', memmap=False)
                    weight2 = hdusw2[0].data.astype(np.float32)
                hdusp2 = fits.open(epoch2 + '.psf.fits', memmap=False)
                headp2 = hdusp2[0].header
                step2 = headp2['PSF_SAMP']
                psfs2 = hdusp2[0].data.astype(np.float32)
                if psfs2.shape[0] == 1 and psfs2.shape[1] == 1:
                    continue

                # Check previous events for current image pair in event repository
                eventlist = sorted(glob.glob(eventdir + "event_" + prefix + \
	                	"_%02d" %e1 + "_%02d" %e2 + "_????.fits"))
     
                # Set up vector of position coordinates
                npos = npair + len(eventlist)
                pos = np.zeros((npos, 2), dtype=float)

                for j in range(npair):
                    pos[j] = np.random.random_sample(2) * (imsize - cutsize) + cutsize / 2.0
                    
                # Get event positions previously computed
                j = npair
                for event in eventlist:
                    print("\x1b[2K", end='\r', flush=True),
                    print("Reading " + event + " ...", end='\r', flush=True),
                    hdus = fits.open(event, memmap=False)
                    head = hdus[0].header
                    # get positions from header
                    pos[j] = np.array([int(head['YPOS']) - 1, int(head['XPOS']) - 1]).astype(np.float)
                    # Randomize slightly positions from previous events
                    # Randomly draw new offset and check that it is falling in the image size.
                    # seems not to be used as dpos is overwritten later on
                    while True:
                        dpos = (0.5 - np.random.random_sample(2)) * cutsize / 4.0
                        tpos = pos[j] + dpos
                        if np.all(tpos >= 0.0) and np.all(tpos < imsize):
                            break
                    hdus.close()
                    j += 1
                print("\x1b[2K", end='\r', flush=True),
                # Copy or augment data at every recorded position

                s = 0
                for j in range(len(pos)):
                    s += 1
                    print("Processing " + epoch1.split('/')[-1] + "    " + epoch2.split('/')[-1] + \
	 	         " pair #%d/%d ..." %(s,npos), end='\r', flush=True),
                    cutout = cutdir + prefix.split('/')[-1] + "_%02d" %e1 + "_%02d" %e2 + "_%04d.fits" %s
                    # get pixels indexes in the image
                    ipos = pos[j].astype(int)
                    
                    # get position with psf oversampling??
                    ppos = (pos[j] * posfac).astype(int)

                    # extract subimage centered on position of the object and store in cima1 and 2
                    # same for weight maps
                    iposrange = np.s_[ipos[0] - hcutsize[0] : ipos[0] + hcutsize[0], \
	                              ipos[1] - hcutsize[1] : ipos[1] + hcutsize[1]]
                    cima1 = ima1[iposrange]
                    cima2 = ima2[iposrange]

                    if useweight:
                        cweight1 = weight1[iposrange]
                        cweight2 = weight2[iposrange]

                    psf1 = psfs1[ppos[0],ppos[1]]
                    psf2 = psfs2[ppos[0],ppos[1]]
                    # define object postion at epoch1. Avoid to place it at the edge
                    # using the typical fwhm. Can be placed at pixels > fwhm from the edge. but only for edges x=0, and y=0, so half of the edges, why not all of them?
                    # might use np.random.uniform(low=fwhm, high=cutsize[0] - fwhm, size=(2,)) instead
                    pos1 = np.random.random_sample(2) * (cutsize - fwhm) + fwhm
 
                    # randomly set label to true or false, 50/50
                    label = True if random.random() < 0.5 else False
                    # if label true, shift object from epoch1 to 2 by dpos
                    # can be only between -maxshift and +maxshift and must be located in image
                    if label:
                        # draw the difference in magnitude betwen epoch1 and epoch2
                        if type == 'variable':
                            dmag=0
                            while (dmag > -dmag_min) & (dmag < dmag_min):
                                dmag = random.gauss(0.0, dmagmax)
                        else:
                            dmag = random.gauss(0.0, dmagmax)

                        if type != 'moving':
                            dpos = [0, 0]

                        """while True:
                            dpos = 2.0 * (np.random.random_sample(2) - 0.5) * maxshift
                            if np.sum(dpos * dpos) >= minshift * minshift and \
                       	      np.sum(dpos * dpos) < maxshift * maxshift and \
	                      np.all((pos1 + dpos) >= 0.0) and \
                 	      np.all((pos1 + dpos) < (cutsize - fwhm)):
                                break
                        """
                        # empty = False means object has been placed in subimage
                        empty = False
                    # if label false, same position at epoch1 and 2

                    else:
                        dpos = [0.0, 0.0]
                        dmag = 0
                        """if useweight:
                            empty = False if np.count_nonzero(cweight1) > ngoodpix and \
 		                        np.count_nonzero(cweight2) > ngoodpix and \
                                         random.random() < 0.5 else True
                        else:
                            empty = False if random.random() < 0.5 else True"""
                        empty = True
                        #empty = False if random.random() < 0.5 else True

                    # no object placed in image if empty = True
                    if empty:
                        if cleandata:
                            cima1[cima1<0]=0
                            cima2[cima2<0]=0
                        if data_compression:
                            cube[0]=np.arcsinh((cima1)/10)
                            cube[1]=np.arcsinh((cima2)/10)
                            #cube[0] = cima_obj1/(2**16-1)
                            #cube[1] = cima_obj2/(2**16-1)
                            #cube[0] = cima_obj1
                            #cube[1] = cima_obj2
                        else:
                            cube[0] = cima1
                            cube[1] = cima2

                        pos1 = [0.0, 0.0]
                        pos2 = [0.0, 0.0]
                        mag = 99.0
                        #dmag = 0.0
                    else:

                        pos2 = pos1 + dpos
                        # step1 is psf_samp parameter from psfex, used in mat1 and 2 to rescale psf to image
                        # pos1 is the x position of object in image1, rescaled with psf oversampling in mat1 and 2
                        # psf1 is from psfex, and has different sampling than image1.
                        # Pb is that psf1 seems to be a float and so has no .shape attribute... Probably misunderstood something
                        mat1 = np.array([[step1, 0.0, pos1[0] - psf1.shape[0]*step1 / 2.0], \
 		                	[0.0, step1, pos1[1] - psf1.shape[1]*step1 / 2.0]])
                        mat2 = np.array([[step2, 0.0, pos2[0] - psf2.shape[0]*step2 / 2.0], \
 		                	[0.0, step2, pos2[1]  - psf2.shape[1]*step2 / 2.0]])
                        # transformation of the PSF, resampling + relocalisation. need to understand better
                        # the transformation matrix defined above
                        cpsf1 = cv2.warpAffine(psf1, mat1, cpsf1.shape, flags=cv2.INTER_LANCZOS4)
                        cpsf2 = cv2.warpAffine(psf2, mat2, cpsf2.shape, flags=cv2.INTER_LANCZOS4)
                        # draw random number between 0 and 1
                        r2 = np.random.random()
                        # define the object magnitude using this random number and predefined ranges
                        mag = magmax - r2 * r2 * magrange
                        
                        # convert the magnitude in ADU using the zeropoint magnitude.
                        # Note that the zeropoint magnitude is define as 30, so did not care of the exact value
                        # simply needed to draw random magnitudes. We could estimate the proper one for our telescopes
                        
                        if type == 'appears':
                            if label == True:
                                if random.random() < 0.5:
                                    amp1 = 0
                                    amp2 = np.exp(0.921034 * (magzp - mag))
                                else:
                                    amp1 = np.exp(0.921034 * (magzp - mag))
                                    amp2 = 0
                            else:
                                amp1 = np.exp(0.921034 * (magzp - mag))
                                amp2 = np.exp(0.921034 * (magzp - mag ))
                        else :
                            amp1 = np.exp(0.921034 * (magzp - mag))
                            amp2 = np.exp(0.921034 * (magzp - mag - dmag))
  
                        # Apply Poisson Noise to simulated object
                        noisy_object1 = cpsf1 * amp1#np.random.poisson(cpsf1 * amp1)
                        noisy_object2 = cpsf2 * amp2#np.random.poisson(cpsf2 * amp2)
                        noisy_object1[noisy_object1<0] = 0
                        noisy_object2[noisy_object2<0] = 0
                        noisy_object1 = np.random.poisson(noisy_object1)
                        noisy_object2 = np.random.poisson(noisy_object2)
                        noisy_object1 = noisy_object1 / 3.17
                        noisy_object2 = noisy_object2 / 3.17
                        # Why arcsinh? because it increases the dynamic range for low ADU counts?
                        if useweight:
                            if data_compression:
                                cube[0] = np.arcsinh((cima1 + noisy_object1 * (cweight1 > 0.01))/10 )
                                cube[1] = np.arcsinh((cima2 + noisy_object2 * (cweight2 > 0.01))/10 )
                            else:
                                cube[0] = cima1 + noisy_object1 * (cweight1 > 0.01)
                                cube[1] = cima2 + noisy_object2 * (cweight2 > 0.01)
                        else:
                            cima_obj1 = cima1 + noisy_object1 
                            cima_obj2 = cima2 + noisy_object2
                            if cleandata:
                                cima_obj1[cima_obj1<0]=0
                                cima_obj2[cima_obj2<0]=0
                            if data_compression:
                                cube[0]=np.arcsinh((cima_obj1)/10)
                                cube[1]=np.arcsinh((cima_obj2)/10)
                                #cube[0] = cima_obj1/(2**16-1) 
                                #cube[1] = cima_obj2/(2**16-1) 
                                #cube[0] = cima_obj1
                                #cube[1] = cima_obj2   
                            else:
                                cube[0] = cima_obj1
                                cube[1] = cima_obj2
                            #cube[0]=cima_obj1
                            #cube[1]=cima_obj2	
                    # mask negative value / background
                    #mask = ((cube[0]<3) | (cube[1]<3))
                    #cube[0][mask]=0
                    #cube[1][mask]=0
                    # write the cubes, for epoch1 and epoch2 in a single fits file.
                    # labels, mag, dmag and positions are stored as tables in the header
                    hdu = fits.PrimaryHDU(cube)
                    hdr = hdu.header
                    hdr['EMPTY'] = empty
                    hdr['EVENT'] = label
                    hdr['SIMMAG'] = mag
                    hdr['SIMDMAG'] = dmag
                    hdr['XPOS1'] = pos1[1]
                    hdr['YPOS1'] = pos1[0]
                    hdr['XPOS2'] = pos2[1]
                    hdr['YPOS2'] = pos2[0]
                    hdul = fits.HDUList([hdu])
                    hdul.writeto(cutout, overwrite=True)
                    #if label == True and mag < 17:
                    #    print ("\n\n",cutout)
                    #    hdusi2.close()
                    #    hdusp2.close()
                    #    stop
                hdusi2.close()
                if useweight: hdusw2.close()
                hdusp2.close()
                print("\x1b[2K", end='\r', flush=True),
                print("Processed " + epoch1 + "/" + epoch2 + \
	         	" with %d positions" %npos, end='\r', flush=True)

            hdusi1.close()
            if useweight: hdusw1.close()
            hdusp1.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Stacking images.')

    parser.add_argument('--path_stacks',
                        dest='path_stacks',
                        required=True,
                        type=str,
                        help='Path to stack images')

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

    sim(args.path_stacks,args.telescope, useweight=args.useweight)
