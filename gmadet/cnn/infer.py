#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @file		infer.py
# @brief	Apply a previously trained CNN for detecting moving sources
# @date		29/10/2018
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
#if len(sys.argv) < 5:
#  print('syntax: ' + sys.argv[0] + ' <CNN_model> <prob_ratio> <event_start> <stack_prefix_1> [<stack_prefix2> ...]')
#  sys.exit(0)

import errno, glob, os
import numpy as np
from astropy.io import fits
import keras
import argparse
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from vis.visualization import visualize_cam
from vis.utils import utils
from keras import activations
from gmadet.utils import getpath, rm_p, mkdir_p
import matplotlib.pyplot as plt

def infer(eventDir, telescope, modelname, probratio):
    """Apply a previously trained CNN"""

    path_gmadet = getpath()

    trainedModelDir = path_gmadet + '/cnn/CNN_training/' + telescope + '/'

    model_name = trainedModelDir + modelname + '.h5'
    probratio = 1.0 / float(probratio)
    sstart = int(0)

    edge_shift = 64

    # Get all the images
    filenames = sorted(glob.glob(eventDir + '/*.fits'))
    print("Loading model " + model_name + " ...", end='\r', flush=True)

    model = keras.models.load_model(model_name)
    model.summary()

    # Load same model in another variable to be used and modified for cam visualisation
    #model_cam = keras.models.load_model(model_name)
    # Search for the layer index corresponding to the 2 last convolutional layer
    #layer_idx = utils.find_layer_idx(model_cam, 'conv2d_4')
    #layer_idx2 = utils.find_layer_idx(model_cam, 'conv2d_3')
    #model_cam.layers[layer_idx].activation = activations.softmax

    # Type of gradient modifier to use
    grad = 'small_values'
    # Type of backpropagation to use
    backprop = None
    # This is the output node we want to maximize.
    filter_idx = None

    cube = []
    for ima in filenames:
        hdus = fits.open(ima, memmap=False)
        head = hdus[0].header
        data = hdus[0].data
        #ima1 = ima1[edge_shift:-edge_shift,edge_shift:-edge_shift]
        cube.append(data)
        hdus.close()

    # Convert lists to B.I.P. NumPy arrays
    cube = np.asarray(cube, dtype=np.float32)
    print (cube.shape)
    if cube.ndim < 4:
        cube = np.reshape(cube, [cube.shape[0], cube.shape[1], cube.shape[2], 1])
    else:
        cube = np.moveaxis(cube, 1,-1)

    p = model.predict(cube)
    p2 = p[:,1]
    #print (p)
    #print (p2)

    #label[j] = p / (p + (1.0 - p) * probratio)

    RA_list = []
    Dec_list = []

    for i in range(len(p)):
        hdus = fits.open(filenames[i], memmap=False)
        head = hdus[0].header
        RA_list.append(head['RA'])
        Dec_list.append(head['DEC'])

        # Format image cube to run the visualisation_cam
        #img = np.expand_dims(hdus[0].data, axis=0)
        #img = np.moveaxis(img, 1,-1)
        #im_res = visualize_cam(model_cam, layer_idx=layer_idx,filter_indices=filter_idx,
        #                       seed_input=img, penultimate_layer_idx=layer_idx2, \
        #                       backprop_modifier=backprop, grad_modifier=grad)

        #plt.imshow(im_res)
        #plt.show()
        #print (filenames[i], head['RA'], head['DEC'], p[i])
    table = Table([filenames, RA_list, Dec_list, p[:,0], p[:,1]], names=['filename', 'RA', 'Dec', 'label0', 'label1'])
    #table.show_in_browser()
    table.write(eventDir+'/infer_results.dat', format='ascii.commented_header', overwrite=True)


    """
    size = 64
    cutsize = np.array([size, size], dtype=np.int32)
    hcutsize = cutsize // 2
    step = hcutsize
    cubebsq = np.zeros((3, cutsize[1], cutsize[0]), dtype=np.float32)
    tot=0
    for prefix in prefixes:
        imalists = sorted(glob.glob(dir_data + prefix + '*_??.fits'))
        epochs = []
        for imalist in imalists:
            epochs += [os.path.splitext(imalist)[0]]
        e1 = 0
        # take only 2 epochs
        epochs = epochs[:5]
        print(epochs)
        for epoch1 in epochs:
            e1 += 1
            hdusi1 = fits.open(epoch1 + '.fits', memmap=False)
            headi1 = hdusi1[0].header
            ima1 = hdusi1[0].data.astype(np.float32)
            ima1 = ima1[edge_shift:-edge_shift,edge_shift:-edge_shift]
            imasize = ima1.shape
            labelsize = (imasize - cutsize) // step + 1
            label = np.zeros(labelsize)
            e2 = e1
            for epoch2 in epochs[e1:]:
                e2 += 1
                hdusi2 = fits.open(epoch2 + '.fits', memmap=False)
                headi2 = hdusi2[0].header
                ima2 = hdusi2[0].data.astype(np.float32)
                ima2 = ima2[edge_shift:-edge_shift,edge_shift:-edge_shift]
                header = headi1.copy()
                header['CD1_1'] *= step[0]
                header['CD1_2'] *= step[1]
                header['CD2_1'] *= step[0]
                header['CD2_2'] *= step[1]
                header['CRPIX1'] /= step[0]
                header['CRPIX2'] /= step[1]
               
                dtime = (headi2['MJD-OBS'] - headi1['MJD-OBS']) * 24.0
                header['EXPTIME'] = dtime
                #print (epoch1,epoch2,label.shape, label)
                hdu = fits.PrimaryHDU(label)
                hdu.header = header
                cubebip = np.zeros((labelsize[1],cutsize[0], cutsize[1],3), dtype=np.float32)
                j = 0
                s = sstart
                for ypos in range(hcutsize[0], imasize[0] - hcutsize[0] + 1, step[0]):
                    i = 0
                    print("Processing " + epoch1.split('/')[-1] + " / " + epoch2.split('/')[-1] + " (%.1f h)" %dtime + \
		          " at line %d: %d events ..." %(ypos,s-sstart), end='\r', flush=True),
                    for xpos in range(hcutsize[1], imasize[1] - hcutsize[1] + 1, step[1]):
                        iposrange = np.s_[ypos - hcutsize[0] : ypos + hcutsize[0], \
		                          xpos - hcutsize[1] : xpos + hcutsize[1]]
                        cubebsq[0] = ima1[iposrange]
                        cubebsq[1] = ima2[iposrange]
                        if cleandata:
                            cubebsq[0][cubebsq[0]<0]=0
                            cubebsq[1][cubebsq[1]<0]=0
                        if data_compression:
                            cubebip[i] = np.arcsinh(np.moveaxis(cubebsq, 0,-1) / 10.0)
                        else:
                            cubebip[i] = np.moveaxis(cubebsq, 0,-1)
                        i += 1

                    p = model.predict(cubebip[:,:,:,:2])[:,1]
                    label[j] = p / (p + (1.0 - p) * probratio)
                    #label[j] = p
                    i = 0
                    for xpos in range(hcutsize[1], imasize[1] - hcutsize[1] + 1, step[1]):
                        if label[j,i] > 0.5:
                            s += 1
                            iposrange = np.s_[ypos - hcutsize[0] : ypos + hcutsize[0], \
		                              xpos - hcutsize[1] : xpos + hcutsize[1]]
                            cubebsq[0] = ima1[iposrange]
                            cubebsq[1] = ima2[iposrange]
                            # Format image cube to run the visualisation_cam
                            img = np.expand_dims(cubebsq[:2], axis=0)
                            img = np.moveaxis(img, 1,-1)
                            im_res = visualize_cam(model_cam, layer_idx=layer_idx,filter_indices=filter_idx,
                                                    seed_input=img, penultimate_layer_idx=layer_idx2, \
                                                    backprop_modifier=backprop, grad_modifier=grad)

                            cubebsq[2] = im_res
                            eventout = eventdir + "event_" + prefix.split('/')[-1] + "_%1d" %e1 + "_%1d" %e2 + \
		                       "_%04d.fits" %s
                            w1 = WCS(headi1)
                            ra, dec = w1.all_pix2world(xpos + 1 + edge_shift, ypos + 1 + edge_shift, 0)
                            print (xpos + 1 + edge_shift, ypos + 1 + edge_shift, ra, dec)
                            hduc = fits.PrimaryHDU(cubebsq)
                            hdrc = hduc.header
                            hdrc['XPOS'] = xpos + 1 + edge_shift
                            hdrc['YPOS'] = ypos + 1 + edge_shift
                            hdrc['EQUINOX'] = headi1['EQUINOX'] 
                            hdrc['MJD-OBS'] = headi1['MJD-OBS']
                            hdrc['RADESYS'] = headi1['RADESYS']
                            hdrc['CTYPE1'] = headi1['CTYPE1']
                            hdrc['CUNIT1'] = headi1['CUNIT1']
                            hdrc['CRVAL1'] = float(ra)
                            hdrc['CRPIX1'] = int(size/2)
                            hdrc['CD1_1'] = headi1['CD1_1']
                            hdrc['CD1_2'] = headi1['CD1_2']
                            hdrc['CTYPE2'] = headi1['CTYPE2']
                            hdrc['CUNIT2'] = headi1['CUNIT2']
                            hdrc['CRVAL2'] = float(dec)
                            hdrc['CRPIX2'] = int(size/2)
                            hdrc['CD2_1'] = headi1['CD2_1']
                            hdrc['CD2_2'] = headi1['CD2_2']
                            hducl = fits.HDUList([hduc])
                            hducl.writeto(eventout, overwrite=True)
                        i += 1
                    j += 1
                hdu.writeto(prefix + '_%d_%d.event.fits' %(e1,e2), overwrite=True)
                hdusi2.close()
                print("\x1b[2K", end='\r', flush=True),
                print("Processed " + epoch1.split('/')[-1] + " / " + epoch2.split('/')[-1] + " (%.1f h, %d events)" %(dtime, s-sstart), flush=True)
                tot +=s-1
            hdusi1.close()
            print('nb events total: ',tot)
    """
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Apply a previously trained CNN.')

    parser.add_argument('--path',
                        dest='path',
                        required=True,
                        type=str,
                        help='Path to images')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--probratio',
                        dest='probratio',
                        required=True,
                        type=float,
                        help='Proba ratio???')

    parser.add_argument('--modelname',
                        dest='modelname',
                        required=True,
                        type=str,
                        help='Name of the trained model')

    args = parser.parse_args()

    infer(args.path, args.telescope, args.modelname, args.probratio)

