#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @file		infer.py
# @brief	Apply a previously trained CNN for detecting moving sources
# @date		29/10/2018
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# 	This file part of:	P9 search scripts
#
# 	Copyright:		(C) 2018 IAP/CNRS/SorbonneU
#
# 	Author:			Emmanuel Bertin (IAP)
#
# 	License:		GNU General Public License
#
# 	Bertinoscopic is free software: you can redistribute it and/or modify
# 	it under the terms of the GNU General Public License as published by
# 	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
# 	Bertinoscopic is distributed in the hope that it will be useful,
# 	but WITHOUT ANY WARRANTY; without even the implied warranty of
# 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# 	GNU General Public License for more details.
# 	You should have received a copy of the GNU General Public License
# 	along with Bertinoscopic. If not, see <http://www.gnu.org/licenses/>.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#    Original scrip modified by: David Corre (IJCLab/CNRS)

import sys

import errno
import glob
import os
import numpy as np
from astropy.io import fits
import keras
import argparse
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table

# from vis.visualization import visualize_cam
from vis.utils import utils
from keras import activations
from gmadet.utils import getpath, rm_p, mkdir_p
import matplotlib.pyplot as plt


def infer(path_cutouts, path_model, probratio):
    """Apply a previously trained CNN"""
    
    model_name = path_model 
    probratio = 1.0 / float(probratio)
    sstart = int(0)

    edge_shift = 64

    # Get all the images
    filenames = sorted(glob.glob(path_cutouts + "/*.fits"))
    print("Loading model " + model_name + " ...", end="\r", flush=True)

    model = keras.models.load_model(model_name)
    model.summary()

    # Load same model in another variable to be used and modified for
    # cam visualisation.
    # model_cam = keras.models.load_model(model_name)
    # Search for the layer index corresponding to the 2 last convolutional
    # layer
    # layer_idx = utils.find_layer_idx(model_cam, 'conv2d_4')
    # layer_idx2 = utils.find_layer_idx(model_cam, 'conv2d_3')
    # model_cam.layers[layer_idx].activation = activations.softmax

    # Type of gradient modifier to use
    grad = "small_values"
    # Type of backpropagation to use
    backprop = None
    # This is the output node we want to maximize.
    filter_idx = None

    cube = []
    for ima in filenames:
        hdus = fits.open(ima, memmap=False)
        head = hdus[0].header
        data = hdus[0].data
        # ima1 = ima1[edge_shift:-edge_shift,edge_shift:-edge_shift]
        cube.append(data)
        hdus.close()

    # Convert lists to B.I.P. NumPy arrays
    cube = np.asarray(cube, dtype=np.float32)
    print(cube.shape)
    if cube.ndim < 4:
        cube = np.reshape(
            cube, [
                cube.shape[0], cube.shape[1], cube.shape[2], 1])
    else:
        cube = np.moveaxis(cube, 1, -1)

    p = model.predict(cube)
    p2 = p[:, 1]
    # print (p)
    # print (p2)

    # label[j] = p / (p + (1.0 - p) * probratio)

    RA_list = []
    Dec_list = []
    original_file = []
    Xpos_list = []
    Ypos_list = []
    Cand_ID = []
    mag = []
    magerr = []
    FWHM = []
    FWHM_PSF = []

    for i in range(len(p)):
        hdus = fits.open(filenames[i], memmap=False)
        head = hdus[0].header
        RA_list.append(head["RA"])
        Dec_list.append(head["DEC"])
        original_file.append(head["FILE"])
        Xpos_list.append(head["XPOS"])
        Ypos_list.append(head["YPOS"])
        Cand_ID.append(head["CANDID"])
        mag.append(head["MAG"])
        magerr.append(head["MAGERR"])
        FWHM.append(head["FWHM"])
        FWHM_PSF.append(head["FWHMPSF"])

        # Format image cube to run the visualisation_cam
        # img = np.expand_dims(hdus[0].data, axis=0)
        # img = np.moveaxis(img, 1,-1)
        # im_res = visualize_cam(model_cam, layer_idx=layer_idx,
        #                        filter_indices=filter_idx,
        #                        seed_input=img,
        #                        penultimate_layer_idx=layer_idx2,
        #                        backprop_modifier=backprop,
        #                        grad_modifier=grad)

        # plt.imshow(im_res)
        # plt.show()
        # print (filenames[i], head['RA'], head['DEC'], p[i])
    table = Table(
        [
            filenames,
            RA_list,
            Dec_list,
            original_file,
            Xpos_list,
            Ypos_list,
            mag,
            magerr,
            FWHM,
            FWHM_PSF,
            Cand_ID,
            p[:, 0],
            p[:, 1],
        ],
        names=[
            "filename",
            "RA",
            "Dec",
            "originalFile",
            "Xpos",
            "Ypos",
            "mag",
            "magerr",
            "FWHM",
            "FWHMPSF",
            "candID",
            "label0",
            "label1",
        ],
    )
    # table.show_in_browser()
    table.write(path_cutouts + "/infer_results.dat",
                format="ascii.commented_header",
                overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Apply a previously trained CNN.")

    parser.add_argument(
        "--path_cutouts",
        dest="path_cutouts",
        required=True,
        type=str,
        help="Path to cutouts."
    )

    parser.add_argument(
        "--probratio",
        dest="probratio",
        required=True,
        type=float,
        help="Proba ratio, not used at the moment.",
    )

    parser.add_argument(
        "--path_model",
        dest="path_model",
        required=True,
        type=str,
        help="Path to the trained model, including its name and extension.",
    )

    args = parser.parse_args()

    infer(args.path_cutouts, args.path_model, args.probratio)
