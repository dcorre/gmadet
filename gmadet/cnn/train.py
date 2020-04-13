#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @file		train.py
# @brief	Train a CNN for detecting moving point-sources
# @date		01/11/2018
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

import sys, os, errno

import numpy as np
import matplotlib.pyplot as plt
import keras
import argparse
from keras.utils import multi_gpu_model
from gmadet.utils import getpath, rm_p, mkdir_p


def train(telescope, cubename, modelname):
    """Train CNN with simulated data"""

    gpus = -1 
    path_gmadet = getpath()
    inputdir = path_gmadet + '/cnn/datacube/' + telescope + '/'
    outdir = path_gmadet + '/cnn/CNN_training/' +  telescope + '/'
    mkdir_p(outdir)

    # Fraction of data used for the validation test
    fract = 0.1
    # define dropout percentage of each dropout
    dprob = np.array([0.3, 0.3, 0.3])
    # define padding
    padding = 'same'  # valid, same
    # number of epochs
    epochs = 10
    # name of the cube with data
    npz = cubename + '.npz'
    # outputname for the trained model
    model_name = outdir + "%s.h5" % modelname
     
    print("Loading " + inputdir +  npz + " ...", end='\r', flush=True)
    data = np.load(inputdir + npz)
    ima = data["cube"]
    lab = keras.utils.to_categorical(data["labels"])
    mag = data["mags"]
    errmag = data["errmags"]
    band = data["filters"]
    nclass = lab.shape[1]
    n = ima.shape[0]
    nt = int(n * fract)

    print("Shuffling data ...", end='\r', flush=True)
    randomize = np.arange(n)
    np.random.shuffle(randomize)
    ima = ima[randomize]
    lab = lab[randomize]
    mag = mag[randomize]
    errmag = errmag[randomize]
    band = band[randomize]

    print("Splitting dataset ...", end='\r', flush=True)
    imal = ima[nt:]
    labl = lab[nt:]
    magl = mag[nt:]
    errmagl = errmag[nt:]
    bandl = band[nt:]

    imat = ima[:nt]
    labt = lab[:nt]
    magt = mag[:nt]
    errmagt = errmag[:nt]
    bandt = band[:nt]

    model = keras.models.Sequential()

    model.add(keras.layers.Conv2D(64, (3, 3), activation='elu', padding=padding, input_shape=ima.shape[1:]))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.AveragePooling2D(pool_size=(2, 2)))
    #model.add(keras.layers.MaxPooling2D(pool_size=(2, 2)))
    model.add(keras.layers.Dropout(dprob[0]))
    model.add(keras.layers.Conv2D(128, (3, 3), activation='elu', padding=padding))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.MaxPooling2D(pool_size=(2, 2)))
    model.add(keras.layers.Dropout(dprob[1]))
    model.add(keras.layers.Conv2D(256, (3, 3), activation='elu', padding=padding))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.MaxPooling2D(pool_size=(2, 2)))
    model.add(keras.layers.Dropout(dprob[1]))
    model.add(keras.layers.Conv2D(256, (3, 3), activation='elu', padding=padding))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.MaxPooling2D(pool_size=(2, 2)))
    model.add(keras.layers.Dropout(dprob[2]))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(512, activation='elu'))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.Dropout(dprob[2]))
    model.add(keras.layers.Dense(32, activation='elu'))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.Dropout(dprob[2]))
    model.add(keras.layers.Dense(nclass, activation='softmax'))

    model.summary()
 
    if gpus > 0:
        parallel_model = multi_gpu_model(model, gpus=gpus)

        parallel_model.compile(loss='categorical_crossentropy',
                               optimizer=keras.optimizers.Adam(lr=0.001),
                               #optimizer=keras.optimizers.Nadam(),
                               metrics=['accuracy'])

        parallel_model.fit(imal, labl,
                           batch_size=1024,
                           epochs=epochs,
                           verbose=1,
                           validation_data=(imat, labt))

        score = parallel_model.evaluate(imat, labt, verbose=0)
        # save does not work on multi_gpu_model
        #parallel_model.save(model_name)
        labp = parallel_model.predict(imat)
   
    else: 

        model.compile(loss='categorical_crossentropy',
                      optimizer=keras.optimizers.Adam(lr=0.001),
                      #optimizer=keras.optimizers.Nadam(),
                      metrics=['accuracy'])
        #log = keras.callbacks.ModelCheckpoint('callbacks.h5', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto', period=1)
        #log = keras.callbacks(TensorBoard(log_dir='./logs', histogram_freq=5, batch_size=1024, write_graph=True, write_grads=False, write_images=False, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None, embeddings_data=None, update_freq='epoch'))
        model.fit(imal, labl,
                  batch_size=1024,
                  epochs=epochs,
                  verbose=1,
                  validation_data=(imat, labt))
        score = model.evaluate(imat, labt, verbose=0)
        labp = model.predict(imat)

    model.save(model_name)

    trange = np.arange(0.5,1.0,0.0001)
    fig, ax = plt.subplots()
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    mag_min = np.min(magt[magt != 99])
    mag_max = np.max(magt[magt != 99])
    for maglim in np.linspace(mag_min, mag_max, 6):
        labpm = labp[magt < maglim]
        labtm = labt[magt < maglim]
        labpf = labpm[labtm[:,1] <= 0.5]
        labpt = labpm[labtm[:,1] > 0.5]
        tpr = [np.mean(labpt[:,1] > t) for t in trange]
        fpr = [np.mean(labpf[:,1] > t) for t in trange]
        plt.plot(fpr,tpr, label = "mag < %.1f" %maglim)
    legend = ax.legend(loc='lower right')
    plt.savefig('ROC_mag.png')
    #plt.show()

    # ROC with dmag
    fig, ax = plt.subplots()
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    errmag_min = np.min(abs(errmagt[errmagt != 0]))
    errmag_max = np.max(abs(errmagt[errmagt != 0]))
    for errmaglim in np.linspace(errmag_min, errmag_max, 6):
        labpm = labp[errmagt < errmaglim]
        labtm = labt[errmagt < errmaglim]
        labpf = labpm[labtm[:,1] <= 0.5]
        labpt = labpm[labtm[:,1] > 0.5]
        tpr = [np.mean(labpt[:,1] > t) for t in trange]
        fpr = [np.mean(labpf[:,1] > t) for t in trange]
        plt.plot(fpr,tpr, label = "errmag < %.1f" %errmaglim)
    legend = ax.legend(loc='lower right')
    plt.savefig('ROC_errmag.png')
    #plt.show()

    fig, ax = plt.subplots()
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    for band in ['g', 'r', 'i', 'z']:
        labpm = labp[bandt == band]
        labtm = labt[bandt == band]
        labpf = labpm[labtm[:,1] <= 0.5]
        labpt = labpm[labtm[:,1] > 0.5]
        tpr = [np.mean(labpt[:,1] > t) for t in trange]
        fpr = [np.mean(labpf[:,1] > t) for t in trange]
        plt.plot(fpr,tpr, label = "%s" % band)
    legend = ax.legend(loc='lower right')
    plt.savefig('ROC_band.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Formatting data before stacking.')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    parser.add_argument('--cubename',
                        dest='cubename',
                        required=True,
                        type=str,
                        help='Name of the datacube')

    parser.add_argument('--modelname',
                        dest='modelname',
                        required=True,
                        type=str,
                        help='Name of the trained model')

    args = parser.parse_args()

    train(args.telescope, args.cubename, args.modelname)

