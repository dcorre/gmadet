import numpy as np
import matplotlib.pyplot as plt
import keras
import argparse
from gmadet.utils import getpath

def hist(modelname, telescope):

    path_gmadet = getpath()

    trainedModelDir = path_gmadet + '/cnn/CNN_training/' + telescope + '/'
    model_name = trainedModelDir + modelname + '.h5'

    print("Loading " + modelname + " ...", end='\r', flush=True)
    data = np.load(modelname)
    ima = data["cube"]
    lab = keras.utils.to_categorical(data["labels"])
    mag = data["mags"]
    errmag = data["errmags"]
    band = data["filter"]
    errmag_min = np.min(errmag)
    errmga_max = np.max(errmag)
    mag_min = np.min(mag])
    mag_max = np.max(mag)

    plt.figure()
    hist_mag = plt.hist(mag,range = (mag_min, mag_max), bins = 20)
    plt.savefig('hist_mag.png')
    plt.figure()
    hist_dmag = plt.hist(dmag,range = (errmag_min,errmag_max), bins = 20, edgecolor = 'red')
    plt.savefig('hist_errmag.png')
    plt.figure()
    hist2D = plt.hist2d(mag[mag < 18.5],errmag[mag < 18.5],bins = 100)
    plt.savefig('hist2d.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Plot the distributions of some data used for the training.')

    parser.add_argument('--modelname',
                        dest='modelname',
                        required=True,
                        type=str,
                        help='Name of the trained model')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    args = parser.parse_args()

    hist(args.modelname, args.telescope)
