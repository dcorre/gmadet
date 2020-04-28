import numpy as np
import matplotlib.pyplot as plt
import keras
import argparse
from gmadet.utils import getpath

def hist(modelname, telescope):

    path_gmadet = getpath()

    cubeDir = path_gmadet + '/cnn/datacube/' + telescope + '/'
    model_name = cubeDir + modelname + '.npz'

    print("Loading " + model_name + " ...", end='\r', flush=True)
    data = np.load(model_name)
    ima = data["cube"]
    lab = keras.utils.to_categorical(data["labels"])
    mag = data["mags"]
    errmag = data["errmags"]
    band = data["filters"]
    errmag_min = np.min(errmag)
    errmga_max = np.max(errmag)
    mag_min = np.min(mag)
    mag_max = np.max(mag)

    plt.figure()
    hist_mag = plt.hist(mag,range = (mag_min, mag_max), bins = 20)
    plt.savefig(cubeDir+modelname+'_hist_mag.png')
    plt.figure()
    hist_dmag = plt.hist(dmag,range = (errmag_min,errmag_max), bins = 20, edgecolor = 'red')
    plt.savefig(cubeDir+modelname+'_hist_errmag.png')
    plt.figure()
    hist2D = plt.hist2d(mag[mag < 18.5],errmag[mag < 18.5],bins = 100)
    plt.savefig(cubeDir+modelname+'_hist2d.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Plot the distributions of some data used for the training.')

    parser.add_argument('--modelname',
                        dest='modelname',
                        required=True,
                        type=str,
                        help='Name of the datacube')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    args = parser.parse_args()

    hist(args.modelname, args.telescope)
