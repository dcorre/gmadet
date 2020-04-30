###Draw the roc curve
###Two way of doing it 
###1) macro/micro way (just smooth rpztation)
###2)with diff of magnitude 


import numpy as np
import matplotlib.pyplot as plt
import keras

from sklearn.metrics import roc_curve, auc
from scipy import interp
from itertools import cycle
import sys
import argparse
from gmadet.utils import getpath, mkdir_p

def draw(cubename, modelname, telescope):

    path_gmadet = getpath()

    inputdir = path_gmadet + '/cnn/datacube/' + telescope + '/'
    npz = cubename + '.npz' 

    trainedModelDir = path_gmadet + '/cnn/CNN_training/' + telescope + '/'
    model_name = trainedModelDir + modelname + '.h5'

    outdir = trainedModelDir + 'plots/'
    mkdir_p(outdir)

    print("Loading " + inputdir +  npz + " ...", end='\r', flush=True)
    data = np.load(inputdir + npz)
    ima = data["cube"]
    lab = keras.utils.to_categorical(data["labels"])
    mag = data["mags"]
    errmag = data["errmags"]
    band = data["filters"]
    nclass = lab.shape[1]
    n = ima.shape[0]

    model = keras.models.load_model(model_name)
    
    labp = model.predict(ima)

    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for i in range(nclass):
        fpr[i], tpr[i], _ = roc_curve(lab[:, i], labp[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    fpr["micro"], tpr["micro"], _ = roc_curve(lab.ravel(), lab.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(nclass)]))
    # Then interpolate all ROC curves at these points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(nclass):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= nclass

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    lw = 2
    plt.figure(1)
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
    for i, color in zip(range(nclass), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Some extension of Receiver operating characteristic to multi-class')
    plt.legend(loc="lower right")
    #plt.show()

    plt.savefig(outdir+'ROC.png')
    fig, ax = plt.subplots()
    trange = np.arange(0.5,1.0,0.0001)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    errmag_min = np.min(errmag[errmag < 0.5])
    errmag_max = np.max(errmag[errmag < 0.5])
    print (errmag_min, errmag_max)
    print (np.linspace(errmag_min, errmag_max, 5))
    for errmaglim in np.linspace(errmag_min, errmag_max, 5):
        labpm = labp[errmag < errmaglim]
        labm = lab[errmag < errmaglim]
        labpf = labpm[labm[:,1] <= 0.5]
        labpt = labpm[labm[:,1] > 0.5]
        tpr = [np.mean(labpt[:,1] > t) for t in trange]
        fpr = [np.mean(labpf[:,1] > t) for t in trange]
        plt.plot(fpr,tpr, label = "errmag < %3f" % errmaglim)
    
    legend = ax.legend(loc='lower right')
    #plt.show()
    plt.savefig(outdir + 'ROC_errmag.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Drawing roc_curve.')

    parser.add_argument('--datacube',
                        dest='datacube',
                        required=True,
                        type=str,
                        help='Path to get npz datacube')

    parser.add_argument('--modelname',
                        dest='modelname',
                        required=True,
                        type=str,
                        help='name of the cnn trained model')

    parser.add_argument('--telescope',
                        dest='telescope',
                        required=True,
                        type=str,
                        help='Telescope identifier')

    args = parser.parse_args()

    draw(args.datacube, args.modelname, args.telescope)

