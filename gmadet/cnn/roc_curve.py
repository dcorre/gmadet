###Draw the roc curve
###Two way of doing it 
###1) macro/micro way (just smooth rpztation)
###2)with diff of magnitude(Bertin way)


import numpy as np
import matplotlib.pyplot as plt
import keras

from sklearn.metrics import roc_curve, auc
from scipy import interp
from itertools import cycle
import sys
import argparse

def draw(npz, model_name, name_curve):
    fract = 0.1
    dprob = 0.3

    

    print("Loading " + npz + " ...", end='\r', flush=True)
    data = np.load(npz)
    ima = data["cube"]
    lab = keras.utils.to_categorical(data["labels"])
    mag = data["mags"]
    dmag = data["dmags"]
    nclass = lab.shape[1]
    n = ima.shape[0]
    nt = int(n * fract)

    print("Shuffling data ...", end='\r', flush=True)
    randomize = np.arange(n)
    np.random.shuffle(randomize)
    ima = ima[randomize]
    lab = lab[randomize]
    mag = mag[randomize]
    dmag = dmag[randomize]
    print("Splitting dataset ...", end='\r', flush=True)
    imal = ima[nt:]
    labl = lab[nt:]
    magl = mag[nt:]
    dmagl = dmag[nt:]

    imat = ima[:nt]
    labt = lab[:nt]
    magt = mag[:nt]
    dmagt = dmag[:nt]
    model = keras.models.load_model(model_name)
    
    
    
    
    labp = model.predict(imat)

    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for i in range(nclass):
        fpr[i], tpr[i], _ = roc_curve(labt[:, i], labp[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    fpr["micro"], tpr["micro"], _ = roc_curve(labt.ravel(), labp.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(nclass)]))
    # Then interpolate all ROC curves at this points
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
    plt.show()

    plt.savefig(name_curve+'.png')
    trange = np.arange(0.5,1.0,0.0001)
    fig, ax = plt.subplots()
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    dmag_min = np.min(abs(dmagt[dmagt != 0]))
    dmag_max = np.max(abs(dmagt[dmagt != 0]))
    for dmaglim in np.linspace(dmag_min, dmag_max, 15):
        labpm = labp[dmagt < dmaglim]
        labtm = labt[dmagt < dmaglim]
        labpf = labpm[labtm[:,1] <= 0.5]
        labpt = labpm[labtm[:,1] > 0.5]
        tpr = [np.mean(labpt[:,1] > t) for t in trange]
        fpr = [np.mean(labpf[:,1] > t) for t in trange]
        plt.plot(fpr,tpr, label = "dmag < %.1f" %dmaglim)
    
    legend = ax.legend(loc='lower right')
    plt.show()
    plt.savefig('ROC_dmag.png')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Drawing roc_curve.')

    parser.add_argument('--path_npz',
                        dest='path_npz',
                        required=True,
                        type=str,
                        help='Path to get npz')

    parser.add_argument('--path_model',
                        dest='path_model',
                        required=True,
                        type=str,
                        help='path model')

    parser.add_argument('--destination_png',
                        dest='destination_png',
                        required=True,
                        type=str,
                        help='destination png')

    args = parser.parse_args()

    draw(args.path_npz,args.path_model,args.destination_png)










