#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: David Corre (IJCLab/CNRS)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from gmadet.utils import mkdir_p

class SortRes():

    def __init__(self, path_crossmatch, path_infer):
        self.data1 = ascii.read(path_crossmatch)
        self.data2 = ascii.read(path_infer)

    def filter_prob(self, data, prob, colname):
        mask = data[colname] >= prob
        return data[mask]

    def combine_match_cnn(self, prob, colname):
        idx_list = []
        label0 = []
        label1= []
        fwhm = []
        fwhm_psf = []
        mag = []
        magerr = []
        candID = []
        for row in self.data2:
            mask = (self.data1['closest_candID'] == row['cand_ID']) & (row[colname] > prob )
            if len(self.data1[mask]) == 1:
                idx_list.append(self.data1['idx'][mask][0])
                label0.append(row['label0'])
                label1.append(row['label1'])
                fwhm.append(row['FWHM'])
                fwhm_psf.append(row['FWHMPSF'])
                mag.append(row['mag'])
                magerr.append(row['magerr'])
                candID.append(row['cand_ID'])
        newtab = self.data1[idx_list]
        newtab['FWHM'] = fwhm
        newtab['FWHMPSF'] = fwhm_psf
        newtab['mag2'] = mag
        newtab['mag2err'] = magerr
        newtab['label0'] = label0
        newtab['label1'] = label1
        newtab['cand_ID'] = candID
        return newtab

    def filter_pos(self, data, RA, Dec, radius):
        offset = (data['RA'] - RA)**2 + (data['Dec'] - Dec)**2
        offset = np.sqrt(offset)
        data['offset'] = offset
        mask = data['offset'] < radius / 3600

        return data[mask]

    def hist(self, data, col, bins=20):
        plt.hist(data[col], bins,  density=False, facecolor='C0', alpha=0.75, log=True)


def makestats(path_plots, path_crossmatch, path_infer, maglim,
              CNNproblim, FWHM_ratio_lower, FWHM_ratio_upper):
    """ Create some statistics on the CNN training."""

    path_plots = os.path.join(path_plots, 'CheckInfer')
    mkdir_p(path_plots)
    res = SortRes(path_crossmatch, path_infer)

    # label0: probability that event is false
    # label1: probability that event is true
    # Combine information from the crossmatch perfromed on the simulated data
    # and the information from the CNN training.
    simID_list=res.combine_match_cnn(prob=0., colname='label1')

    # Create plots

    mask1 = simID_list['mag'] > maglim[2]
    mask2 = (simID_list['mag'] < maglim[2]) &   (simID_list['mag'] > maglim[1])
    mask3 = (simID_list['mag'] < maglim[1]) &   (simID_list['mag'] > maglim[0])

    plt.figure()
    bins = np.linspace(0,1,40)
    _,_,_ = plt.hist(simID_list['label1'][mask3],
                        bins=bins,
                        label='sim: %.1f<mag<%.1f '% (maglim[0], maglim[1]),
                        color='C0',
                        alpha=0.3)

    _,_,_ = plt.hist(simID_list['label1'][mask2],
                     bins=bins,
                     label='sim: %.1f<mag<%.1f' % (maglim[1], maglim[2]),
                     color='C1',
                     alpha=0.5)
    _,_,_ = plt.hist(simID_list['label1'][mask1],
                     bins=bins,
                     label='sim: mag>%.1f' % maglim[2],
                     color='C2',
                     alpha=0.75)

    plt.legend(loc='best')
    plt.xlabel('CNN prob')
    plt.title('Distribution of CNN prob for simulated events')
    plt.savefig(os.path.join(path_plots, 'distrib_CNN_prob_sim.png'))


    plt.figure()
    _,_,_ = plt.hist(res.data2['label1'],
                        bins=bins,
                        log=True,
                        color='C0',
                        alpha=0.2,
                        label='All: sim + real')

    _,_,_ = plt.hist(simID_list['label1'][mask3],
                        bins=bins,
                        label='sim: %.1f<mag<%.1f '% (maglim[0], maglim[1]),
                        color='C1',
                        alpha=0.3)
    _,_,_ = plt.hist(simID_list['label1'][mask2],
                     bins=bins,
                     label='sim: %.1f<mag<%.1f' % (maglim[1], maglim[2]),
                     color='C2',
                     alpha=0.5)
    _,_,_ = plt.hist(simID_list['label1'][mask1],
                     bins=bins,
                     label='sim: mag>%.1f' % maglim[2],
                     color='C3',
                     alpha=0.75)
    plt.legend(loc='best')
    plt.xlabel('CNN prob')
    plt.title('Distribution of CNN prob for all events')
    plt.savefig(os.path.join(path_plots, 'distrib_CNN_prob_all.png'))

    plt.figure()
    mask_sim = np.isin(res.data2['cand_ID'], simID_list['cand_ID'])
    tot,bins2,_ = plt.hist(res.data2['label1'],
                     bins=bins,
                     log=True,
                     color='C0',
                     alpha=0.3,
                     label='All sources (simulated included)')
    real,_,_ = plt.hist(res.data2['label1'][~mask_sim],
                     bins=bins,
                     log=True,
                     color='C1',
                     alpha=0.5,
                     label='Real sources')
    mask= simID_list['mag'] < maglim[-1]
    sim,_,_ = plt.hist(simID_list['label1'][mask],
                     bins=bins,
                     label='Simulated sources with mag < %.2f' % maglim[-1],
                     color='C2',
                     alpha=0.3)
    plt.xlabel('CNN prob')
    plt.legend(loc='best')
    plt.savefig(os.path.join(path_plots, 'distrib_CNNprob_all_and_sim.png'))

    x=(bins2[1:]+bins2[:-1])/2

    plt.figure()
    plt.plot(x, sim/np.sum(sim), color='C0', label='Simulated sources')
    plt.plot(x, real/np.sum(real), color='C1', label='Real sources')
    plt.plot(x, tot/np.sum(tot), color='C2', label='All sources')
    plt.xlabel('CNN prob')
    plt.ylabel('%')
    plt.legend(loc='best')
    plt.title('Fraction of sources per CNN prob bins')
    plt.savefig(os.path.join(path_plots, 'fraction_per_bin.png'))

    colmask = 'label1'
    coldistrib = 'FWHM'
    mask1 = simID_list[colmask] > CNNproblim[2]
    mask2 = (simID_list[colmask] < CNNproblim[2]) &   (simID_list[colmask] > CNNproblim[1])
    mask3 = (simID_list[colmask] < CNNproblim[1]) &   (simID_list[colmask] > CNNproblim[0])

    plt.figure()

    bins = np.linspace(0,20,60)
    _,_,_ = plt.hist(simID_list[coldistrib][mask3]/simID_list['FWHMPSF'][mask3],
                        bins=bins,
                        label='%.2f<%s<%.2f '% (CNNproblim[0], colmask,CNNproblim[1]),
                        color='C0',
                        alpha=0.5,
                        log=True)
    _,_,_ = plt.hist(simID_list[coldistrib][mask2]/simID_list['FWHMPSF'][mask2],
                     bins=bins,
                     label='%.2f<%s<%.2f' % (CNNproblim[1],colmask,CNNproblim[2]),
                     color='C1',
                     alpha=0.5)
    _,_,_ = plt.hist(simID_list[coldistrib][mask1]/simID_list['FWHMPSF'][mask1],
                     bins=bins,
                     label='%s>%.2f' % (colmask,CNNproblim[2]),
                     color='C2',
                     alpha=0.5)
    plt.axvline(FWHM_ratio_lower, color='red', ls='--')
    plt.axvline(FWHM_ratio_upper, color='red', ls='--')
    plt.legend(loc='best')
    plt.xlabel('FWHM / FWHM_PSF')
    plt.title('Distribution of CNN prob for simulated events')
    plt.savefig(os.path.join(path_plots, 'distrib_FWHM_simulation.png'))

    bins = np.linspace(0,20,60)
    plt.figure()
    """
    _,_,_ = plt.hist(res.data2['FWHM']/res.data2['FWHMPSF'],
                          bins=bins,
                          log=True,
                          color='C0',
                          alpha=0.4,
                          label='All sources (simulated included)')
    """
    _,_,_ = plt.hist(res.data2['FWHM'][~mask_sim]/res.data2['FWHMPSF'][~mask_sim],
                          bins=bins,
                          log=True,
                          color='C0',
                          alpha=0.4,
                          label='Real sources')
    mask= simID_list['label1'] >0.0
    _,_,_ = plt.hist(simID_list['FWHM'][mask]/simID_list['FWHMPSF'][mask],
                     bins=bins,
                     label='Simulated sources',
                     color='C1',
                     alpha=0.5)

    plt.axvline(0.5, color='red', ls='--')
    plt.axvline(4.2, color='red', ls='--')
    plt.legend()
    plt.xlabel('FWHM / FWHM_PSF')
    plt.legend(loc='best')
    plt.title('Distribution of CNN prob for all events.')
    plt.savefig(os.path.join(path_plots, 'distrib_FWHM_real_vs_sim.png'))
