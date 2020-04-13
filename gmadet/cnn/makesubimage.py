#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: David Corre (IJCLab/CNRS)
"""

import errno, glob, math, os, shutil, subprocess, sys
import argparse
from gmadet.utils import getpath, make_sub_image, mv_p, rm_p, mkdir_p
from astropy.io import ascii, fits
from astropy.table import Table 
import numpy as np

def getCandPos(path, pattern='transient_candidates'):
    resfiles = glob.glob(path + '/**/gmadet_results/*%s*' % pattern, recursive=True)
    filelist = []
    origfilelist = []
    Xpos_list = []
    Ypos_list = []
    RA_list = []
    Dec_list = []
    mag_list = []
    magerr_list = []
    band_list = []
    for resfile in resfiles:
        data = ascii.read(resfile)
        mask = data['FlagSub'] == 'Y'
        RA_list.extend(data['_RAJ2000'][mask])
        Dec_list.extend(data['_DEJ2000'][mask])
        Xpos_list.extend(data['Xpos'][mask])
        Ypos_list.extend(data['Ypos'][mask])
        mag_list.extend(data['mag_calib'][mask])
        magerr_list.extend(data['mag_calib_err'][mask])
        band_list.extend(data['filter_cat'][mask])
        filelist.extend(data['filenames'][mask])
        _, origfile = os.path.split(data['OriginalIma'][mask][0])
        origfilelist.extend([origfile]*len(data[mask]))
    cand_id = np.arange(len(RA_list)) 
    table = Table([cand_id, RA_list, Dec_list, Xpos_list, Ypos_list, band_list, mag_list, magerr_list, filelist, origfilelist],
                   names=['ID', 'RA', 'Dec', 'Xpos', 'Ypos', 'Band', 'mag', 'magerr', 'filename', 'OriginalIma'])
    table.write(path+'/candidates_list.dat', format='ascii.commented_header', overwrite = True)

    return table


def subimage(path, telescope, training, size=64, radius=2):
    """ Extract a sub-image centered on the candidate position """

    path_gmadet = getpath()

    candidates_list = getCandPos(path)
    if training:
        sim_list = ascii.read(path + 'simulated_objects.list')
        newfilename = []
        for sim in sim_list:
            newfilename.append(sim['filename'].split('/')[-1])
        sim_list['filename2'] = newfilename

        resdir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/'
        mkdir_p(resdir)
        truedir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/true/'
        mkdir_p(truedir)
        falsedir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/false/'
        mkdir_p(falsedir)

    else:
        resdir = path + '/candidates/'
        mkdir_p(resdir)

    for i, cand in enumerate(candidates_list[:1000]):
        print ("processing candidate %d/%d ..." % (i, len(candidates_list)), end='\r', flush=True)

        OT_coords = [cand['RA'], cand['Dec']]
        
        if training:
            # Check whether it is a simulated object
            mask1 = sim_list['filename2'] == cand['OriginalIma']
            mask2 = (sim_list['RA'] - cand['RA'])**2 + (sim_list['Dec'] - cand['Dec'])**2 < (radius/3600)**2
            mask = np.bitwise_and(mask1,mask2)
            if len(sim_list[mask]) == 1:
                outdir = truedir
            else:
                outdir = resdir
            inputname = path_gmadet + '/cnn/data/sim/%s/images/gmadet_results/substraction/' % telescope + cand['filename'].split('/')[-1]
        else:
            outdir = resdir
            inputname = cand['filename']
        
        outname = outdir + 'candidate_%d.fits' % i
        make_sub_image(inputname, OT_coords, coords_type='world',
                       output_name=outname, size=[size,size],
                       fmt='fits', addheader=False)
        # add information to header
        hdus = fits.open(outname, memmap=False)
        hdr = hdus[0].header
        hdr['MAG'] = cand['mag']
        hdr['MAGERR'] = cand['magerr']
        hdr['FILTER'] = cand['Band']
        hdr['RA'] = cand['RA']
        hdr['Dec'] = cand['Dec']
        hdr['Xpos'] = cand['Xpos']
        hdr['Ypos'] = cand['Ypos']
        hdr['FILE'] = cand['filename']
        hdr['CANDID'] = cand['ID']
        hdus.writeto(outname, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Compute sub-images centered on candidates position.')

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

    parser.add_argument('--training',
                        dest='training',
                        action='store_true',
                        help='If set, training mode')

    args = parser.parse_args()

    subimage(args.path, args.telescope, args.training)
