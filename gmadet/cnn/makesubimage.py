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
import copy

def getCandPos(path, pattern='.alldetections'):
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
        # Select only detection in the substracted images and not
        # close to the edge.
        mask = (data['FlagSub'] == 'Y') & (data['edge'] == 'N')
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

def crossmatch_detections(path, candidates_list, radius=1):
    """ Crossmatch detections with simulated events positions """

    # Load simulated events summary file
    sim_list = ascii.read(path + '/simulated_objects.list')
    newfilename = []
    for sim in sim_list:
        filename = sim['filename'].split('/')[-1]
        newfilename.append(filename)
    sim_list['filename2'] = newfilename

    # Check whether the simulated events have been detected
    Ndet = []
    candID_list = []
    closest_candID_list = []

    Nsim = len(sim_list)
    for i, event in enumerate(sim_list):
        print ("processing simulated event  %d/%d ..." % (i, Nsim), end='\r', flush=True)
        event_coords = [event['RA'], event['Dec']]
        # Check whether it is a simulated object
        mask1 = candidates_list['OriginalIma'] == event['filename2']
        candidates= copy.deepcopy(candidates_list[mask1])
        # Compute the sepaation with detections and sort by ascending values
        offset = (candidates['RA'] - event['RA'])**2 + (candidates['Dec'] - event['Dec'])**2
        offset = np.sqrt(offset)
        candidates['offset'] = offset
        candidates.sort('offset')
        mask2 = candidates['offset'] < radius / 3600
        Nmatch = len(candidates[mask2])
        Ndet.append(Nmatch)
        if len(candidates[mask2]) > 0:
            text = ''
            for j in list(candidates['ID'][mask2]):
                text += '%s' % j
                if j < len(candidates[mask2]):
                    text += ','
            candID_list.append(text)
            closest_candID_list.append(candidates['ID'][mask2][0])
        else:
            candID_list.append(None)
            closest_candID_list.append(None)

    # Add a column 'detected' to the simulated events list table
    sim_list['Nmatches'] = Ndet
    sim_list['closest_candID'] = closest_candID_list
    sim_list['all_candIDs'] = candID_list
    sim_list.write(path+'/crossmatch.dat', format='ascii.commented_header', overwrite = True)
    return sim_list


def subimage(path, telescope, training, size=32, radius=1):
    """ Extract a sub-image centered on the candidate position """

    path_gmadet = getpath()

    # size of the extracted image
    cutsize = (size, size)

    print ('Combine the detections from all simulated images.')
    candidates_list = getCandPos(path)


    if training:
        print ('Crossmatch simulated events with detections')
        sim_list = crossmatch_detections(path, candidates_list, radius=radius)

        resdir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/'
        mkdir_p(resdir)
        truedir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/true/'
        mkdir_p(truedir)
        falsedir = path_gmadet + '/cnn/data/sim/' + telescope + '/candidates/false/'
        mkdir_p(falsedir)

    else:
        resdir = path + '/candidates/'
        mkdir_p(resdir)

    for i, cand in enumerate(candidates_list):
        print ("processing candidate %d/%d ..." % (i, len(candidates_list)), end='\r', flush=True)

        OT_coords = [cand['RA'], cand['Dec']]
        
        if training:
            # Check if corresponds to a simulated event
            mask = sim_list['closest_candID'] == cand['ID']
            # Check whether it is a simulated object
            #mask1 = sim_list['filename2'] == cand['OriginalIma']
            #mask2 = (sim_list['RA'] - cand['RA'])**2 + (sim_list['Dec'] - cand['Dec'])**2 < (radius/3600)**2
            #mask = np.bitwise_and(mask1,mask2)
            
            if len(sim_list[mask]) == 1:
                outdir = truedir
            else:
                outdir = resdir
            inputname = path_gmadet + '/cnn/data/sim/%s/images/gmadet_results/substraction/' % telescope + cand['filename'].split('/')[-1]
        else:
            outdir = resdir
            inputname = cand['filename']

        outname = outdir + 'candidate_%d.fits' % i
        # If the inputname can not be found for some reasons
        # Make sure the code is not crashing
        if True:
            # Get initial image size
            hdr_input = fits.getheader(inputname)
            Naxis1 = hdr_input['NAXIS1']
            Naxis2 = hdr_input['NAXIS2']

            # Extract small image
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
            # Whether it is close to the edge of the image
            # If yes the image will not be size x size in pixels
            if (cand['Xpos'] > Naxis1 - size) or (cand['Xpos'] < size) or (cand['Ypos'] > Naxis2 - size) or (cand['Ypos'] < size):
                hdr['edge'] = 'True'
                print ('Edge ', outname)
            else:
                hdr['edge'] = 'False'
            hdus.writeto(outname, overwrite=True)

        else:
            print ('Could not extract candidate in %s' % inputname)

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

    parser.add_argument('--size',
                        dest='size',
                        required=False,
                        default=32,
                        type=int,
                        help='Size in pixels of the extracted images. Default: 32')

    parser.add_argument('--radius',
                        dest='radius',
                        required=False,
                        default=2,
                        type=float,
                        help='Radius for crossmatching detected sources with simulated events. Default: 2 arcseconds')

    args = parser.parse_args()

    subimage(args.path, args.telescope, args.training, size=args.size, radius=args.radius)
