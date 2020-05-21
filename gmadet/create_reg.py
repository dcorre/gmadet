#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: David Corre
# email: corre@lal.in2p3.fr

import errno, glob, os, shutil, subprocess, sys
import numpy as np
from astropy.io import ascii, fits
from astropy.table import vstack, Table

path = 'test/test/gmadet_results/substraction/'
#path = 'test/test/gmadet_results/'
files = glob.glob(path + '*.magwcs', recursive=False)

regfile = open('sub.reg', 'w')
#regfile = open('notsub.reg', 'w')
for file in files:
    data = ascii.read(file)

    for row in data:

        regfile.write('icrs;point(%.7f, %.7f) ' % \
                      (row['RA'], row['DEC']) + \
                      '# color=red ' + \
                      'text={%.2f %.2f %.2f %.2f %d} ' % (row['psf_chi2'], row['model_chi2'], row['FWHM'], row['FWHMPSF'], row['flag_psf']) + \
                      'font="times 15 bold"\n')

regfile.close()

