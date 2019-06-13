#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, subprocess, glob, math, shutil

# Create login.cl at execution of the script if flag set to true
if sys.argv[3] == True:
    proc = subprocess.Popen(['mkiraf'], stdin = subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    outs, errs = proc.communicate('y\nxterm')

from catalogues import *
from sources_det import psf

from pyraf import iraf
from pyraf.iraf import daophot

#filename = sys.argv[1]
#filename without .fits suffix


def get_photometry(filename,fwhmpsf):


    iraf.noao
    iraf.digiphot
    iraf.daophot

    iraf.unlearn('phot')
    iraf.unlearn('datapars')
    iraf.unlearn('photpars')
    iraf.unlearn('daopars')

    iraf.datapars.epadu = 1
    iraf.datapars.readnoi = 6
    iraf.datapars.datamin = "INDEF"
    iraf.datapars.datamax = "50000"

    # fwhm tarot 1.5   oaj 3.5
    iraf.datapars.fwhm = fwhmpsf
    THRESH = 3

    IMstatres = iraf.imstat(filename,Stdout=1)
    IMmean = IMstatres[1].split()[2]
    iraf.datapars.sigma = math.sqrt(float(IMmean) * float(iraf.datapars.epadu) + float(iraf.datapars.readnoi)**2) / float(iraf.datapars.epadu)

    iraf.daofind(filename,output="default",verify="no", threshold=5)
    iraf.datapars.datamax = "INDEF"
    iraf.daophot.phot(image=filename,coords="default", output="default", interactive="no", sigma="INDEF", airmass="AIRMASS", exposure="EXPOSURE", filter="FILTER", obstime="JD", calgorithm="gauss", verify="no", verbose="yes")


def psfex(filename):
    psf(filename) 
    psffilename = filename.split('.')[0] + '.psf.fits'


def select_good_stars(filename,limiting_mag_err):
    list_of_magfiles = glob.glob("./*.mag.*")
    magfile = list_of_magfiles[0]
    resmaggile = filename+".rmg"
    f1 = open(magfile, "r")
    f2 = open(resmaggile,"w")

    for kk in range(1,77):
        lajna = f1.readline()

    while lajna:
        
        lajna = f1.readline()
        xpos = lajna.split()[0]
        ypos = lajna.split()[1]
        lajna = f1.readline()
        lajna = f1.readline()
        lajna = f1.readline()
        mag = lajna.split()[4]
        merr = lajna.split()[5]
        errmessage = lajna.split()[7]

        
    
        if (errmessage == "NoError") or (mag != "INDEF") or (merr != "INDEF"):
            if (float(merr) < limiting_mag_err):
                print(xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                f2.write(xpos+" "+ypos+" "+mag+" "+merr+"\n")
            else:
                print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
        else:
            print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
        lajna = f1.readline()

    f1.close()
    f2.close()

def convert_xy_radec(filename):
    magfile = filename+".rmg"
    magfilewcs = filename+".rwc"
    iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")


def crosscheck_with_catalogues(filename,degrad):
    magfilewcs = filename+".rwc"
    transients = filename+".oc"
    print(magfilewcs)
    
    f3 = open(magfilewcs,"r")
    lajna = f3.readline()
    lajna = f3.readline()
    lajna = f3.readline()
    lajna = f3.readline()
    
    f4 = open(transients,"w")
    
    while lajna:
        print(lajna)
        ra = lajna.split()[0]
        dec = lajna.split()[1]
        try:
            wohoo = USNO_B1_query(ra,dec,degrad)
            print(wohoo)

        except:
            print("NEW TRANSIENT")
            f4.write(ra+" "+dec+"\n")
        
        
        
        
        #wohoo['RAJ2000']
        lajna=f3.readline()
    
    f3.close()
    f4.close()
    



psfex(sys.argv[1])

#get_photometry(sys.argv[1],sys.argv[2])
#select_good_stars(sys.argv[1],0.5)
#convert_xy_radec(sys.argv[1])
#crosscheck_with_catalogues(sys.argv[1],3*0.000907203)
#print('baf')





