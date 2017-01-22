#!/usr/bin/env python

############################################
#			  spec_readspec.py 			   #
# script to read a spectrum from a 1D fits #
# with multiple extenstions 			   #
############################################

import pyfits
import numpy as np

def spec_readspec_multiext(file,header='NO OUTPUT HEADER',ext_num=0):
    wave=[]
    flux=[]
    hdr=[]
    split = file.split('.')
    exten = split[-1]
    if (exten == 'fits') or (exten == 'fit'):
        hdu = pyfits.open(file)
        hdr = hdu[ext_num].header
        if 'crpix1' in hdu[ext_num].header:
            flux = hdu[ext_num].data[0,:]  # Added for FRODO spectra          
            wave = readlambda_multiext(flux,hdu,ext_num=ext_num)
        else:
            print('!!!	Wavelength keyword not found in FITS HEADER 	!!!')
            return
    else:
        print("Not yet supported!")
        return
    hdu.close	
    if header == 'NO OUTPUT HEADER':
        return wave,flux
    else:
        return wave,flux,hdr


def readlambda_multiext(spec, hdu_sp,ext_num=0):
    # For FROSOspec files, keywords are in capitals
    
    crpix1 = hdu_sp[ext_num].header['CRPIX1']
    #/ value of ref pixel
    crval1 = hdu_sp[ext_num].header['CRVAL1']
    #/ delta per pixel
    if 'cd1_1' in hdu_sp[ext_num].header:
        cd1_1 = hdu_sp[ext_num].header['CD1_1']
    #cd1_1 is sometimes called cdelt1.
    else:
        cd1_1 = hdu_sp[ext_num].header['CDELT1']

    print cd1_1
    print crval1
    print crpix1

    if cd1_1 == 0:
        print("NOT WORKING")
        return

    n_lambda = len(spec)
    wave = np.zeros(n_lambda)
    for l  in xrange(n_lambda):
        wave[l] = (l+1.0-crpix1)*cd1_1+crval1
        #Use pixel number starting with 0 if no lambda information is found.
	if (np.min(wave)+np.max(wave) == 0.0):
		print('No lambda information found: used pixel number starting with 0')
		for l  in xrange(n_lambda):
			wave[l] = l	
   
    print wave,len(wave),n_lambda
    
    return wave


