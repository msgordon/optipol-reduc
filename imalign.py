#!/usr/bin/env python
import argparse
import os.path
from os import makedirs
from astropy.io import fits
from astropy.time import Time
import numpy as np
from functools import partial
from joblib import Parallel, delayed
from image_registration import chi2_shift
from image_registration.fft_tools import shift
from shutil import copy
import warnings

def register(ref,toshift,outdir,overwrite=False):
    '''Register and shift images'''
    hdu_ref = fits.open(ref)
    hdu_shift = fits.open(toshift)

    # get offsets
    xoff,yoff,_,_ = chi2_shift(hdu_ref[0].data,hdu_shift[0].data,
                      boundary='constant')

    # shift
    hdu_shift[0].data = shift.shiftnd(hdu_shift[0].data,(-yoff,-xoff))
    hdu_shift[0].header.add_history('%s - %s' % (os.path.basename(__file__),Time(Time.now(),format='fits')))
    hdu_shift[0].header.add_history('%s - aligned to %s' % (os.path.basename(__file__),os.path.basename(ref)))
    hdu_shift[0].header.add_history('%s - shift_x: %.3f, shift_y: %.3f' % (os.path.basename(__file__),-xoff,-yoff))
    hdu_shift[0].header['ALIGNED'] = (True,'Image aligned to reference')
    hdu_shift[0].header['REFALGND'] = (ref,'Reference file for alignment')
    hdu_shift[0].header['SXOFF'] = (-xoff,'X shift')
    hdu_shift[0].header['SYOFF'] = (-yoff,'Y shift')

    outfile = os.path.join(outdir,os.path.basename(toshift))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        hdu_shift.writeto(outfile,overwrite=overwrite)
    return outfile
    
    

def main():
    parser = argparse.ArgumentParser(description='Register & align images')
    parser.add_argument('filenames',nargs='+',help='List of target files to register. Images are aligned to first in list.')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    # create output directory
    if args.odir not in ['','.']:
        makedirs(args.odir,exist_ok=True)

    # align all images to first filename
    ref = args.filenames[0]
    align = args.filenames[1:]

    imref = partial(register,ref=ref,outdir=args.odir,overwrite=args.c)
    
    outfiles = Parallel(n_jobs=args.njobs,verbose=11)(delayed(imref)(toshift=a) for a in align)

    # Write ref to outdir
    refnew = os.path.join(os.path.basename(ref))
    copy(ref,refnew)

    outfiles.append(refnew)

    print('Wrote %i files to %s' % (len(outfiles), args.odir))

if __name__ == '__main__':
    main()
