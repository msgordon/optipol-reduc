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
from astropy.io.fits.verify import VerifyWarning
import warnings
from astropy.nddata.utils import Cutout2D
from photutils.centroids import centroid_2dg as centfunc
#from astropy.modeling import models, fitting
#fit_g = fitting.LevMarLSQFitter()

def point_align(ref,toshift):
    '''Align on point source'''
    ref_cent = centfunc(ref)
    shift_cent = centfunc(toshift)

    # find centers of gaussians
    #ginit = models.Gaussian2D()
    #y,x = np.indices(ref.shape)
    #refG = fit_g(ginit,x,y,ref)
    #print(refG)
    #shiftG = fit_g(ginit,x,y,toshift)
    #print(shiftG);exit()

    xoff = shift_cent[0] - ref_cent[0]
    yoff = shift_cent[1] - ref_cent[1]
    return xoff,yoff

def extended_align(ref,toshift):
    '''Align with fft xcorr'''
    xoff,yoff,_,_ = chi2_shift(ref,toshift,boundary='constant')
    return xoff,yoff

def register(ref,toshift,method,outdir=None,center=None,size=None,overwrite=False):
    '''Register and shift images'''
    hdu_ref = fits.open(ref)
    hdu_shift = fits.open(toshift)

    if center is None:
        # use center of image
        center = [int(x/2.) for x in hdu_ref[0].shape]
    
    # resize
    if size is not None:
        refdat = Cutout2D(hdu_ref[0].data,position=center,size=size).data
        shiftdat = Cutout2D(hdu_shift[0].data,position=center,size=size).data
    else:
        refdat = hdu_ref[0].data
        shiftdat = hdu_shift[0].data
        
    if method == 'point':
        xoff,yoff = point_align(refdat,shiftdat)
    elif method == 'extended':
        xoff,yoff = extended_align(refdat,shiftdat)
    else:
        raise(NotImplementedError("method %s unknown"%method))

    # shift
    hdu_shift[0].data = shift.shiftnd(hdu_shift[0].data,(-yoff,-xoff))
    hdu_shift[0].header.add_history('%s - %s' % (os.path.basename(__file__),Time(Time.now(),format='fits')))
    hdu_shift[0].header.add_history('%s - aligned to %s' % (os.path.basename(__file__),os.path.basename(ref)))
    hdu_shift[0].header.add_history('%s - shift_x: %.3f, shift_y: %.3f' % (os.path.basename(__file__),-xoff,-yoff))
    hdu_shift[0].header['ALIGNED'] = (True,'Image aligned to reference')
    hdu_shift[0].header['REFALGND'] = (os.path.basename(ref),'Reference file for alignment')
    hdu_shift[0].header['SXOFF'] = (-xoff,'X shift')
    hdu_shift[0].header['SYOFF'] = (-yoff,'Y shift')

    if outdir:
        #write file and return filename
        outfile = os.path.join(outdir,os.path.basename(toshift))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',VerifyWarning)
            hdu_shift.writeto(outfile,overwrite=overwrite,output_verify='silentfix')
    
        return outfile
    else:
        # if outdir is None, return hdu
        return hdu_shift
    
    

def main():
    parser = argparse.ArgumentParser(description='Register & align images')
    parser.add_argument('filenames',nargs='+',help='List of target files to register. Images are aligned to first in list.')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('-m',metavar='method',choices=('point','extended'),default='extended',help='Specify alignment method (point or extended); default=extended.')
    parser.add_argument('-xy',nargs=2,type=float,default=None,help='Specify approximate "x y" pixel coordinate of object to centroid on.  Required for point mode; useful for extended mode (default=center of image).')
    parser.add_argument('-box',nargs=2,type=int,default=None,help='Specify box size (w h) to restrict alignment search.  Useful for both point & extended modes (default=full size of array).')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    if args.m == 'point' and args.xy is None:
        parser.error("-m point requires -xy coordinate")

    # create output directory
    if args.odir not in ['','.']:
        makedirs(args.odir,exist_ok=True)

    # align all images to first filename
    ref = args.filenames[0]
    align = args.filenames[1:]

    imref = partial(register,ref=ref,outdir=args.odir,
                    method=args.m,center=args.xy,size=args.box,
                    overwrite=args.c)
    
    outfiles = Parallel(n_jobs=args.njobs,verbose=11)(delayed(imref)(toshift=a) for a in align)

    # Write ref to outdir
    refnew = os.path.join(args.odir,os.path.basename(ref))
    copy(ref,refnew)

    outfiles.append(refnew)
    print('Wrote %i files to %s' % (len(outfiles), args.odir))

if __name__ == '__main__':
    main()
