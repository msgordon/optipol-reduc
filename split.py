#!/usr/bin/env python
import argparse
import os.path
from os import makedirs
from astropy.io import fits
from astropy.time import Time
import numpy as np
from functools import partial
from joblib import Parallel, delayed
from stack import read_mask

def split_by_mask(filename,mask,outdir,overwrite=False):
    '''Read fits file, split into Ordinary (A) and Extraordinary (B) beams'''
    hdu = fits.open(filename)
    hdu[0].header['OFNAME'] = (os.path.basename(filename),'Original file')
    hdu[0].header.add_history('%s - %s' % (os.path.basename(__file__),Time(Time.now(),format='fits')))
    hdu[0].header.add_history('%s - split into wolly beams' % os.path.basename(__file__))
    
    # mask is (top,bottom,total,shape_top,shape_bot,shape_total)
    A = hdu[0].data[~mask[0]].reshape(mask[3])
    B = hdu[0].data[~mask[1]].reshape(mask[4])

    # measure size difference and pad to same shape
    diff = np.diff(np.vstack([B.shape,A.shape]),axis=0)[0]

    if diff[0] < 0:
        # A is shorter than B
        A = np.pad(A,((-diff[0],0),(0,0)),
                   mode='constant',constant_values=np.nan)
    else:
        # B is shorter than A
        B = np.pad(B,((diff[0],0),(0,0)),
                   mode='constant',constant_values=np.nan)

    if diff[1] < 0:
        # B is wider than A
        A = np.pad(A,((0,0),(-diff[1],0)),
                   mode='constant',constant_values=np.nan)
    else:
        # A is wider than B
        B = np.pad(B,((0,0),(diff[1],0)),
                   mode='constant',constant_values=np.nan)

    hduA = fits.HDUList(fits.PrimaryHDU(data=A,header=hdu[0].header))
    hduB = fits.HDUList(fits.PrimaryHDU(data=B,header=hdu[0].header))

    hduA[0].header['WOLLY'] = ('A','Image half // A = Or, B = Ex')
    hduB[0].header['WOLLY'] = ('B','Image half // A = Or, B = Ex')

    base = os.path.splitext(os.path.basename(filename))
    oA = os.path.join(outdir,''.join([base[0],'_A',base[1]]))
    oB = os.path.join(outdir,''.join([base[0],'_B',base[1]]))

    hduA.writeto(oA,overwrite=overwrite)
    hduB.writeto(oB,overwrite=overwrite)
    return oA,oB

def main():
    parser = argparse.ArgumentParser(description='Split images using wollaston mask')
    parser.add_argument('filenames',nargs='+',help='List of target files to split')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('-maskfile',required=True,type=str,help='Specify wollaston mask.')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    # create output directory
    if args.odir not in ['','.']:
        makedirs(args.odir,exist_ok=True)

    # read in mask file, get two masks back
    mask = read_mask(args.maskfile,args.filenames[0],
                     wolly=True,return_shape=True)

    sfunc = partial(split_by_mask,mask=mask,outdir=args.odir,overwrite=args.c)
    
    outfiles = Parallel(n_jobs=args.njobs,verbose=11)(delayed(sfunc)(fname) for fname in args.filenames)

    print('Wrote %i files to %s' % (len(outfiles)*2, args.odir))

if __name__ == '__main__':
    main()
