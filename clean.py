#! /usr/bin/env python
import argparse
from ccdproc import CCDData
from joblib import Parallel, delayed
import os.path
from astropy.time import Time
try:
    from astroscrappy import detect_cosmics
    from ccdproc import cosmicray_lacosmic
    CLNMTHD = 'ccdproc.cosmicray_lacosmic'
except ImportError:
    from cosmics_wrapper import cosmicray_lacosmic
    CLNMTHD = 'cosmics.cosmicray_lacosmic'
    

def main():
    parser = argparse.ArgumentParser(description='Perform LACosmic cleaning of images')
    parser.add_argument('filenames',nargs='+',help='List of files to clean.')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    #parser.add_argument('-mode',choices=['lacosmic','median'],default='lacosmic',help='Specify mode of operation (default=lacosmic)')
    parser.add_argument('-sclip',metavar='sigclip',type=float,default=5,help='Laplacian-to-noise limit for cosmic ray detection. Lower values will flag more pixels as cosmic rays (default=5).')
    parser.add_argument('-sfrac',metavar='sigfrac',type=float,default=0.3,help='Fractional detection limit for neighboring pixels. For cosmic ray neighbor pixels, a Laplacian-to-noise detection limit of sigfrac * sigclip will be used. (default=0.3).')
    parser.add_argument('-objlim',type=float,default=5,help='Minimum contrast between Laplacian image and the fine structure image. Increase this value if cores of bright stars are flagged as cosmic rays (default=5).')
    parser.add_argument('-satlevel',type=float,default=65535,help='Saturation level of the image (electrons). This value is used to detect saturated stars and pixels at or above this level are added to the mask (default=65535)')
    parser.add_argument('-niter',type=int,default=5,help='umber of iterations of the LA Cosmic algorithm to perform (default=5).')
    #parser.add_argument('-thresh',metavar='threshold',type=float,default=5,help='Threshold for detecting cosmic rays [median] (default=5).')
    #parser.add_argument('-mbox',type=float,default=11,help='Median box for detecting cosmic rays [mbox] (default=11).')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')

    args = parser.parse_args()

    ccds = (CCDData.read(fname,unit='adu') for fname in args.filenames)
    with Parallel(args.njobs,verbose=11) as parallel:
        cleaned = parallel(delayed(cosmicray_lacosmic)(ccd,sigclip=args.sclip,sigfrac=args.sfrac,niter=args.niter,objlim=args.objlim,satlevel=args.satlevel) for ccd in ccds)

    outfiles = (os.path.join(args.odir,os.path.basename(fname)) for fname in args.filenames)
    for hdu,outfile in zip(cleaned,outfiles):
        header = hdu[0].header
        header.add_history('clean.py - %s' % Time(Time.now(),format='fits'))
        header['CLEANED'] = (True,'Cleaned with LACosmics')
        header['CLNMTHD'] = (CLNMTHD,'Method used to clean')
        hdu.writeto(outfile,overwrite=args.c)
        
    
    


if __name__ == '__main__':
    main()
