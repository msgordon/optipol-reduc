#! /usr/bin/env python
import argparse
from ccdproc import CCDData, cosmicray_lacosmic
from joblib import Parallel, delayed
import os.path

def main():
    parser = argparse.ArgumentParser(description='Perform LACosmic cleaning of images')
    parser.add_argument('filenames',nargs='+',help='List of files to clean.')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('-sclip',metavar='sigclip',type=float,default=5,help='Laplacian-to-noise limit for cosmic ray detection. Lower values will flag more pixels as cosmic rays (default=5).')
    parser.add_argument('-sfrac',metavar='sigfrac',type=float,default=0.3,help='Fractional detection limit for neighboring pixels. For cosmic ray neighbor pixels, a Laplacian-to-noise detection limit of sigfrac * sigclip will be used. (default=0.3).')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "0" is all CPUs (default=1).')

    args = parser.parse_args()

    with Parallel(args.njobs) as parallel:
        ccds = parallel(delayed(CCDData.read)(fname,unit='adu') for fname in args.filenames)
        cleaned = parallel(delayed(cosmicray_lacosmic)(ccd,sigclip=args.sclip,sigfrac=args.sfrac) for ccd in ccds)

    print(cleaned)
    
    


if __name__ == '__main__':
    main()
