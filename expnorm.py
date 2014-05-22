#! /usr/bin/env python
import pyfits
import argparse
import numpy as np
import os


def normalize(name,outdir=None):
    f = pyfits.open(name)
    hdr = f[0].header
    exptime = hdr['EXPTIME']

    data = f[0].data / exptime

    hdr['EXPNORM'] = (True, 'Normalized by exptime')

    if outdir is None:
        outname = name
    else:
        outname = os.path.join(outdir,name)

    print 'Writing to %s' % outname
    pyfits.writeto(outname,data, header=hdr,clobber=True)
    return outname
    


def main():
    parser = argparse.ArgumentParser(description='Normalize all images by exptime.')
    parser.add_argument('filelist',nargs='+',help='List of filenames')
    parser.add_argument('-o',default=None,dest='outdir',type=str,help='Specify output directory for data files after normalization.')

    args = parser.parse_args()

    for name in args.filelist:
        normalize(name,args.outdir)


if __name__ == '__main__':
    main()
