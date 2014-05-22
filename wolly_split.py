#! /usr/bin/env python
import pyfits
import os
import argparse
import numpy as np


def get_filenum(name,prefix):
    stripped = os.path.basename(name).split(prefix)
    num = stripped[1].split('_')[0]
    return int(num)

def split(filename, outdir, prefix):
    basename, ext = os.path.splitext(filename)
    basename = os.path.basename(basename)
    basename = os.path.join(outdir,basename)

    f = pyfits.open(filename)

    ydim = f[0].header['NAXIS2']/2

    Adat = f[0].data[ydim:,:]   #top
    Bdat = f[0].data[0:ydim,:]  #bot

    Afile = ''.join([basename,'_A',ext])
    Bfile = ''.join([basename,'_B',ext])

    hdr = f[0].header
    hdr['WOLLY'] = ('A','Image half')
    hdr['FILENUM'] = (get_filenum(filename,prefix),'Observation number')
    
    pyfits.writeto(Afile,Adat,header=hdr,clobber=True)

    hdr['WOLLY'] = ('B','Image half')
    pyfits.writeto(Bfile,Bdat,header=hdr,clobber=True)

    return (Afile, Bfile)
    


def main():
    parser = argparse.ArgumentParser(description='Split OptiPol images into A and B components.')

    parser.add_argument('file',nargs='+',help='Input file(s)')
    parser.add_argument('-o',metavar='outdir',dest='outdir',required=True,help='Output directory')
    parser.add_argument('-prefix',required=True,help='Filename prefix before filenumber.')

    args = parser.parse_args()

    for filename in args.file:
        print split(filename,args.outdir,args.prefix)

    print 'Split %i files to %s' % (len(args.file),args.outdir)
            

if __name__ == '__main__':
    main()
