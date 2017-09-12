#! /usr/bin/env python
import pyfits
import argparse
import numpy as np
from stsci.convolve import correlate2d
from scipy.ndimage.interpolation import shift
import os.path
#import matplotlib.pyplot as plt

def xcorr(data1,data2,size=None,center=None):
    '''Cross-correlate spectra'''

    # If strings, read in arrays
    if isinstance(data1,str):
        data1 = pyfits.getdata(data1)
    if isinstance(data2,str):
        data2 = pyfits.getdata(data2)

    if not center:
        center = [int(x)/2 for x in data1.shape]

    # if size is specified, make swatches around center
    if size:
        data1 = data1[center[1]-size/2:center[1]+size/2,
                      center[0]-size/2:center[0]+size/2]

        data2 = data2[center[1]-size/2:center[1]+size/2,
                      center[0]-size/2:center[0]+size/2]

    xcorr_im = correlate2d(data2,data1, mode='constant',fft=True)

    return xcorr_im

def find_shift(xcorr_im):
    '''Find maximum of cross-correlation and calculate shift'''
    y, x = np.unravel_index(np.argmax(xcorr_im), xcorr_im.shape)
    shiftx, shifty = xcorr_im.shape[1]/2 - x, xcorr_im.shape[0]/2 - y

    return shiftx,shifty
    

def main():
    parser = argparse.ArgumentParser(description='Cross correlate images and return shift necessary to align image2 to image1')
    parser.add_argument('image1',type=str, help='FITS file of image1')
    parser.add_argument('image2',type=str, help='FITS file of image2')
    parser.add_argument('-s',metavar='size',type=int, default=None, help='Specify box size for correlation. Default is the full image, which can be very slow')
    parser.add_argument('-c',metavar=('x_cen', 'y_cen'),type=int,nargs=2, default=None,help="If '-size' specified, optionally include a center for the box region. Default is the center of image1.")
    parser.add_argument('-o',type=str,nargs='?',metavar='outfile',const='-1',default=None,help="If '-o' specified, shift image2 and write to [image2].shft.fits.  If '-o [filename]', shift image2 and write to [filename].")

    args = parser.parse_args()

    print 'Cross-correlating\n\t%s\n\t%s' % (args.image1,args.image2)
    xcorr_im = xcorr(args.image1,args.image2,size=args.s,center=args.c)

    print 'Calculating shift'
    shiftx, shifty = find_shift(xcorr_im)

    print '\t(%i, %i)' % (shiftx, shifty)

    # if outfile specified, perform shift of second image
    if args.o:
        outfile = args.o if args.o != '-1' else \
                  os.path.splitext(args.image2)[0] + '.shft.fits'

        image2, header = pyfits.getdata(args.image2, header=True)
        image2 = shift(image2, (shifty,shiftx), cval=np.nan)
        header['SHFT_REF'] = (args.image1, 'Reference image of shift')
        header['SHFT_X'] = (shiftx, 'X shift pix')
        header['SHFT_Y'] = (shifty, 'Y shift pix')

        print 'Performing shift on %s' % args.image2
        print '\tWriting to %s' % outfile
        pyfits.writeto(outfile,image2,header=header,clobber=True)
    
    return 0


if __name__ == '__main__':
    main()
