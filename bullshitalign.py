#! /usr/bin/env python
import numpy as np
import argparse
import pyfits
from scipy.ndimage.interpolation import shift as sc_shift
#from shift import shift as fft_shift
import os
from pyds9 import pydisplay


def shift(filename, xs, ys, refFile,noShift=False):
    f = pyfits.open(filename)

    header = f[0].header
    header['REF_FILE'] = (os.path.basename(refFile),'Reference file')
    header['PRE_FILE'] = (os.path.basename(filename),'Filename before shift')
    header['XSHIFT'] = (xs,'X shift from ref_file')
    header['YSHIFT'] = (ys,'Y shift from ref_file')

    newName = os.path.splitext(filename)
    newName = ''.join([newName[0],'_s',newName[1]])

    #return newName
    if noShift:
        newDat = f[0].data
    else:
        #newDat = fft_shift(f[0].data,xs,ys)
        newDat = sc_shift(f[0].data,[ys,xs])

    print 'Writing to %s' % newName
    pyfits.writeto(newName,newDat,header=header,clobber=True)

    return newName


def main():
    parser = argparse.ArgumentParser(description='Shift images to align objects in input file')

    parser.add_argument('file',help='Input file with coordinates')
    
    args = parser.parse_args()

    data = np.genfromtxt(args.file,names=['fname','x','y'],dtype=['a100','f8','f8'],autostrip=True)


    # Copy reference file
    checkMe = []
    ref = data[0]
    checkMe.append(shift(ref['fname'],0,0,ref['fname'],noShift=False))
    
    for dat in data[1:]:
        xs = ref['x'] - dat['x']
        ys = ref['y'] - dat['y']
        checkMe.append(shift(dat['fname'],xs,ys,ref['fname']))


    #pydisplay(checkMe)




if __name__ == '__main__':
    main()










