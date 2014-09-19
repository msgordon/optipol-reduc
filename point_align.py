#! /usr/bin/env python
import pyfits
from scipy.ndimage.interpolation import shift
import numpy as np
from centroid import centroid
import argparse
import os

def imshift(filename,shifts,center,refFile,name_ext='.al',clobber=False):
    f = pyfits.open(filename)

    header = f[0].header
    header['REF_FILE'] = (os.path.basename(refFile),'Reference file')
    header['PRE_FILE'] = (os.path.basename(filename),'Filename before shift')
    header['XSHIFT'] = (shifts[0],'X shift from ref_file')
    header['YSHIFT'] = (shifts[1],'Y shift from ref_file')
    header['XCEN'] = (center[0],'X shift from ref_file')
    header['YCEN'] = (center[1],'Y shift from ref_file')
    header['PALIGN'] = (True,'Aligned')

    newName = os.path.splitext(filename)
    newName = ''.join([newName[0],name_ext,newName[1]])

    if shifts[0] != 0 and shifts[1] != 0:
        newDat = shift(f[0].data,(shifts[0],shifts[1]))
    else:
        newDat = f[0].data
        
    print filename
    print '\tShifting (%.2f,%.2f) pixels' % (shifts[0],shifts[1])
    print '\tWriting to %s' % newName
    pyfits.writeto(newName,newDat,header=header,clobber=clobber)

    return newName

    
def align(filelist,coords=None,rad=None,aligned=False,manual=False):
    if len(filelist) < 2:
        raise InputError('At least 2 input images required')
    
    if not coords:
        h = pyfits.getheader(filelist[0])
        coords = [int(h['XCEN']),int(h['YCEN'])]

    centers = []
    if not aligned:
        # Centroid on all images
        print 'Centroiding...'
        for filename in filelist:
            center = centroid(filename,coords,rad,manual=manual)
            print center
            print '%s: (%.2f,%.2f)' % (filename,center[0],center[1])
            centers.append(center)

    else:
        # Read centers
        for filename in filelist:
            h = pyfits.getheader(filename)
            centers.append([int(h['XCEN']),int(h['YCEN'])])
        

    # Shift to center of first image
    centers = np.array(centers)
    shifts = centers-centers[0]
    if aligned:
        shifts = [(-x[1],-x[0]) for x in shifts]
    print
    print 'ref: %s (%.2f,%.2f)' % (filelist[0],centers[0][0],centers[0][1])
    print

    newNames = []
    for filename,shift in zip(filelist,shifts):
        if aligned:
            name_ext = ''
        else:
            name_ext = '.al'
        newName = imshift(filename,shift,centers[0],filelist[0],name_ext,clobber=True)
        newNames.append(newName)

    return zip(newNames,centers,shifts)
        


def main():
    parser = argparse.ArgumentParser(description='Align input images given an initial guess.')
    
    parser.add_argument('filelist',nargs='+',help='List of filenames for alignment.')
    parser.add_argument('-coords',nargs=2,type=int,help='Initial guess of coordinates: x y (default is XCEN,YCEN keywords in header)')
    parser.add_argument('--r',type=int,default=40,help='Search radius (default is 40 pixels)')
    parser.add_argument('--aligned',action='store_true',help='Specify if images have already been aligned. Align to XCEN,YCEN.')
    parser.add_argument('--man',action='store_true',help='Select centers manually')
    
    args = parser.parse_args()
    align(args.filelist,args.coords,args.r,args.aligned,args.man)

    return 0

if __name__ == '__main__':
    main()
