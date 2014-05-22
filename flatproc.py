#! /usr/bin/env python
import pyfits
import argparse
import os
import numpy as np
from pyds9 import pydisplay


def get_HWP(filename):
    #grab HWP from filter field
    hdr = pyfits.getheader(filename)
    return float(hdr['FILTER'])

def normalize(arr):
    # Normalize two halves separately
    #   mask out black parts
    ydim = np.shape(arr)[1]/2
    arr[0:ydim,:] = arr[0:ydim,:]/np.mean(arr[20:400,100:]) #B bot
    arr[ydim:,:] = arr[ydim:,:]/np.mean(arr[591:960,100:])  #A top
    
    return arr

def mean_combine(fList):
    dat = []
    for f in fList:
        datum = pyfits.getdata(f)
        dat.append(normalize(datum))

    dat = np.mean(dat,axis=0)
    
    return normalize(dat)

    
def flatcombine(filenames,outdir=None):
    fList0 = [x for x in filenames if get_HWP(x) == 0.0]
    fList45 = [x for x in filenames if get_HWP(x) == 45.0]
    fList22 = [x for x in filenames if get_HWP(x) == 22.5]
    fList67 = [x for x in filenames if get_HWP(x) == 67.5]

    f0 = mean_combine(fList0)
    f45 = mean_combine(fList45)
    f22 = mean_combine(fList22)
    f67 = mean_combine(fList67)

    retNames = []
    
    for im,HWP in zip([f0,f45,f22,f67],[0,45,22.5,67.5]):
        HWP_s = str(HWP) if HWP == 0 or HWP == 45 else str(HWP)[:2]
        if outdir is None:
            newName = os.path.dirname(filenames[0])
            newExt = os.path.splitext(filenames[0])[1]
            newName = os.path.join(newName,'MasterFlat_')
            newName = ''.join([newName,HWP_s,newExt])
        else:
            newExt = os.path.splitext(filenames[0])[1]
            newName = os.path.join(outdir,'MasterFlat_')
            newName = ''.join([newName,HWP_s,newExt])

        hdr = pyfits.getheader(filenames[0])
        hdr['FILTER'] = HWP
        hdr['HWP'] = HWP
        hdr['COMMENT'] = 'Master flatfield'

        retNames.append(newName)
        print 'Writing to %s' % newName
        
        pyfits.writeto(newName,im,header=hdr,clobber=True)
    

    return retNames


def flatfield(flatF,dataF,outdir=None):
    flat = pyfits.getdata(flatF)
    data = pyfits.getdata(dataF)

    if get_HWP(flatF) != get_HWP(dataF):
        exit('HWP positions of flat and data do not match!')
        
    newData = data / flat

    hdr = pyfits.getheader(dataF)
    hdr['FLATFILE'] = (flatF,'Flat applied')
    hdr['FLATPROC'] = (True, 'Flat-fielded')
    HWP = get_HWP(dataF)
    hdr['FILTER'] = HWP
    hdr['HWP'] = HWP
        

    if outdir is None:
        newName = dataF
    else:
        newName = os.path.join(outdir,os.path.basename(dataF))

    print 'Writing to %s' % newName
    pyfits.writeto(newName,newData,header=hdr,clobber=True)
    return newName

def flatSub(flatFiles, dataFiles, outdir=None):
    dList0 = [x for x in dataFiles if get_HWP(x) == 0.0]
    dList45 = [x for x in dataFiles if get_HWP(x) == 45.0]
    dList22 = [x for x in dataFiles if get_HWP(x) == 22.5]
    dList67 = [x for x in dataFiles if get_HWP(x) == 67.5]

    try:
        f0 = [x for x in flatFiles if get_HWP(x) == 0.0][0]
        f45 = [x for x in flatFiles if get_HWP(x) == 45.0][0]
        f22 = [x for x in flatFiles if get_HWP(x) == 22.5][0]
        f67 = [x for x in flatFiles if get_HWP(x) == 67.5][0]

    except:
        exit('Flats for each HWP position required')

    print 'Applying 0 flat'
    fList0 = [flatfield(f0,x,outdir) for x in dList0]
    print 'Applying 45 flat'
    fList45 = [flatfield(f45,x,outdir) for x in dList45]
    print 'Applying 22 flat'
    fList22 = [flatfield(f22,x,outdir) for x in dList22]
    print 'Applying 67 flat'
    fList67 = [flatfield(f67,x,outdir) for x in dList67]

    return (fList0,fList45,fList22,fList67)
    

def main():
    parser = argparse.ArgumentParser(description='Combine flat images or apply flats to input images.')
    parser.add_argument('flats',nargs='+',help='List of zero filenames')
    parser.add_argument('--data',nargs='+',default=None,help='List of data files to be flat-fielded.')
    parser.add_argument('-fout',default=None,type=str,dest='foutdir',help='Specify output master flat directory.')
    parser.add_argument('-dout',default=None,type=str,dest='doutdir',help='Specify output directory for data files after flat-fielding.')
    

    args = parser.parse_args()

    # if only one flat image, nothing to be done
    if len(args.flats) == 1 and args.data is None:
        print 'Need at least 2 flat images to combine'

    # if data is present, assume that 'flats' are master
    if args.data is not None and len(args.flats) == 4:
        print 'Applying flat field to images'
        flatSub(args.flats,args.data,args.doutdir)

    else:
        # combine flats
        print 'Combining flat frames'
        flatcombine(args.flats,args.foutdir)




if __name__ == '__main__':
    main()
