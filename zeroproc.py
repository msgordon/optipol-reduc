#! /usr/bin/env python
import pyfits
import argparse
import os
import numpy as np

def zerocombine(filenames,outfile=None):
    zeroIm = []
    for fname in filenames:
        zeroIm.append(pyfits.getdata(fname))

    zeroIm = np.mean(zeroIm,axis=0)
    zeroHdr = pyfits.getheader(filenames[0])
    zeroHdr['NUMBIAS'] = (len(filenames),'Num of averaged bias frames')
    zeroHdr['COMMENT'] = 'MasterBias frame'

    if outfile is None:
        zeroPath = os.path.dirname(filenames[0])
        zeroExt = os.path.splitext(filenames[0])[1]
        if zeroPath[-1] != '/':
            zeroPath = ''.join([zeroPath,'/'])
        zeroName = ''.join([zeroPath,'MasterBias',zeroExt])

    else:
        zeroName = outfile
    
    print 'Writing master bias frame to %s' % zeroName
    pyfits.writeto(zeroName,zeroIm,header=zeroHdr,clobber=True)

    return zeroName

def zeroSub(zeroFilename, listOfImages,outdir=None):
    zeroIm = pyfits.getdata(zeroFilename)

    outnames = []
    for fname in listOfImages:
        data = pyfits.getdata(fname)
        newData = data - zeroIm

        if outdir is None:
            newName = os.path.splitext(fname)
            newName = ''.join([newName[0],'_zsub',newName[1]])
        else:
            newName = os.path.join(outdir,os.path.basename(fname))

        hdr = pyfits.getheader(fname)
        hdr['zeroproc'] = (True,'Bias subtracted')
        outnames.append(newName)

        print 'Writing to %s' % newName
        pyfits.writeto(newName,newData,header=hdr,clobber=True)
        
    return outnames
    

def main():
    parser = argparse.ArgumentParser(description='Combine zero-bias images or apply zero to input images.')
    parser.add_argument('zeros',nargs='+',help='List of zero filenames')
    parser.add_argument('--data',nargs='+',default=None,help='List of data files to be zero-subtracted.')
    parser.add_argument('-zout',default=None,type=str,dest='outfile',help='Specify output master bias.')
    parser.add_argument('-dout',default=None,type=str,dest='outdir',help='Specify output directory for data files after subtraction.')
    

    args = parser.parse_args()

    # if only one zero image, nothing to be done
    if len(args.zeros) == 1 and args.data is None:
        print 'Need at least 2 zero images to combine'

    # if only one zero image, apply to data
    elif len(args.zeros) == 1 and args.data is not None:
        print 'Applying bias frame to images'
        zeroSub(args.zeros[0],args.data,args.outdir)

    else:
        # combine zeros
        print 'Combining bias frames'
        zerocombine(args.zeros,args.outfile)




if __name__ == '__main__':
    main()
