#! /usr/bin/env python
import argparse
from ccdproc import CCDData, ccd_process, subtract_bias, flat_correct
from astropy.io import fits
import astropy.units as u
from sys import exit
from functools import partial
import os.path
from os import makedirs
from astropy.time import Time
import numpy as np
from astropy import log
log.setLevel('WARNING') # CCDData throws too much info at you

def imreduce(target,bias=None,dark=None,flat=None,expkey='EXPTIME',expunit=u.s,normkey='NORMALZD',biaskey='BIASSUB',maskfile=None,fixpix=False):
    # must read bias and dark in here since CCDData not picklable
    bias = CCDData.read(bias,unit='adu') if bias else None
    dark = CCDData.read(dark,unit='adu') if dark else None
    
    if dark:
        # if dark is already normalized, don't normalize again
        if normkey in dark.header and dark.header[normkey]:
            print(dark.header[normkey])
            dark_exposure = None
            dark_scale = False
        else:
            dark_exposure = dark.header[expkey]*expunit
            dark_scale = True

        # if bias frame, bias subtract darks first
        if bias:
            # if dark frame is already bias subtracted, skip
            if biaskey in dark.header and dark.header[biaskey]:
                pass
            else:
                dark = subtract_bias(dark,bias)
    else:
        dark_exposure = None
        dark_scale = False

    if flat:
        if isinstance(flat,str):
            # flat passed in as string, must read to CCDData
            ffile = flat
            flat = CCDData.read(ffile,unit='adu')
            flat.header['OFNAME'] = ffile

        # ccdproc has issues with nans when flatfielding
        #  -> replace nans with 1
        mask = np.isnan(flat.data)
        flat.data[mask] = 1.        

    if isinstance(target,str):
        # target passed in as string, must read to CCDData
        tfile = target
        target = CCDData.read(tfile,unit='adu')
        target.header['OFNAME'] = tfile
        
    data_exposure = target.header[expkey]*expunit

    # read in bad pixel mask if specified
    if maskfile:
        target.header.add_history('%s - masked points from %s' % (os.path.basename(__file__),maskfile))
        maskfile = fits.getdata(maskfile)


    # do the ccd reduction!
    ccd = ccd_process(target,master_bias=bias,dark_frame=dark,master_flat=flat,
                      #exposure_key=expkey,exposure_unit=expunit,
                      data_exposure=data_exposure,dark_exposure=dark_exposure,
                      dark_scale=dark_scale,bad_pixel_mask=maskfile)

    # fix bad pixels through grid interpolation
    if fixpix:
        if maskfile is None:
            warn(RuntimeWarning('No maskfile specified: --fixpix ignored'))
        else:
            from photutils.utils import interpolate_masked_data
            from photutils.utils import ShepardIDWInterpolator as idw
            from scipy.interpolate import griddata
            mask = np.asarray(maskfile,dtype=bool)
            goodpoints = np.where(~mask)
            goodvals = ccd.data[goodpoints]
            badpoints = np.where(mask)
            # interpolate and replace values
            res = griddata(goodpoints,goodvals,badpoints)
            ccd.data[badpoints] = res
            '''
            mask = np.asarray(maskfile,dtype=bool)
            goodcoord = np.where(~mask)
            print(goodcoord,len(goodcoord))
            goodval = ccd.data[goodcoord]
            print(len(goodval));exit()
            iz = idw(goodcoord,goodval)
            fixed = iz(np.where(mask))
            ccd.data[mask] = fixed
            
            #fixed,_,_ = interpolate_masked_data(ccd.data,np.asarray(maskfile,dtype=bool))
            #ccd.data = fixed
            '''
            ccd.header.add_history('%s - interpolated over bad pixels' % os.path.basename(__file__))

    return ccd.to_hdu(hdu_mask=None,hdu_uncertainty=None)


def main():
    parser = argparse.ArgumentParser(description='Perform ccd reduction on target images, honoring correct HWP position for flats, if specified.')
    parser.add_argument('filenames',nargs='+',help='List of target files to process.')
    parser.add_argument('-bias',default=None,type=str,help='Bias frame')
    parser.add_argument('-dark',default=None,type=str,help='Dark frame')
    parser.add_argument('-flat',nargs='+',default=None,help='Flat frame. If included, requires flats matching the HWP positions of target files.')
    parser.add_argument('-expkey',default='EXPTIME',type=str,help='exposure time header key (default=EXPTIME)')
    parser.add_argument('-normkey',default='NORMALZD',type=str,help='norm header key (default=NORMALZD)')
    parser.add_argument('-biaskey',default='BIASSUB',type=str,help='bias subtracted header key (default=BIASSUB)')
    parser.add_argument('-hwpkey',default='HWP',type=str,help='HWP header key (default=HWP), but sometimes needs to be FILTER')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('-maskfile',default=None,type=str,help='Specify bad pixel mask.')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('--fixpix',action='store_true',help='Interpolate over bad pixel mask.')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    # read in files and store fname in header
    targets = [CCDData.read(fname,unit='adu') for fname in args.filenames]
    for t,f in zip(targets,args.filenames):
        t.header['OFNAME'] = f
        
    #bias = CCDData.read(args.bias,unit='adu') if args.bias else None
    #dark = CCDData.read(args.dark,unit='adu') if args.dark else None

    # create output directory
    if args.odir not in ['','.']:
        makedirs(args.odir,exist_ok=True)

    flats = [CCDData.read(fname,unit='adu') for fname in args.flat] if args.flat else None

    # if flats, verify that the HWP positions are represented
    if flats:
        for flat,fname in zip(flats,args.flat):
            flat.header['OFNAME'] = fname
            
        # MAXIM sometimes writes 'HWP' as FILTER, and sometimes as 'HWP'
        #  this overrides the command line arguments for now
        if args.hwpkey in targets[0].header:
            thkey = args.hwpkey
        else:
            thkey = 'FILTER'
        targetHWP = [float(ccd.header[thkey]) for ccd in targets]

        # just in case, do same check for flats
        if args.hwpkey in flats[0].header:
            fhkey = args.hwpkey
        else:
            fhkey = 'FILTER'
        flatHWP = [float(ccd.header[fhkey]) for ccd in flats]

        # if sets are not equal, then you are missing flats
        if set(targetHWP) != set(flatHWP):
            exit('ERROR: Flats must have HWP angles matching target files')

        # match flats to target
        flatDict = {hwp:flat for hwp,flat in zip(flatHWP,flats)}

        # tuples of (targetCCD, flatCCD)
        targets = ((ccd,flatDict[thwp]) for ccd,thwp in zip(targets,targetHWP))

    else:
        # tuples of (targetCCD, None) if no flats
        targets = ((ccd,None) for ccd in targets)
        
    # finally process images
    #  generate partial function for ease of deploying to joblib
    #   bias and dark must stay as files.  CCDData not picklable
        
    rfunc = partial(imreduce,bias=args.bias,dark=args.dark,
                    expkey=args.expkey,normkey=args.normkey,
                    biaskey=args.biaskey,
                    maskfile=args.maskfile,fixpix=args.fixpix)
    
    if args.njobs != 1:
        from joblib import Parallel, delayed
        #reduced = Parallel(n_jobs=args.njobs,verbose=11)(delayed(rfunc)(target,flat=flat) for target,flat in targets)
        reduced = Parallel(n_jobs=args.njobs,verbose=51)(delayed(rfunc)(target.header['OFNAME'],flat=flat.header['OFNAME']) for target,flat in targets)

        for hdu in reduced:
            basename = os.path.basename(hdu[0].header['OFNAME'])
            basename = '.PROC'.join(os.path.splitext(basename))
            outfile = os.path.join(args.odir,basename)

            hdu[0].header['HWP'] = (float(hdu[0].header[thkey]),'HWP position angle')

            #hdu = reduced.to_hdu()
            hdu[0].header.add_history('%s - %s' % (__file__,Time(Time.now(),format='fits')))
            hdu.writeto(outfile,overwrite=args.c)

    else:
        from astropy.utils.console import ProgressBar
        for t,f in ProgressBar(list(targets)):
            hdu = rfunc(t,flat=f)
            basename = os.path.basename(t.header['OFNAME'])
            basename = '.PROC'.join(os.path.splitext(basename))
            outfile = os.path.join(args.odir,basename)

            hdu[0].header['HWP'] = (float(hdu[0].header[thkey]),'HWP position angle')

            #hdu = reduced.to_hdu()
            hdu[0].header.add_history('%s - %s' % (__file__,Time(Time.now(),format='fits')))
            hdu.writeto(outfile,overwrite=args.c)
        
    print('Wrote %i files to %s' % (len(args.filenames),args.odir))

if __name__ == '__main__':
    main()
