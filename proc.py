#! /usr/bin/env python
import argparse
from ccdproc import CCDData, ccd_process
from joblib import Parallel, delayed
import astropy.units as u

def reduce(target,bias=None,dark=None,flat=None,expkey='EXPTIME',expunit=u.s,normkey='NORMALZD'):
    if dark.header['NORMALZD']:
        dark_exposure = None
        dark_scale = False
    else:
        dark_exposure = dark.header[expkey]
        dark_scale = True

    data_exposure = target.header[expkey]*expunit
    
    ccd = ccd_process(target,master_bias=bias,dark_frame=dark,master_flat=flat,
                      exposure_key=expkey,exposure_unit=expunit,
                      data_exposure=data_exposure,dark_exposure=dark_exposure,
                      dark_scale=dark_scale)
    print(ccd)

def main():
    parser = argparse.ArgumentParser(description='Perform ccd reduction on target images.')
    parser.add_argument('filenames',nargs='+',help='List of target files to process.')
    parser.add_argument('-bias',default=None,type=str,help='Bias frame')
    parser.add_argument('-dark',default=None,type=str,help='Dark frame')
    parser.add_argument('-flat',nargs='+',default=None,help='Flat frame. If included, requires flats matching the HWP positions of target files.')
    parser.add_argument('-expkey',default='EXPTIME',type=str,help='exposure time header key (default=EXPTIME)')
    parser.add_argument('-normkey',default='NORMALZD',type=str,help='norm header key (default=NORMALZD)')
    parser.add_argument('-hwpkey',default='HWP',type=str,help='HWP header key (default=HWP)')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    targets = (CCDData.read(fname,unit='adu') for fname in args.filenames)
    bias = CCDData.read(args.bias,unit='adu') if args.bias else None
    dark = CCDData.read(args.dark,unit='adu') if args.dark else None

    flats = (CCDData.read(fname,unit='adu') for fname in args.flat) if args.flat else None
    #with Parallel(args.njobs) as parallel:
        
    #bias = CCDData.read(args.bias,unit='adu') if args.bias else None
    #    dark = CCDData.read(args.dark,unit='adu') if args.dark else None

        #flats = parallel(delayed(CCDData.read)(fname,unit='adu') for fname in args.flat) if args.flat else None

    # if flats, verify that the HWP positions are represented
    targetHWP = set((float(ccd.header[args.hwpkey]) for ccd in targets))
    flatHWP = set((float(ccd.header[args.hwpkey]) for ccd in flats))
        
    exit()

if __name__ == '__main__':
    main()
