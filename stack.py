#! /usr/bin/env python
import argparse
import os.path
from os import makedirs
from ccdproc import ImageFileCollection, Combiner, CCDData
from astropy.io import fits
from astropy.time import Time
import numpy as np
from functools import partial
from warnings import warn
from joblib import Parallel, delayed

def read_mask(maskfile,fitsfile,wolly=False):
    '''Parse mask file and return masked region. If 'wolly' is True, return two masks.'''
    import pyregion
    hdu = fits.open(fitsfile)
    reg = pyregion.open(maskfile)

    # total mask
    mask = ~pyregion.get_mask(reg,hdu[0])
    
    if wolly:
        # return three masks: (top, bottom, total)
        if len(reg) != 2:
            raise AttributeError("Wolly masks must have two regions")

        # make total mask
        reg = [pyregion.ShapeList([r]) for r in reg]
        pmask = [~pyregion.get_mask(r,hdu[0]) for r in reg]
        pmask.append(mask)
        mask = pmask
    return mask

def get_files(filenames,location='.',imagetyp=None,filter=None,fexptime='*'):
    '''Gather files from input directories, choosing by image type filter
    if specified'''

    # for summary display only
    keywords = ('object','date-obs','IMAGETYP','FILTER','EXPTIME')
    
    collection = ImageFileCollection(filenames=filenames,location=location,
                                     keywords=keywords)
    
    if imagetyp in ['light','bias','dark']:
        imagetyp = '%s frame' % imagetyp
    else:
        imagetyp = '*'

    if fexptime is None:
        fexptime = '*'

    if imagetyp == 'bias frame':
        fexptime = 0.

    # filter and reread
    if filter is None:
        filelist = collection.files_filtered(IMAGETYP=imagetyp,EXPTIME=fexptime,include_path=True)
    else:
        filelist = collection.files_filtered(IMAGETYP=imagetyp,EXPTIME=fexptime,FILTER=filter,include_path=True)

    if not filelist:
        raise RuntimeError('No matching files found.')

    location = os.path.join(collection.location,os.path.dirname(filenames[0]))

    collection = ImageFileCollection(filenames=filelist,location=location,
                                     keywords=keywords)

    collection.summary.pprint()
    return collection

def get_exptime(header,key='EXPTIME'):
    '''Extract exptime from header'''
    t = header[key]
    if t == 0:
        # disallow normalization for 0 sec bias frames
        t = 1.
    return t

def combine(filelist,method,clip=None,clip_lo=None,clip_hi=None,norm=False,mask=None,normw=False):
    '''Combine input files with specified method'''
    if isinstance(filelist,ImageFileCollection):
        filelist = filelist.ccds(ccd_kwargs={'unit':'adu'})
    elif isinstance(filelist,list):
        filelist = (CCDData.read(fname,unit='adu') for fname in filelist)

    # must convert again. ccd_reader does not process masks
    if normw:
        # wolly masks have form: (top, bottom, total)
        tmask = mask[-1]
    else:
        tmask = mask

    filelist = [CCDData(f,unit=f.unit,mask=tmask) for f in filelist]

    c = Combiner(filelist)

    header = filelist[0].header
    header.add_history('stack.py - %s' % Time(Time.now(),format='fits'))

    # divide by exptime, unless normw specified
    if norm:
        if normw:
            #warn(RuntimeWarning('--norm parameter ignored in --normw mode.  Wolly halves will be normalized to 1.'))
            pass
        else:
            c.scaling = [1./get_exptime(f.header) for f in filelist]
            header['NORMALZD'] = (True,'Images normalized by exposure time')
            header.add_history('         - normalized by exptime')

    if clip == 'ccdclip':
        c.clip_extrema(int(clip_lo),int(clip_hi))
    elif clip == 'minmax':
        c.minmax_clipping(clip_lo,clip_hi)
    elif clip == 'sigclip':
        c.sigma_clipping(clip_lo,clip_hi)
    elif clip is None:
        pass
    else:
        raise NotImplementedError('clip method "%s" not implemented'%clip)

    
    if method == 'mean':
        ccd = c.average_combine()
    elif method == 'median':
        ccd = c.median_combine()
    elif method == 'sum':
        ccd = c.sum_combine()
    else:
        raise NotImplementedError('imcombine method %s not implemented'%method)
        
    header['NCOMBINE'] = (len(filelist),'Num images combined')
    header['NMETHOD'] = (method,'Image stack method')
    header.add_history('         - stacked %i images' % len(filelist))
    
    if clip:
        header['NCLIPMET'] = (clip,'Image clip method')
        header['NCLIPLO'] = (clip_lo,'Clip lo')
        header['NCLIPHI'] = (clip_hi,'Clip hi')
        header.add_history('         - clipped with method %s, clip_lo=%.3f, clip_hi=%.3f' % (clip,clip_lo,clip_hi))

    if normw:
        # normalize by wolly masks separately
        masktop = np.ma.array(ccd.data, mask = mask[0])
        maskbot = np.ma.array(ccd.data, mask = mask[1])
        ccd.data[~mask[0]] /= np.ma.mean(masktop)
        ccd.data[~mask[1]] /= np.ma.mean(maskbot)

        # set all unmasked values to nan (wolly overscan)
        ccd.data[mask[2]] = np.nan
        ccd.mask = None
        
        header['NORMALWS'] = (True,'Images normalized to 1 by Wolly split')
        header.add_history('         - normalized to 1 by Wolly split')
        #header['IMAGETYP'] = 'Flat Frame'
        
    header.add_history('         - %s combined images' % method)

    ccd.header = header
    #return ccd
    return ccd.to_hdu()


        


def main():
    mchoices = ('mean','median','sum')
    cchoices = (None,'ccdclip','minmax','sigclip')
    ichoices = (None,'light','bias','dark')

    parser = argparse.ArgumentParser(description='Stack target frames using specified mode.')
    parser.add_argument('filenames',nargs='+',help='List of files to combine.')
    parser.add_argument('-o',metavar='outfile',required=True,type=str,help='Output file')
    parser.add_argument('-datadir',default='.',type=str,help='Optionally specify data directory (default: %(default)s)')
    parser.add_argument('-m',metavar='method',choices=mchoices,default='mean',help='Specify combine method (default: %(default)s)')
    parser.add_argument('-clip',choices=cchoices,default=None,help='Specify clipping method (default: %(default)s)')
    parser.add_argument('-cval',nargs=2,type=float,default=(None,None),metavar=('clip_lo','clip_hi'),help='Specify low/high clipping values')
    parser.add_argument('-imgtype',default=None,choices=ichoices,help='Specify image type to filter files (default: %(default)s)')
    parser.add_argument('-filter',default=None,type=str,help='Specify filter or HWP position to filter files (default: %(default)s)')
    parser.add_argument('-exptime',default=None,type=float,help='Optionally filter by exposure time.')
    parser.add_argument('-maskfile',default=None,type=str,help='Specify maskfile.')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('--norm',action='store_true',help='Normalize by exptime before combining')
    parser.add_argument('--normw',action='store_true',help='Normalize to 1 by wolly split')
    parser.add_argument('--wolly',action='store_true',help='In wolly mode, HWP positions are treated separately.')
    parser.add_argument('-njobs',type=int,default=1,help='Process wolly groups in parallel. "-1" is all CPUs (default=1).')

    args = parser.parse_args()

    if args.maskfile:
        # process mask file
        mask = read_mask(args.maskfile,args.filenames[0],wolly=args.wolly)
    else:
        mask = None


    if args.wolly:
        # process HWP separately
        if args.maskfile is None:
            parser.error("--wolly requires -maskfile (e.g. masks/wolly_mask.reg")

        if args.filter:
            warn(RuntimeWarning('-filter parameter ignored in --wolly mode.  MaxIm DL currently uses the FITS "FILTER" keyword for the HWP position.'))

        # initialize partial functions for ease of parallelization later
        get_files_HWP = partial(get_files,
                                filenames=args.filenames,
                                location=args.datadir,
                                imagetyp=args.imgtype,
                                fexptime=args.exptime)
        proc_files = partial(combine,method=args.m,clip=args.clip,clip_lo=args.cval[0],clip_hi=args.cval[1],norm=args.norm,mask=mask,normw=args.normw)

        filters = ('0','45','22.5','67.5')

        # if njobs != 1, execute in parallel
        with Parallel(args.njobs) as parallel:
            # get files by HWP pos
            collections = parallel(delayed(get_files_HWP)(filter=filt) for filt in filters)
            # combine by HWP pos
            ccds = parallel(delayed(proc_files)(c) for c in collections)

            for ccd in ccds:
                #update header to include HWP pos
                ccd[0].header['HWP'] = (np.float(ccd[0].header['FILTER']),
                                        'HWP position')
        
        # create output directory
        mdir = os.path.dirname(args.o)
        if mdir not in ['','.']:
            makedirs(mdir,exist_ok=True)
        else:
            mdir = '.'

        # generate each wolly file
        base,ext = os.path.splitext(args.o)
        ofiles = ('%s_%s%s'%(base,suf,ext) for suf in ('00','45','22','67'))
        try:
            for ofile, ccd in zip(ofiles,ccds):
                ccd.writeto(ofile,overwrite=args.c)
        except OSError as e:
            raise OSError("File '%s' already exists.  Re-run with --c flag to overwrite existing files." % args.o) from e
        print('Images written to %s' % mdir)

    else:
        # not wolly mode, stack normally
        if args.njobs != 1:
            warn(RuntimeWarning('-njobs parameter ignored when not in --wolly mode. Parallel processing disabled.'))
        collection = get_files(args.filenames,location=args.datadir,
                               imagetyp=args.imgtype,filter=args.filter,
                               fexptime=args.exptime)
        
        # combine images using args
        ccd = combine(collection,method=args.m,clip=args.clip,clip_lo=args.cval[0],clip_hi=args.cval[1],norm=args.norm,mask=mask,normw=args.normw)
        

        # create output directory
        mdir = os.path.dirname(args.o)
        if mdir not in ['','.']:
            makedirs(os.path.dirname(args.o),exist_ok=True)
        try:
            #ccd.write(args.o,hdu_mask=None,hdu_uncertainty=None,overwrite=args.c)
            ccd.writeto(args.o,overwrite=args.c)
        except OSError as e:
            raise OSError("File '%s' already exists.  Re-run with --c flag to overwrite existing files." % args.o) from e

        print('Image written to %s' % args.o)
        

if __name__ == '__main__':
    main()
