#! /usr/bin/env python
import argparse
import os.path
from os import makedirs
from ccdproc import ImageFileCollection, Combiner, CCDData
from astropy.io import fits
from astropy.time import Time

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

    # filter and reread
    if filter is None:
        filelist = collection.files_filtered(IMAGETYP=imagetyp,EXPTIME=fexptime,include_path=True)
    else:
        filelist = collection.files_filtered(IMAGETYP=imagetyp,EXPTIME=fexptime,FILTER=filter,include_path=True)

    location = os.path.join(collection.location,os.path.dirname(filenames[0]))

    collection = ImageFileCollection(filenames=filelist,location=location,
                                     keywords=keywords)

    collection.summary.pprint()
    return collection

def get_exptime(header,key='EXPTIME'):
    t = header[key]
    if t == 0:
        # disallow normalization for 0 sec bias frames
        t = 1.
    return t

def combine(filelist,method,clip=None,clip_lo=None,clip_hi=None,norm=False):
    '''Combine input files with specified method'''
    if isinstance(filelist,ImageFileCollection):
        filelist = list(filelist.ccds(ccd_kwargs={'unit':'adu'}))
    elif isinstance(filelist,list):
        filelist = [CCDData.read(fname,unit='adu') for fname in filelist]

    c = Combiner(filelist)

    if norm:
        c.scaling = [1./get_exptime(f.header) for f in filelist]

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

    header = filelist[0].header
    header.add_history('stack.py - %s' % Time(Time.now(),format='fits'))
    if norm:
        header['NORMALZD'] = (True,'Images normalized by exposure time')
        header.add_history('         - normalized by exptime')
        
    header['NCOMBINE'] = (len(filelist),'Num images combined')
    header['NMETHOD'] = (method,'Image stack method')
    header.add_history('         - stacked %i images' % len(filelist))
    if clip:
        header['NCLIPMET'] = (clip,'Image clip method')
        header['NCLIPLO'] = (clip_lo,'Clip lo')
        header['NCLIPHI'] = (clip_hi,'Clip hi')
        header.add_history('         - clipped with method %s, clip_lo=%.3f, clip_hi=%.3f' % (clip,clip_lo,clip_hi))
        
    header.add_history('         - %s combined images' % method)

    ccd.header = header
    
    return ccd

        


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
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('--norm',action='store_true',help='Normalize by exptime before combining')

    args = parser.parse_args()

    collection = get_files(args.filenames,location=args.datadir,
                           imagetyp=args.imgtype,filter=args.filter,
                           fexptime=args.exptime)

    ccd = combine(collection,method=args.m,clip=args.clip,clip_lo=args.cval[0],clip_hi=args.cval[1],norm=args.norm)

    makedirs(os.path.dirname(args.o),exist_ok=True)
    try:
        ccd.write(args.o,hdu_mask=None,hdu_uncertainty=None,overwrite=args.c)
    except OSError as e:
        raise OSError("File '%s' already exists.  Re-run with --c flag to overwrite existing files." % args.o) from e
    
    print('Image written to %s' % args.o)
        

if __name__ == '__main__':
    main()
