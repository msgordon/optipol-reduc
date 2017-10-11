#! /usr/bin/env python
# wrapper to cosmics.py
# fallback option in clean.py if astroscappy throws errors
from cosmics import cosmicsimage
from ccdproc import CCDData

#Perform Laplacian cosmic ray subtraction
def cosmicray_lacosmic(im,sigclip=5,sigfrac=0.3,readnoise=1,gain=1,objlim=5,satlevel=65535.0,niter=5,verbose=False):

    if isinstance(im,CCDData):
        cosmic = cosmicsimage(im.data,gain=gain,readnoise=readnoise,sigclip=sigclip,sigfrac=sigfrac,objlim=objlim,satlevel=satlevel,verbose=verbose)
        cosmic.run(maxiter=niter,verbose=verbose)
        im.data = cosmic.cleanarray
        im.mask = cosmic.mask
        im = im.to_hdu()
        
    else:
        cosmic = cosmicsimage(im,gain=gain,readnoise=readnoise,sigclip=sigclip,sigfrac=sigfrac,objlim=objlim,satlevel=satlevel,verbose=verbose)
        cosmic.run(maxiter=niter,verbose=verbose)
        im = cosmic.cleanarray

    return im
