#! /usr/bin/env python
import cosmics


def clean_data(data,maxiter=5,gain=2.2,readnoise=10.0,sigclip=5,
               sigfrac=0.3,objlim=5.0,satlevel=50000,getmask=False,verbose=False):

    # Perform Laplacian cosmic ray subtraction
    cosmic = cosmics.cosmicsimage(data,gain=gain,readnoise=readnoise,sigclip=sigclip,sigfrac=sigfrac,objlim=objlim,satlevel=satlevel,verbose=verbose)
    cosmic.run(maxiter=maxiter,verbose=verbose)

    if getmask:
        return cosmic.cleanarray,cosmic.mask
    else:
        return cosmic.cleanarray
