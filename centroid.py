#! /usr/bin/env python
import numpy as np
import argparse
from astropy.modeling import models, fitting
import pyfits
from scipy.optimize import curve_fit,leastsq
from scipy import stats
import matplotlib.pyplot as plt

fitter = fitting.LevMarLSQFitter()
CEN_MAN_ATTEMPTS = 0
returnMe = {}

def centroid(data, coords, rad = 30, returnFit = False,manual=False):
    global CEN_MAN_ATTEMPTS
    if CEN_MAN_ATTEMPTS > 2:
        CEN_MAN_ATTEMPTS = 0
        return returnMe['key']
    if isinstance(data,str):
        print 'Reading from %s' % data
        data = pyfits.getdata(data)

    # Transpose x and y b/c reasons
    center_y,center_x = coords
    
    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]
    x,y = np.mgrid[0:dslice.shape[0],0:dslice.shape[1]]
    x -= dslice.shape[0]/2.
    y -= dslice.shape[1]/2.

    if manual:
        extent = (np.min(x),np.max(x),np.min(y),np.max(y))# (l,r,b,t)
        centroid_manual(data,dslice,extent,coords,rad,returnFit)
        CEN_MAN_ATTEMPTS = 0
        return returnMe['key']
    
    p_init = models.Gaussian2D(np.max(dslice),0,0,rad,rad)
    p = fitter(p_init,x,y,dslice)

    # Rescale coordinates to match data
    p.x_mean = center_y - p.x_mean
    p.y_mean = center_x - p.y_mean

    dist = np.sqrt((p.x_mean.value-center_y)**2 + (p.y_mean.value-center_x)**2)
    # If distance is too big, prompt to select center by hand
    if dist > rad:
        print 'Cannot locate center.  Select manually.'
        extent = (np.min(x),np.max(x),np.min(y),np.max(y))# (l,r,b,t)
        return centroid_manual(data,dslice,extent,coords,rad,returnFit)

    CEN_MAN_ATTEMPTS = 0
    if returnFit:
        return p.x_mean.value, p.y_mean.value, p
    
    else:
        return p.x_mean.value, p.y_mean.value


def centroid_manual(data,dslice,extent,coords,rad,returnFit):
    def onclick(event):
        global returnMe
        global CEN_MAN_ATTEMPTS
        xc,yc = (event.xdata,event.ydata)
        if not xc or not yc:
            print 'Click on center of target'
        else:
            xr,yr = (xc+coords[0],yc+coords[1])
            print '(%.2f,%.2f) -> (%.2f,%.2f)' %(xc,yc,xr,yr)
            plt.close()
            CEN_MAN_ATTEMPTS += 1
            returnMe['key'] = centroid(data,(xr,yr),rad,returnFit)

    print 'Attempt %i' % (CEN_MAN_ATTEMPTS+1)
    fig = plt.figure()
    plt.imshow(dslice,origin='lower',extent=extent,interpolation='none',vmin=0)
    cid = fig.canvas.mpl_connect('button_press_event',onclick)
    print 'Click on center of target'
    plt.show()
    return returnMe['key']

def centroid_airy(data, coords, rad = 30, returnFit = False):
    if isinstance(data,str):
        data = pyfits.getdata(data)

    # Transpose x and y b/c reasons
    center_y,center_x = coords
    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]

    # Construct a grid of coordinates
    x,y = np.mgrid[0:dslice.shape[0],0:dslice.shape[1]]
    x -= dslice.shape[0]/2.
    y -= dslice.shape[1]/2.
                
    p_init = models.AiryDisk2D(np.max(dslice),0,0,rad)
    p = fitter(p_init,x,y,dslice)

    # Rescale coordinates to match data
    px = center_y + p.y_0
    py = center_x + p.x_0

    if returnFit:
        return px, py, p
    
    else:
        return px, py


def r2(p, x, y):
    return p[4] * (x - p[2])**2 + p[5] * (y - p[3])**2

def moffat(p, rsq, satlevel):
    return np.clip(p[0] + p[1] / np.abs(1 + rsq)**p[6], 0, satlevel)

def errfunc(p, x, y, f, satlevel):
    rsq = r2(p, x, y)   
    return (moffat(p, rsq, satlevel) - f) 

#def moffat(xy,scale,a,b):
#    x,y = xy
#    return scale*(b-1.0)/(np.pi*a**2)*(1.0+(x**2+y**2)/a**2)**(-b)

def centroid_moffat(data,coords, rad=30):
    if isinstance(data,str):
        data = pyfits.getdata(data)

    # Transpose x and y b/c reasons
    center_y,center_x = coords
    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]
    dim = np.max(dslice.shape)

    # Generate x,y grids
    x = np.linspace(0, dim - 1., dim) - dim // 2
    y = np.linspace(0, dim - 1., dim) - dim // 2
    x, y = np.meshgrid(x, y)

    dslice = np.reshape(dslice, (-1))
    x = np.reshape(x, (-1))
    y = np.reshape(y, (-1))

    # leastsq requires a saturation level?
    satlevel = stats.scoreatpercentile(dslice, 99)
    
    #p_init = [vertshift,scaling,xcent,ycent,alpha,beta]
    p_init = [np.min(dslice),np.max(dslice),center_y,center_x,5e-4,1,8]
    p,success = leastsq(errfunc,p_init[:],args=(x,y,dslice,satlevel))
    print p
    return (p[2],p[3])


        
def main():
    parser = argparse.ArgumentParser(description='Return central coordinates of object based on initial guess.')
    parser.add_argument('data',type=str,nargs=1,help='FITS file with object to be located')
    parser.add_argument('coords',type=int,nargs=2,help='x y coordinates of initial guess')
    parser.add_argument('-rad',type=int,default=30,help='Search radius [in pixels] from initial guess.')
    parser.add_argument('--fit',action='store_true',help='Return fitting params')

    args = parser.parse_args()

    center = centroid(args.data[0],args.coords,args.rad,args.fit)

    print '%s\t%f\t%f' % (args.data[0],center[0],center[1])
    if args.fit:
        print center[2]

    return 0

if __name__ == '__main__':
    main()
