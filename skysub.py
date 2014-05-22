#! /usr/bin/env python
import pyfits
import argparse
import os
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pyregion


def get_region_box(reglist):
    xcen,ycen,w,h,theta = reglist
    theta = np.deg2rad(theta)
    
    # botleft,botright,topright,topleft
    box = np.array([[-w/2,-h/2],[w/2,-h/2],
                    [w/2,h/2],[-w/2,h/2],[-w/2,-h/2]])

    rot = np.array([[np.cos(theta),np.sin(theta)],
                    [-np.sin(theta),np.cos(theta)]])

    tbox = [np.dot(x,rot) for x in box]
    tbox = np.array([[x+xcen,y+ycen] for x,y in tbox])


    path = Path(tbox,closed=True)
    return path


def make_reg_masks(regfile,shape):
    r = pyregion.open(regfile)

    if len(r) != 2:
        raise Exception('Exactly two box regions required')

    paths = []
    for reg in r:
        paths.append(get_region_box(reg.coord_list))

    #Always have A be the top half
    if paths[0].get_extents().ymax > paths[1].get_extents().ymax:
        pathA = paths[0]
        pathB = paths[1]
    else:
        pathA = paths[1]
        pathB = paths[2]


    print 'Building skymasks'
    maskA = np.array([True if pathA.contains_point([x,y]) else False for x,y in np.ndindex(shape)])
    maskA = maskA.reshape(shape).T
    
    maskB = np.array([True if pathB.contains_point([x,y]) else False for x,y in np.ndindex(shape)])
    maskB = maskB.reshape(shape).T

    return (~maskA, ~maskB)

def skysub(name,maskA,maskB,outdir=None):
    f = pyfits.open(name)
    data = f[0].data
    hdr = f[0].header

    # Skysub two halves separately
    datA = np.ma.masked_array(data,maskA)
    valA = np.ma.median(datA)
    hdr['SKYVALA'] = (valA,'Median sky of A')
    hdr['SKYVARA'] = (np.ma.std(datA)**2,'Sky variance of A')
    hdr['SKYRMSA'] = (np.sqrt(1.0/datA.count()*np.ma.sum(datA**2)),'RMS sky of A')

    datB = np.ma.masked_array(data,maskB)
    valB = np.ma.median(datB)
    hdr['SKYVALB'] = (valB,'Median sky of B')
    hdr['SKYVARB'] = (np.ma.std(datB)**2,'Sky variance of B')
    hdr['SKYRMSB'] = (np.sqrt(1.0/datB.count()*np.ma.sum(datB**2)),'RMS sky of B')

    ydim = np.shape(data)[1]/2
    data[ydim:,:] = data[ydim:,:] - valA
    data[0:ydim,:] = data[0:ydim,:] - valB

    if outdir is None:
        outname = name
    else:
        outname = os.path.join(outdir,name)

    print 'Applying skymasks to %s' % outname
    pyfits.writeto(outname,data,header=hdr,clobber=True)
    return outname
    



def main():
    parser = argparse.ArgumentParser(description='Skysubtract images.')
    parser.add_argument('filelist',nargs='+',help='List of filenames')
    parser.add_argument('-reg',dest='regfile',type=str,required=True,help='File with two region boxes, one for each side of Wolly')
    parser.add_argument('-o',default=None,type=str,dest='outdir',help='Specify output directory for data files after skysubtracting.')

    args = parser.parse_args()

    # get shape for masks
    shape = pyfits.getdata(args.filelist[0]).shape
    
    maskA, maskB = make_reg_masks(args.regfile,shape)

    '''
    data = pyfits.getdata(args.filelist[0])
    dmask = np.ma.masked_array(data,maskA)
    plt.imshow(dmask,cmap='gray_r',origin='lower',interpolation='none',vmin=0,vmax=1)
    plt.colorbar()
    plt.show()
    exit()
    '''
    
    for name in args.filelist:
        skysub(name,maskA,maskB,args.outdir)


if __name__ == '__main__':
    main()
