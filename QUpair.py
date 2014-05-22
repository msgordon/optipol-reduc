#! /usr/bin/env python
import pyfits
import argparse
import os
from pyds9 import pydisplay
import numpy as np
import matplotlib.pyplot as plt
import pyregion
from matplotlib.path import Path

def get_filenum(name):
    num = pyfits.getheader(name)['FILENUM']
    return int(num)

def get_wolly(name):
    wolly = pyfits.getheader(name)['WOLLY']
    return wolly

def get_position(name):
    pos = pyfits.getheader(name)['HWP']
    return int(pos)


def sub_list(listOfPairs):
    retlist = []
    for f1,f2 in listOfPairs:
        if isinstance(f1,str):
            dat1 = pyfits.getdata(f1)
            dat2 = pyfits.getdata(f2)
        else:
            dat1 = f1
            dat2 = f2
        retlist.append(dat1-dat2)
    return retlist

def add_list(listOfPairs):
    retlist = []
    for f1,f2 in listOfPairs:
        if isinstance(f1,str):
            dat1 = pyfits.getdata(f1)
            dat2 = pyfits.getdata(f2)
        else:
            dat1 = f1
            dat2 = f2
        retlist.append(dat1+dat2)
    return retlist
        

def tinbergen(Atup,Btup):
    paired = []
    for A,B in zip(Atup,Btup):
        A0,A1 = [pyfits.getdata(x[0]) for x in A]
        B0,B1 = [pyfits.getdata(x[0]) for x in B]

        R = np.sqrt(np.abs(A0/B0)/np.abs(A1/B1))
        paired.append((R - 1.0)/(R + 1.0))

    return paired

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


def add_all(listOfFiles):
    summation = []
    for f in listOfFiles:
        summation.append(pyfits.getdata(f))

    return np.sum(summation,axis=0)

def total_intensity(allList):
    Alist,Blist = allList

    A_sum = []
    for i in np.arange(0,len(Alist),4):
        curList = [Alist[i][0],Alist[i+1][0],Alist[i+2][0],Alist[i+3][0]]
        A_sum.append(add_all(curList))

    A_mean = np.mean(A_sum,axis=0)

    B_sum = []
    for i in np.arange(0,len(Blist),4):
        curList = [Blist[i][0],Blist[i+1][0],Blist[i+2][0],Blist[i+3][0]]
        B_sum.append(add_all(curList))

    B_mean = np.mean(B_sum,axis=0)

    return (A_mean, B_mean)

    
def make_reg_mask(regfile,shape):
    r = pyregion.open(regfile)

    if len(r) != 1:
        raise Exception('Exactly one box region required')

    path = get_region_box(r[0].coord_list)

    print 'Building mask'
    mask = np.array([True if path.contains_point([x,y]) else False for x,y in np.ndindex(shape)])
    mask = mask.reshape(shape).T
    
    return ~mask

def get_contour_mask(data,levels):
    shape = (data.shape[1],data.shape[0])
    contours = plt.contour(data,levels=levels)
    for contour in contours.collections:
        contour.set_visible(False)

    maxval = contours.levels[0]
    paths = contours.collections[0].get_segments()

    paths = [Path(path,closed=True) for path in paths]

    path = Path.make_compound_path(*paths)
    path = path.cleaned()

    mask = np.array([True if path.contains_point([x,y]) else False for x,y in np.ndindex(shape)])
    mask = mask.reshape(shape).T

    
    #plt.imshow(mask,interpolation='none',origin='lower')
    #plt.show()
    
    return ~mask
    
    
def main():
    parser = argparse.ArgumentParser(description='Form QU pairs on input data')
    parser.add_argument('filelist', nargs='+',help='Input files')
    parser.add_argument('-reg',dest='regfile',type=str,required=True,help='File with region box.')
    
    args = parser.parse_args()

    stripList = []
    for fname in args.filelist:
        stripList.append((fname, get_filenum(fname), get_wolly(fname), get_position(fname)))

    stripList = sorted(stripList, key = lambda x: (x[1],x[3],x[2]))

    # Get all A,B files
    Alist = [x for x in stripList if x[2] == 'A']
    Blist = [x for x in stripList if x[2] == 'B']


#    Ilist = [Alist,Blist]

    # file pairs
    A0_B0 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 0 and y[3] == 0]
    A0_B0 = sub_list(A0_B0)

    A45_B45 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 45 and y[3] == 45]
    A45_B45 = sub_list(A45_B45)

    A22_B22 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 22 and y[3] == 22]
    A22_B22 = sub_list(A22_B22)

    A67_B67 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 67 and y[3] == 67]
    A67_B67 = sub_list(A67_B67)

    Qpairs = sub_list(zip(A0_B0,A45_B45))
    Qpairs = [x/2.0 for x in Qpairs]
    Q = np.mean(Qpairs,axis=0)

    Upairs = sub_list(zip(A22_B22,A67_B67))
    Upairs = [x/2.0 for x in Upairs]
    U = np.mean(Upairs,axis=0)

    IQpairs = add_list(zip(A0_B0,A45_B45))
    IQpairs = [x/2.0 for x in IQpairs]
    IQ = np.mean(IQpairs,axis=0)
    
    IUpairs = add_list(zip(A22_B22,A67_B67))
    IUpairs = [x/2.0 for x in IUpairs]
    IU = np.mean(IUpairs,axis=0)

    q = Q/IQ
    u = U/IU
    p = np.sqrt(q**2 + u**2)
    theta = 0.5 * np.arctan2(u,q)


    levels = [0.1,0.2,0.4,0.5,1]
    cmask = get_contour_mask(p,levels)
    rmask = make_reg_mask(args.regfile,(p.shape[1],p.shape[0]))
    mask = ~cmask & ~rmask
    mask = ~mask
    
    # Pmask
    pmask = np.ma.masked_array(p,mask)
    figp = plt.figure()
    im = plt.imshow(pmask,interpolation='none',cmap='gray_r',origin='lower',vmax=levels[0])
    plt.colorbar(im)

    Quiver = plt.quiver(pmask*np.cos(theta),pmask*np.sin(theta),pivot='mid', scale=None,headlength=1,headwidth=1,color='r',edgecolor='r',linewidth=(2,),alpha=0.7)

    plt.show()
    exit()
    
    '''
    ###THIS IS GOOD
    #pydisplay(p)
    #make mask
    mask = make_reg_mask(args.regfile,(p.shape[1],p.shape[0]))
    boolmask = np.ma.zeros(p.shape,dtype=bool)
    boolmask[p < 1.0] = True

    #mask = ~mask & boolmask
    #mask = ~mask

    # Pmask
    pmask = np.ma.masked_array(p,mask)

    cmask = get_contour_mask(p)
    
    
    figp = plt.figure()

    im = plt.imshow(pmask,interpolation='none',cmap='gray_r',origin='lower')
    plt.colorbar(im)

    plt.show()
    exit()
    #p[p > 0.1] = 0.0
    
    #Quiver

    
    #Quiver = plt.quiver(pmask*np.cos(theta),pmask*np.sin(theta),pivot='mid', scale=None,headlength=1,headwidth=1,color='r',edgecolor='r',linewidth=(2,),alpha=0.7)

    plt.show()
    exit()

    '''

    #### NEED TO GET I
    #pydisplay([Q,U,IQ,IU,q,u,p])
    #exit()

    ####  image mask
    #idy,idx = np.ogrid[0:p.shape[0],0:p.shape[1]]
    #mask = (idy < 28) & (idx > 3) & (idy > 1)
    #mask = ~mask
    mask = np.ma.ones(p.shape,dtype=bool)
    mask[p < 1.0] = False

    # Quiver
    figp = plt.figure()
    pmask = np.ma.masked_array(p,mask)
    im = plt.imshow(pmask,interpolation='none',cmap='gray_r',origin='lower')
    plt.colorbar(im)

    
    Quiver = plt.quiver(pmask*np.cos(theta),pmask*np.sin(theta),pivot='mid', scale=None,headlength=1,headwidth=1,color='r',edgecolor='r',linewidth=(2,),alpha=0.7)

    plt.show()
    exit()

    A0_45list = zip([x for x in Alist if x[3] == 0],[x for x in Alist if x[3] == 45])
    A22_67list = zip([x for x in Alist if x[3] == 22],[x for x in Alist if x[3] == 67])

    B0_45list = zip([x for x in Blist if x[3] == 0],[x for x in Blist if x[3] == 45])
    B22_67list = zip([x for x in Blist if x[3] == 22],[x for x in Blist if x[3] == 67])

    print Alist
    print Blist


    exit()

    print A0_45list[0],A22_67list[0]
    thing = map(pyfits.getdata,[A0_45list[0][0][0],A0_45list[0][1][0],A22_67list[0][0][0],A22_67list[0][1][0]])
    totalP = np.sum(thing,axis=0)
    pyfits.writeto('totalpTEST.fits',totalP)
    exit()
    
    ####  Q, U
    AQlist = sub_list(A0_45list)
    AUlist = sub_list(A22_67list)
    BQlist = sub_list(B0_45list)
    BUlist = sub_list(B22_67list)

    AQ = np.mean(AQlist,axis=0)
    AU = np.mean(AUlist,axis=0)
    BQ = np.mean(BQlist,axis=0)
    BU = np.mean(BUlist,axis=0)

    Q = np.mean([AQ,BQ],axis=0)
    U = np.mean([AU,BU],axis=0)
    #pyfits.writeto('QTEST.fits',Q)
    #pyfits.writeto('UTEST.fits',U)
    

    ####  q, u
    q = tinbergen(A0_45list,B0_45list)
    u = tinbergen(A22_67list,B22_67list)

    q = np.mean(q,axis=0)
    u = np.mean(u,axis=0)
    pyfits.writeto('qTEST.fits',q)
    pyfits.writeto('uTEST.fits',u)

    ####  p, theta
    p = np.sqrt(q**2 + u**2)
    theta = 0.5 * np.arctan2(u,q)
    #pyfits.writeto('pTEST.fits',p)

    exit()


    ######  regions
    r = pyregion.open('pBoxes.reg')
    paths = []
    for reg in r:
        paths.append(get_region_box(reg.coord_list))

    #print paths[0].get_extents().get_points()
    #exit()

    gridy,gridx = np.mgrid[0:p.shape[0],0:p.shape[1]]
    #for col, row in zip(gridx,gridy):
    #    grid.append(zip(col,row))
    grid = [zip(col,row) for col,row in zip(gridx,gridy)]
    points = paths[0].contains_points()
    print points
    exit()
    mask = np.unravel_index(mask,p.shape)
    pmask = np.ma.masked_array(p,mask)
    plt.imshow(pmask,interpolation='none',cmap='gray_r',origin='lower')
    plt.show()
    

    
    exit()
    bbox = paths[0].get_extents()
    bbox = [bbox.xmin,bbox.xmax,bbox.ymin,bbox.ymax]
    plt.imshow(p,extent=bbox,interpolation='none',cmap='gray_r',origin='lower')
    plt.show()
    


    exit()







    
    ####  image mask
    idy,idx = np.ogrid[0:p.shape[0],0:p.shape[1]]
    mask = (idy < 28) & (idx > 3) & (idy > 1)
    mask = ~mask
    #mask = np.ma.ones(p.shape,dtype=bool)
    #mask[p < 0.12] = False


    
    # Quiver
    figp = plt.figure()
    pmask = np.ma.masked_array(p,mask)
    im = plt.imshow(pmask,interpolation='none',cmap='gray_r',origin='lower')
    plt.colorbar(im)

    #Quiver = plt.quiver(pmask*np.cos(theta),pmask*np.sin(theta),pivot='mid', scale=0.5,headlength=1,headwidth=1,color='r',edgecolor='r',linewidth=(2,),alpha=0.7)

    #plt.quiverkey(Quiver,0,0,0.05,'5%',fontproperties={'weight': 'bold'},coordinates='data',color='r')

    # Ip
    IP = np.sqrt(Q**2 + U**2)
    IPmask = np.ma.masked_array(IP,mask)

    plt.figure()
    im = plt.imshow(IPmask,interpolation='none',cmap='gray_r',origin='lower',vmax = 200)
    plt.colorbar()

    #pyfits.writeto('IpTEST.fits',IP)
    
    
    

    
    plt.show()


    exit()
    
    AQlist = sub_list(A0_45list)
    AUlist = sub_list(A22_67list)

    BQlist = sub_list(B0_45list)
    BUlist = sub_list(B22_67list)

    
    Aqlist = []
    for x,y in zip(sub_list(A0_45list),add_list(A0_45list)):
        Aqlist.append(x/y)

    Aulist = []
    for x,y in zip(sub_list(A22_67list),add_list(A22_67list)):
        Aulist.append(x/y)

    #Aq = np.mean(Aqlist,axis=0)
    #Au = np.mean(Aulist,axis=0)
    Aplist = []
    for q,u in zip(Aqlist,Aulist):
        Aplist.append(np.sqrt(q**2 + u**2))

    #Ap = np.sqrt(Aq**2 + Au**2)
    Ap = np.mean(Aplist,axis=0)

    AQdisp = np.sum(AQlist,axis=0)
    AUdisp = np.sum(AUlist,axis=0)
    
    pydisplay([AQdisp,AUdisp,Ap])
    
    



if __name__ == '__main__':
    main()
