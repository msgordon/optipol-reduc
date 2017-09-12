#! /usr/bin/env python
import pyfits
import argparse
from nudge import Plotter
import os.path
from os import mkdir

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

def save_sub_list(filelist,datalist,outdir,clobber=False):
    outlist = []

    for names,data in zip(filelist,datalist):
        head1 = pyfits.getheader(names[0])
        head2 = pyfits.getheader(names[1])
        
        subp = '%s%i-%s%i_%i' % (head1['WOLLY'],int(head1['HWP']),head2['WOLLY'],int(head2['HWP']),head1['FILENUM'])

        head1['SUBP'] = (subp,'Images subtracted')
        head1['SUBP1'] = (os.path.basename(names[0]), 'Image subtraction op 1')
        head1['SUBP2'] = (os.path.basename(names[1]), 'Image subtraction op 2')
       
        outfile = ''.join([subp,'.fit'])
        print '\tWriting %s' % outfile

        outfile = os.path.join(outdir,outfile)
        pyfits.writeto(outfile,data,header=head1,clobber=clobber)

        outlist.append(outfile)

    return outlist


def main():
    parser = argparse.ArgumentParser(description='Form QU pairs on input data')
    parser.add_argument('filelist', nargs='+',help='Input files')
    parser.add_argument('-o',type=str,default='./nudged/', help='Specify output directory (default=./nudged/)')

    args = parser.parse_args()

    stripList = []
    for fname in args.filelist:
        stripList.append((fname, get_filenum(fname), get_wolly(fname), get_position(fname)))

    stripList = sorted(stripList, key = lambda x: (x[1],x[3],x[2]))

    # Get all A,B files
    Alist = [x for x in stripList if x[2] == 'A']
    Blist = [x for x in stripList if x[2] == 'B']

    # file pairs
    A0_B0 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 0 and y[3] == 0]
    print 'Subtracting A0 and B0...'
    A0_B0_datlist = sub_list(A0_B0)

    A45_B45 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 45 and y[3] == 45]
    print 'Subtracting A45 and B45...'
    A45_B45_datlist = sub_list(A45_B45)

    A22_B22 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 22 and y[3] == 22]
    print 'Subtracting A22 and B22...'
    A22_B22_datlist = sub_list(A22_B22)

    A67_B67 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 67 and y[3] == 67]
    print 'Subtracting A67 and B67...'
    A67_B67_datlist = sub_list(A67_B67)

    
    print 'Saving subtracted images to %s' % args.o
    try:
        mkdir(args.o)
    except:
        pass

    A0_B0_filelist = save_sub_list(A0_B0,A0_B0_datlist,args.o,clobber=True)
    A45_B45_filelist = save_sub_list(A45_B45,A45_B45_datlist,args.o,clobber=True)
    A22_B22_filelist = save_sub_list(A22_B22,A22_B22_datlist,args.o,clobber=True)
    A67_B67_filelist = save_sub_list(A67_B67,A67_B67_datlist,args.o,clobber=True)

    print 'Preparing nudge run'

    # Get Plotter for pair lists
    
    plotter = Plotter(A0_B0_filelist)
    
    

if __name__ == '__main__':
    main()
