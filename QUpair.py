#! /usr/bin/env python
import pyfits
import argparse
from nudge import Plotter

def get_filenum(name):
    num = pyfits.getheader(name)['FILENUM']
    return int(num)

def get_wolly(name):
    wolly = pyfits.getheader(name)['WOLLY']
    return wolly

def get_position(name):
    pos = pyfits.getheader(name)['HWP']
    return int(pos)

def main():
    parser = argparse.ArgumentParser(description='Form QU pairs on input data')
    parser.add_argument('filelist', nargs='+',help='Input files')

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
    A0_B0 = sub_list(A0_B0)

    A45_B45 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 45 and y[3] == 45]
    A45_B45 = sub_list(A45_B45)

    A22_B22 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 22 and y[3] == 22]
    A22_B22 = sub_list(A22_B22)

    A67_B67 = [(x[0],y[0]) for x,y in zip(Alist,Blist) if x[3] == 67 and y[3] == 67]
    A67_B67 = sub_list(A67_B67)


    # Get Plotter for pair lists
    
    
    

if __name__ == '__main__':
    main()
