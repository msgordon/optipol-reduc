#!/usr/bin/env python
import argparse
import os.path
from os import makedirs
from astropy.io import fits
from astropy.time import Time
import numpy as np
from functools import partial
from joblib import Parallel, delayed
from stack import read_mask

def register(

def main():
    parser = argparse.ArgumentParser(description='Register & align images')
    parser.add_argument('filenames',nargs='+',help='List of target files to register')
    parser.add_argument('-odir',metavar='outdir',required=True,type=str,help='Output directory for files.')
    parser.add_argument('--c',action='store_true',help='Clobber (overwrite) on output')
    parser.add_argument('-njobs',type=int,default=1,help='Process images in parallel. "-1" is all CPUs (default=1).')
    
    args = parser.parse_args()

    # create output directory
    if args.odir not in ['','.']:
        makedirs(args.odir,exist_ok=True)

if __name__ == '__main__':
    main()
