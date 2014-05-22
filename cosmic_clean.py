#! /usr/bin/env python
import pyfits
import LAcosmics
import argparse


def main():
    parser = argparse.ArgumentParser(description='LACosmics cleaning on input files.')

    parser.add_argument('file',nargs='+',help='Input file(s)')

    args = parser.parse_args()

    for f in args.file:
        raw = pyfits.open(f)
        raw[0].header['LACosmic'] = (True,'Data cleaned using LACosmic')

        print 'Cleaning %s ...' % f

        clean = LAcosmics.clean_data(raw[0].data)
        pyfits.writeto(f,clean,header=raw[0].header,clobber=True)
    


if __name__ == '__main__':
    main()
