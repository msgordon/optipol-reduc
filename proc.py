#! /usr/bin/env python
import argparse


def main():
    parser = argparse.ArgumentParser(description='Perform ccd reduction on target images.')
    parser.add_argument('filenames',nargs='+',help='List of files to process.')
    parser.add_argument('-o',metavar='outfile',required=True,type=str,help='Output file')



if __name__ == '__main__':
    main()
