#!/usr/bin/env python
import os
import glob
from argparse import ArgumentParser
import h5py 

def main():
    psr = ArgumentParser()
    psr.add_argument('h5file',type=str)
    psr.add_argument('group',type=str)
    args = psr.parse_args()

    with h5py.File(args.h5file) as h5:
        del h5[args.group] 
    
if __name__=="__main__":
    main()

