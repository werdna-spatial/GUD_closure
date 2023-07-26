#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 12:46:00 2022

@author: eric
"""
import netCDF4
from netCDF4 import Dataset
import shutil
import sys
import argparse
from pathlib import Path

def main():
    print(args.filename, file = sys.stdout)
    print(args.DIRnc, file = sys.stdout)
    
    filedir=args.DIRnc
    filename=args.filename
    fileparts= filename.split('.')
    fileparts[1]
    rootgrp = Dataset(Path(filedir,filename), "r", format="NETCDF4")
    descparts=(rootgrp.variables[fileparts[1]].description).split()
    rootgrp.close()
    newname=fileparts[0]+'.'+fileparts[1]+'.'+descparts[0]+'.'+fileparts[2]+'.'+fileparts[3]
    shutil.copy2(Path(filedir,filename), Path(filedir,newname) )
#
if __name__ == "__main__":
    #Initialize
    parser=argparse.ArgumentParser(description="Rename NC files based on variable desc, assumes mitgcm desc")
     
    #Adding optional parameters
    parser.add_argument('-fn',
                        '--filename',
                        help="name of nc file",
                        required=True,
                        type=str)
 
    parser.add_argument('-D',
                        '--DIRnc',
                        help="Dir of netcdf files",
                        required=True,
                        type=str)
 
    args = parser.parse_args()
    main()
    