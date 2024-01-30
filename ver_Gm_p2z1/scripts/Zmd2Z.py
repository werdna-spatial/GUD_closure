# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 16:32:21 2022

@author: eric
"""

import sys
import argparse
import cartopy
import pandas as pd
import xarray as xr
import datetime
import numpy as np
import os 
import geopandas as gpd
from shapely.geometry import Point
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import dask.dataframe as dd
import pathlib
from pylr2 import regress2

import cartopy.feature as cfeature
import cartopy.crs as ccrs
from scipy import stats

def main():
    
    gridfile=args.gridfilename
    print(gridfile, file = sys.stdout)
    inputfile=args.inputfilename
    print(inputfile, file = sys.stdout)
    outputfile=args.outputfilename
    print(outputfile, file = sys.stdout)
    
    #PD_NC=pathlib.Path('./')
    #os.chdir(PD_NC)
    #
    ds_grid=xr.open_mfdataset(gridfile,  combine='by_coords', parallel=False,chunks={'T':1,'Z':1})
    depth_values=np.array(ds_grid.Z.values)
    #pre='3d'
    #post='Z3d'
    I_name=inputfile
    O_name=outputfile
    ds_tracers=xr.open_mfdataset(I_name,  combine='by_coords', parallel=False,chunks={'T':1,'Z':1})
    ds_tracersZ=ds_tracers.rename({str('Zmd000023'):str('Z')})
    ds_tracersZ=ds_tracersZ.assign_coords(Z=depth_values)
    ds_tracersZ.to_netcdf(path=O_name)
    ds_tracers.close()
    ds_tracersZ.close()
    ########################
    #
    print('EXIT  :: EXIT', flush=True, file = sys.stdout)
   
    
   
#
if __name__ == "__main__":
    #Initialize
    #Initialize
    parser=argparse.ArgumentParser(description="Rename Depth attribute in NC file from Zmd to Z.  Need for comforty between tracers and diagnostics")
     
    #Adding optional parameters
    parser.add_argument('-grid',
                        '--gridfilename',
                        help="grid nc file",
                        required=True,
                        type=str)
    parser.add_argument('-fin',
                        '--inputfilename',
                        help="in name of nc file",
                        required=True,
                        type=str)
    parser.add_argument('-fout',
                     '--outputfilename',
                     help="out name of nc file",
                     required=True,
                     type=str)
 
    args = parser.parse_args()
    
    main()
