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

    PD_NC=pathlib.Path('./')
    os.chdir(PD_NC)
    #
    ds_grid=xr.open_mfdataset('grid.v4c.nc',  combine='by_coords', parallel=False,chunks={'T':1,'Z':1})
    depth_values=np.array(ds_grid.Z.values)
    pre='3d'
    post='Z3d'
    files_par=[\
        #'.BP.v4c.nc',\
        '.PP.Primary.v4c.nc',\
        #'.PAR01.v4c.nc'\
        '.TRAC01.DIC.v4c.nc',\
        '.TRAC02.NH4.v4c.nc',\
        '.TRAC03.NO2.v4c.nc',\
        '.TRAC04.NO3.v4c.nc',\
        '.TRAC05.PO4.v4c.nc',\
        '.TRAC06.SiO2.v4c.nc',\
        '.TRAC07.FeT.v4c.nc',\
        '.TRAC08.DOC.v4c.nc',\
        '.TRAC09.DON.v4c.nc',\
        '.TRAC10.DOP.v4c.nc',\
        '.TRAC11.DOFe.v4c.nc',\
        '.TRAC12.POC.v4c.nc',\
        '.TRAC13.PON.v4c.nc',\
        '.TRAC14.POP.v4c.nc',\
        '.TRAC15.POSi.v4c.nc',\
        '.TRAC16.POFe.v4c.nc',\
        '.TRAC17.PIC.v4c.nc',\
        '.TRAC18.ALK.v4c.nc',\
        '.TRAC19.O2.v4c.nc',\
        '.TRAC20.Euk01.v4c.nc',\
        '.TRAC21.Syn01.v4c.nc',\
        '.TRAC22.Zoo01.v4c.nc',\
        '.TRAC23.Zoo02.v4c.nc'\
        ]
    #
    files_dia=[\
        #'.BP.v4c.nc',\
        '.PP.Primary.v4c.nc',\
        #'.PAR01.v4c.nc'\
        '.TRAC01.DIC.v4c.nc',\
        '.TRAC02.NH4.v4c.nc',\
        '.TRAC03.NO2.v4c.nc',\
        '.TRAC04.NO3.v4c.nc',\
        '.TRAC05.PO4.v4c.nc',\
        '.TRAC06.SiO2.v4c.nc',\
        '.TRAC07.FeT.v4c.nc',\
        '.TRAC08.DOC.v4c.nc',\
        '.TRAC09.DON.v4c.nc',\
        '.TRAC10.DOP.v4c.nc',\
        '.TRAC11.DOFe.v4c.nc',\
        '.TRAC12.POC.v4c.nc',\
        '.TRAC13.PON.v4c.nc',\
        '.TRAC14.POP.v4c.nc',\
        '.TRAC15.POSi.v4c.nc',\
        '.TRAC16.POFe.v4c.nc',\
        '.TRAC17.PIC.v4c.nc',\
        '.TRAC18.ALK.v4c.nc',\
        '.TRAC19.O2.v4c.nc',\
        '.TRAC20.Euk01.v4c.nc',\
        '.TRAC21.Syn01.v4c.nc',\
        '.TRAC22.Zoo01.v4c.nc'\
    ]
    files=files_dia
    
    for F in files:
        print(F,PD_NC)
        I_name=pre+F
        print(I_name)
        O_name=post+F
        print(O_name)
        ds_tracers=xr.open_mfdataset(I_name,  combine='by_coords', parallel=False,chunks={'T':1,'Z':1})
        ds_tracersZ=ds_tracers.rename({str('Zmd000023'):str('Z')})
        ds_tracersZ=ds_tracersZ.assign_coords(Z=depth_values)
        ds_tracersZ.to_netcdf(path=O_name)
        ds_tracers.close()
        ds_tracersZ.close()
    ########################
    #
    
    #
    print('EXIT  :: EXIT', flush=True, file = sys.stdout)
   
    
   
#
if __name__ == "__main__":
    #Initialize
    
    main()
