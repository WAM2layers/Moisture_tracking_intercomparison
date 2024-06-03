#### Necessary libraries ####
import numpy as np              # Numpy is the fundamental package for scientific computing in Python.
import netCDF4 as nc            # NetCDF is the data format of the meteorological data that we use.
import xarray as xr
import matplotlib.pyplot as plt # Matplotlib is a scientific plotting package.
import datetime                 # Datetime is a package to deal with dates.
import cartopy.crs as crs
import cartopy
from cmocean import cm
import os
#import regionmask

from Functions import *

from combine_data import read_data

ds_scotland = read_data('./', 'Scotland')

ens_names = ds_scotland.keys()


mask=xr.open_dataset("mask_scotland_box_lon-180to180.nc")




print("--> Plotting sources")
### MAKE PLOT OF DIFFERENT ENSEMBLE MEMBERS (not corrected for different sums) ###

plotting_sources(ds_scotland, mask,ens_names, figwidth=24, figheight=14, vmax=5, central_longitude=0, figrows=5, figcols=5, map_lons_extend=[-85, 40], map_lats_extend=[10,80], fname="Scotland_case_v3.png")

print("--> PLotting fraction")
plotting_frac(ds_scotland, mask, ens_names,  figwidth=24, figheight=14, vmax=0.02, central_longitude=0, figrows=5, figcols=5, map_lons_extend=[-85, 40], map_lats_extend=[10,80], fname="Scotland_fracV3.png")
