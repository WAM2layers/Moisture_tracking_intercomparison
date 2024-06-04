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

basedir = "./"
ds_scotland = read_data(basedir, "Scotland")

ens_names = ds_scotland.keys()


mask=xr.open_dataset(basedir+"/Scotland/mask_scotland_box_lon-180to180.nc")


print("--> Plotting single method")
method = "WAM2layers"
fig, ax = plt.subplots(1, 1, figsize=(15,10),subplot_kw={'projection': crs.PlateCarree()},sharey=False)
plot_single_map(ax, ds_scotland[method], mask, method, [-85, 40], [10, 80], vmax=5, xlabel="Longitude", ylabel="Latitude")
fig.savefig(f"Scotland_{method}_v3.png")

print("--> Plotting sources")
### MAKE PLOT OF DIFFERENT ENSEMBLE MEMBERS (not corrected for different sums) ###

plotting_sources(ds_scotland, mask,ens_names, "Scotland", figwidth=24, figheight=14, vmax=5, central_longitude=0, figrows=5, figcols=5, fname="Scotland_case_v3.png")

print("--> PLotting fraction")
ds_scotland_frac = xr.Dataset(coords=ds_scotland.coords, attrs=ds_scotland.attrs)
for ens in ds_scotland.keys():
    ds_scotland_frac[ens] = calc_fractional_sources(ds_scotland[ens])
plotting_sources(ds_scotland_frac, mask, ds_scotland_frac.keys(), "Scotland", figwidth=24, figheight=14, vmax=0.02, central_longitude=0, figrows=5, figcols=5, fname="Scotland_fracV3.png")