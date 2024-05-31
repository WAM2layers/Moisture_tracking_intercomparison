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


ens_names={
        "wam2layers":"WAM2layers",
        "Vigo_e1_Stohl":"UVigo E1",
        "Vigo_e2_Sodemann": "UVigo E2",
        "utrack_e1": "UTRACK E1",
        "utrack_e2": "UTRACK E2",
        "utrack_e3": "UTRACK E3",
        "utrack_e4": "UTRACK E4",
        "utrack_e5": "UTRACK E5",
        "ghent_e1":"HASMTER E1",
        "ghent_e2":"HASMTER E2",
        "ghent_e3":"HASMTER E3",
        "ghent_e4":"HASMTER E4",
        "ghent_e5":"HASMTER E5",
        "TRACMASS":"TRACMASS",
        "flexpart_tatfancheng_Ens1":"FLEXPART-Watersip E1",
        "flexpart_tatfancheng_Ens2":"FLEXPART-Watersip E2",
        "flexpart_tatfancheng_Ens3":"FLEXPART-Watersip E3",
        "flexpart_xu":"FLEXPART-Xu",
        "lagranto_CHc":"Hc LAGRANTO",
        "flexpart_univie":"FLEXPART_UViena",
        "2ldrm":"2LDRM",
}



ens_outputs=["wam2layers", "Vigo_e1_Stohl", "Vigo_e2_Sodemann", "utrack_e1", "utrack_e2", "utrack_e3", "utrack_e4", "utrack_e5", "ghent_e2", "ghent_e2", "ghent_e3", "ghent_e4", "ghent_e5","TRACMASS", "flexpart_tatfancheng_Ens1","flexpart_tatfancheng_Ens2","flexpart_tatfancheng_Ens3", "flexpart_xu", "lagranto_CHc", "flexpart_univie","2ldrm"]




ds_scotland = read_data('./', 'Scotland')


mask=xr.open_dataset("mask_scotland_box_lon-180to180.nc")


for ens in ens_outputs:
       vars()[ens+"_frac"] = calc_fractional_sources(ds_scotland[ens])

print("--> Plotting sources")
### MAKE PLOT OF DIFFERENT ENSEMBLE MEMBERS (not corrected for different sums) ###

my_projection = crs.PlateCarree(central_longitude=0)
vmax=5
rows=5
cols=5
# Make figure
fig, axs = plt.subplots(rows, cols, figsize=(24, 14),subplot_kw={'projection': crs.PlateCarree()},sharey=False)
i=0
j=0
for ens in ens_outputs:


        print("------  Plotting", ens)

        ds_scotland[ens].plot(ax=axs[i,j],vmin=0,vmax=vmax,robust=False,cmap=cm.rain,
                     cbar_kwargs=dict(fraction=0.05, shrink=0.5,label=None),)
        axs[i,j].set_title(ens_names[ens], loc="left")
        axs[i,j].contour(mask['lon'].values, mask['lat'].values, mask['mask'].values[0,:],colors=["r"])
        axs[i,j].add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        axs[i,j].add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        axs[i,j].set_xticks(np.arange(-180, 181, 20), crs=my_projection)
        axs[i,j].set_yticks(np.arange(-90, 91, 10), crs=my_projection)
        axs[i,j].set_xlim(-85, 40)
        axs[i,j].set_ylim(10, 80)

        #Dismiss label of y-axis, except for left most column

        if(j > 0):
            axs[i,j].set_ylabel("")
        else:
            axs[i,j].set_ylabel("Latitude")


        if i==rows-1:
            axs[i,j].set_xlabel("Longitude")
        else:
            axs[i,j].set_xlabel("")

        if j<cols-1:
                i=i
                j=j+1
        else:
                i=i+1
                j=0


fig.savefig("Scotland_case_v2.png",dpi=300,  bbox_inches="tight")



print("--> PLotting fraction")
#General plotting options
vmin = 0
vmax = 0.02
contour_levels = [0.001, 0.01]

### MAKE PLOT OF DIFFERENT ENSEMBLE MEMBERS (Fractional) ###
my_projection = crs.PlateCarree(central_longitude=0)

# Make figure
rows=5
cols=5
fig, axs = plt.subplots(rows, cols, figsize=(24, 14),subplot_kw={'projection': crs.PlateCarree()},sharey=False)


i=0
j=0
for ens in ens_outputs:


        print("------  Plotting frac", ens)

        vars()[ens+"_frac"].plot(ax=axs[i,j],vmin=0,vmax=vmax,robust=False,cmap=cm.rain,
                        cbar_kwargs=dict(fraction=0.05, shrink=0.5,label=None),)
        #srcs_wam2layers.plot.contour(ax=axs[0,0], levels=[0.1, 1], colors=["lightgrey", "grey"])
        axs[i,j].set_title(ens_names[ens], loc="left")
        axs[i,j].contour(mask['lon'].values, mask['lat'].values, mask['mask'].values[0,:], colors=["r"])
        axs[i,j].add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        axs[i,j].add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

        axs[i,j].set_xticks(np.arange(-180, 181, 20), crs=my_projection)
        axs[i,j].set_yticks(np.arange(-90, 91, 10), crs=my_projection)
        axs[i,j].set_xlim(-85, 40)
        axs[i,j].set_ylim(10, 80)

        #Dismiss label of y-axis, except for left most column
        if(j > 0):
            axs[i,j].set_ylabel("")
        else:
            axs[i,j].set_ylabel("Latitude")


        if i==rows-1:
            axs[i,j].set_xlabel("Longitude")
        else:
            axs[i,j].set_xlabel("")

        if j<cols-1:
                i=i
                j=j+1
        else:
                i=i+1
                j=0

fig.savefig("Scotland_frac.png",dpi=300, bbox_inches="tight")





