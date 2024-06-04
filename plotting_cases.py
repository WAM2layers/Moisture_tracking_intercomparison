#### Necessary libraries ####
import numpy as np              # Numpy is the fundamental package for scientific computing in Python.
import netCDF4 as nc            # NetCDF is the data format of the meteorological data that we use.
import xarray as xr
import matplotlib.pyplot as plt # Matplotlib is a scientific plotting package.
import datetime                 # Datetime is a package to deal with dates.
import cartopy.crs as crs
import cartopy
import matplotlib.ticker as mticker
import os
#import regionmask

from Functions import *

from combine_data import read_data

masks={
	"Pakistan":"Mask_PakistanFlood_box_lon-180to180.nc",
	"Scotland":"mask_scotland_box_lon-180to180.nc",
	"Australia":"mask_australia_box_lon0to360.nc",
	}


maps_features={
	#case: [central_longitude, lonmin, lonmax, latmin, latmax]
	"Scotland":[0, -85, 40, 10, 80],
	"Pakistan":[0, 0, 140, -40, 50],
	'Australia':[180, 120, 240, -60, 10]
	}

fig_features={
	#case:[figwidth, figheight, figrows, figcols , vmax, vmax_frac]
	"Scotland": [24,14,5,5,5, 0.02],
	"Pakistan": [24,14,4,5,20, 0.02],
	"Australia": [25,14,4,5,25, 0.02],
	}


fig_features_one={
	#Figure details to plot individual cases
	#case:[figwidth, figheight, figrows, figcols , vmax, vmax_frac]
	"Scotland": [12,10,1,1,5, 0.02],
	"Pakistan": [12,10,1,1,20, 0.02],
	"Australia": [12,10,1,1,25, 0.02],
	}


cases=['Pakistan', 'Australia',"Scotland"]


means={}
stds={}
for case in cases:


	ds_data = read_data('./', case)

	mean = ds_data.to_array(dim='mean').mean('mean')
	std= ds_data.to_array(dim='std').std('std')

	means[case]= mean
	stds[case] = std

	ens_names = ds_data.keys()

	mask=xr.open_dataset(case+"/"+masks[case])

	print("--> Plotting sources")
	plotting_sources_cases(ds_data, mask, ds_data.keys(), figwidth=fig_features[case][0], figheight=fig_features[case][1], vmax=fig_features[case][4], central_longitude=maps_features[case][0], figrows=fig_features[case][2], figcols=fig_features[case][3], map_lons_extend=[maps_features[case][1], maps_features[case][2]], map_lats_extend=[maps_features[case][3],maps_features[case][4]], fname="Plot_"+case+"_case.png")

	#quit()

	print("--> PLotting fraction")
	ds_data_frac = xr.Dataset(coords=ds_data.coords, attrs=ds_data.attrs)
	for ens in ds_data.keys():
		ds_data_frac[ens] = calc_fractional_sources(ds_data[ens])

	plotting_sources_cases(ds_data_frac, mask, ds_data.keys(), figwidth=fig_features[case][0], figheight=fig_features[case][1], vmax=fig_features[case][5], central_longitude=maps_features[case][0], figrows=fig_features[case][2], figcols=fig_features[case][3], map_lons_extend=[maps_features[case][1], maps_features[case][2]], map_lats_extend=[maps_features[case][3],maps_features[case][4]], fname="Plot_frac_"+case+".png")


	#to plot one case: UNCOMMENT THESE LINES FOR PLOTTING AN SPECIFIC MODEL
	print("\n--> Plotting one cases")

	ens_names_case=['WAM2layers']  #name of the model. EX: ens_names_case=['WAM2layers'].

	plotting_sources_one_case(ds_data, mask, ens_names_case, figwidth=fig_features_one[case][0], figheight=fig_features_one[case][1], vmax=fig_features_one[case][4], central_longitude=maps_features[case][0], figrows=fig_features_one[case][2], figcols=fig_features_one[case][3], map_lons_extend=[maps_features[case][1], maps_features[case][2]], map_lats_extend=[maps_features[case][3],maps_features[case][4]], fname="Plot_"+case)


	plotting_sources_one_case(ds_data_frac, mask, ens_names_case,  figwidth=fig_features_one[case][0], figheight=fig_features_one[case][1], vmax=fig_features_one[case][5], central_longitude=maps_features[case][0], figrows=fig_features_one[case][2], figcols=fig_features_one[case][3], map_lons_extend=[maps_features[case][1], maps_features[case][2]], map_lats_extend=[maps_features[case][3],maps_features[case][4]], fname="Plot_frac_"+case)

	print()

print("\nPlotting_mean_global_map")
#proj=1: PlateCarree projection
#proj=2: Robinson projection
proj=2
combined_means = xr.Dataset(means)
plotting_global(combined_means, ["Scotland", "Pakistan", "Australia"], [cm.matter, cm.rain, cm.dense], vmaxs= [fig_features["Scotland"][4], fig_features["Pakistan"][4], fig_features["Australia"][4]], masks=masks, plot_vline=True, proj = proj, fname="global_mean_map.png")


print("\nPlotting_std_global_map")
combined_stds = xr.Dataset(stds)
plotting_global(combined_stds, ["Scotland", "Pakistan", "Australia"], [plt.cm.RdPu, plt.cm.YlGn, plt.cm.Blues], vmaxs=[10, 10, 10], masks=masks, plot_vline=True, proj = proj, fname="global_std_map.png")

