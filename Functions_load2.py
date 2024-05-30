import numpy as np
import netCDF4 as nc            # NetCDF is the data format of the meteorological data that we use.
import xarray as xr
import datetime                 # Datetime is a package to deal with dates.
import os
import glob

from Functions import * 

### Model options
'''
Appled only to Australia Cases
'''

def load_model(model_name):
    if model_name ==  "univie FLEXPART":
        print('using boundary layer only')
        ds = xr.open_dataset('./results univie FLEXPART/australia_univie.nc')
        srcs = ds.moisture_uptakes_bl.sum(dim='time')
        return srcs
        
    elif model_name == "CHc LAGRANTO":
        ds = xr.open_dataset('./results CHc LAGRANTO/Australia_2022_CHc_eventtotal_ens1.nc')
        
        srcs = xr.Dataset(data_vars = {'moisture_source':(["latitude","longitude"], ds.isel(time=0, dimz_N=0).N.values)},
                 coords = {
                     "latitude": (["latitude"], np.linspace(-90,90,721)),
                     "longitude": (["longitude"], np.linspace(-180,180,1441)),
                 })
        srcs = srcs.moisture_source.rename(latitude='lat', longitude='lon')
        
        return srcs
            
    elif model_name ==  "Ru_Xu_FLEXPART":
        ds = xr.open_dataset('./results Ru_Xu_FLEXPART/aus_e_daily.nc')
        srcs = ds.sum(dim='time')
        return srcs
            
    elif model_name ==  "FLEXPART_WaterSip_TatFanCheng":
        directory_str = "results FLEXPART_WaterSip_TatFanCheng/WaterSip_moisture_source_Australia_"
            
        # Ens1
        ds = xr.open_dataset(directory_str + "20220222-20220228_Ens1.nc")
        #convert to -180 to 180 lon
        ds['lon'] = (ds['lon'] + 180) % 360 - 180
        ds = ds.sortby(ds.lon)
        srcs1 = ds["Cb"]
            
        # Ens2
        ds = xr.open_dataset(directory_str + "20220222-20220228_Ens2.nc")
        #convert to -180 to 180 lon
        ds['lon'] = (ds['lon'] + 180) % 360 - 180
        ds = ds.sortby(ds.lon)
        srcs2 = ds["Cb"]
            
        # Ens3
        ds = xr.open_dataset(directory_str + "20220222-20220228_Ens3.nc")
        #convert to -180 to 180 lon
        ds['lon'] = (ds['lon'] + 180) % 360 - 180
        ds = ds.sortby(ds.lon)
        srcs3 = ds["Cb"]
            
        return srcs1, srcs2, srcs3
            
    elif model_name ==  "WAM2layers":
        directory_str = "results WAM2layers/backtrack_*.nc"
        srcs = xr.open_mfdataset(directory_str).e_track.sum(dim='time').load()
        srcs = srcs.rename(latitude='lat', longitude='lon')
        return srcs
            
    elif model_name ==  "Utrack Arie Staal":
            directory_str = "results Utrack Arie Staal/utrack_back_australia_case_2022-02-22--28_mixing"
            directory =  os.fsencode(directory_str)

            # 48h_dt025h_100p
            ds = xr.open_mfdataset(directory_str + '48h_dt025h_100p.nc',combine = 'nested', concat_dim='time')
            a_gridcell_new, l_ew_gridcell, l_mid_gridcell = get_grid_info_new(ds.lat.values, ds.lon.values) #Calcluate grid cell area
            srcs1 = ds["moisture_source"].sum("time")*a_gridcell_new/10**6.0 #1000.0
            
            # 24h_dt025h_100p
            ds = xr.open_mfdataset(directory_str + '24h_dt025h_100p.nc',combine = 'nested', concat_dim='time')
            srcs2 = ds["moisture_source"].sum("time")*a_gridcell_new/10**6.0 #1000.0
            
            # 12h_dt025h_100p
            ds = xr.open_mfdataset(directory_str + '12h_dt025h_100p.nc',combine = 'nested', concat_dim='time')
            srcs3 = ds["moisture_source"].sum("time")*a_gridcell_new/10**6.0 #1000.0
            
            # 24h_dt025h_1000p
            ds = xr.open_mfdataset(directory_str + '24h_dt025h_1000p.nc',combine = 'nested', concat_dim='time')
            srcs4 = ds["moisture_source"].sum("time")*a_gridcell_new/10**6.0 #1000.0
            
            # 24h_dt010h_100p.nc
            ds = xr.open_mfdataset(directory_str + '24h_dt010h_100p.nc',combine = 'nested', concat_dim='time')
            srcs5 = ds["moisture_source"].sum("time")*a_gridcell_new/10**6.0 #1000.0

            return srcs1, srcs2, srcs3, srcs4, srcs5
            
    elif model_name ==  "TRACMASS Dipanjan Dey":            
            print('Does not have spatial information yet')
            
    elif model_name ==  "Uvigo":
            srcs1 = xr.open_dataset("results Uvigo/ERA5_APA22_reg.nc")["E_P"]
            srcs2 = xr.open_dataset("results Uvigo/ERA5_SJ05_reg.nc")["E_P"]
            return srcs1, srcs2
            
    elif model_name ==  "B-TrIMS":
            for day in range(22,29):
                ds = xr.open_dataset('./results_B-TrIMS/bt.202202_%02d_processed.nc' % day)
                if day == 22:
                    ods = ds.wvcont_mm_daily
                else:
                    ods += ds.wvcont_mm_daily
            srcs = xr.Dataset(data_vars = {'wvcont_mm_daily':(["latitude","longitude"], ods.values)},
                 coords = {
                     "latitude": (["latitude"], ds.latitude[:,0].values),
                     "longitude": (["longitude"], ds.longitude[0,:].values),
                 })
            srcs.coords['longitude'] = (srcs.coords['longitude'] + 180) % 360 - 180
            srcs = srcs.sortby(srcs.longitude)
            srcs = (srcs.wvcont_mm_daily).rename(latitude='lat', longitude='lon')
            return srcs
            
    elif model_name ==  "UGhent HAMSTER":      
            print('5 ensemble members')
            
            directory_str = "results UGhent HAMSTER/Australia simulations/bias_corrected/"
            # E1: Sodemann #
            temp = xr.open_mfdataset(directory_str + "/ens1_sod08/*.nc")
            srcs1 = temp["E2P_BC"].sum("time")  # sum over time
            
            # E2: FAS19 (Fremme + Sodemann, 2019) #
            temp = xr.open_mfdataset(directory_str + "/ens2_fas19/*.nc")
            srcs2 = temp["E2P_BC"].sum("time")  # sum over time
            
            # E3: RH20 #
            temp = xr.open_mfdataset(directory_str + "/ens3_rh20/*.nc")
            srcs3 = temp["E2P_BC"].sum("time")  # sum over time
            
            # E4: ALLABL #
            temp = xr.open_mfdataset(directory_str + "/ens4_allabl/*.nc")
            srcs4 = temp["E2P_BC"].sum("time")  # sum over time
            
            # E5: RHADAP80 #
            temp = xr.open_mfdataset(directory_str + "/ens5_rhadap80/*.nc")
            srcs5 = temp["E2P_BC"].sum("time")  # sum over time
            
            return srcs1, srcs2, srcs3, srcs4, srcs5
        
    elif model_name ==  "2LDRM":  
        directory_str = "results 2LDRM/2LDRM_Australia_case_gl.nc"
        srcs = xr.open_dataset(directory_str).sum(dim='time').moisture_source
        srcs.coords['longitude'] = (srcs.coords['longitude'] + 180) % 360 - 180
        srcs = srcs.sortby(srcs.longitude)
        srcs = (srcs.T).rename(latitude='lat', longitude='lon')
        return srcs
