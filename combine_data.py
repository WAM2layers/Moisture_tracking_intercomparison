#### Necessary libraries ####
import numpy as np  # Numpy is the fundamental package for scientific computing in Python.
import netCDF4 as nc  # NetCDF is the data format of the meteorological data that we use.
import xarray as xr
import datetime  # Datetime is a package to deal with dates.
import os
from Functions import *


def read_data_pakistan(basedir):
    """Examples of data loading of moisture sources

    Data is loaded per individual ensemble member for each model.
    All data is converted to mm evaporative sources over the whole time period

    WRF-WVT is not included yet, as well CHc Lagranto and univie FLEXPART

    There might also be a few (new) ensemble members not included yet:
    - For the HAMSTER model (Ughent) there is an additional ensemble member not included yet
    - Extra ensemble members of Flexpart-Watersip produced by Fandy
    """

    ########################################################
    ## WRF-WVT                                            ##
    ########################################################

    # Just csv file because of nature of simulations #

    ########################################################
    ## WAM2layers                                         ##
    ########################################################
    directory_str = basedir + "results WAM2layers/"
    directory = os.fsencode(directory_str)
    n = 0

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".nc"):
            if n == 0:
                temp = xr.open_dataset(os.path.join(directory_str, filename))
                a_gridcell, lx, ly = get_grid_info(temp)
                srcs_wam2layers = temp["e_track"] * 1000 / a_gridcell[:, None]
                n += 1
            else:
                temp = xr.open_dataset(os.path.join(directory_str, filename))
                a_gridcell, lx, ly = get_grid_info(temp)
                srcs_wam2layers += temp["e_track"] * 1000 / a_gridcell[:, None]
                n += 1
            continue
        else:
            continue
    srcs_wam2layers = srcs_wam2layers.rename(latitude="lat", longitude="lon")

    # Data loading this way gives an error message for me while plotting, but in principe it should work
    dsall = xr.open_mfdataset(
        basedir + "results WAM2layers/backtrack_*T00-00.nc",
        combine="nested",
        concat_dim="time",
    ).rename(latitude="lat", longitude="lon")
    lat = dsall.lat.values
    lon = dsall.lon.values
    a_gridcell_new, l_ew_gridcell, l_mid_gridcell = get_grid_info_new(lat, lon)
    E_track_totalmm = (dsall / a_gridcell_new) * 1000  # mm
    srcs_wam2layers_new = E_track_totalmm["e_track"].sum("time")  # mm

    ########################################################
    ## University of Vigo                                 ##
    ########################################################
    srcs_Vigo_e1_Stohl = xr.open_dataset(basedir + "results Uvigo/ERA5_SJ05_reg.nc")[
        "E_P"
    ]
    srcs_Vigo_e2_Sodemann = xr.open_dataset(
        basedir + "results Uvigo/ERA5_APA22_reg.nc"
    )["E_P"]

    ########################################################
    ## UTRACK - Arie Staal                                ##
    ## - Based on the communication with Arie             ##
    ## - the results should be area corrected             ##
    ## - given the values I assumed that the area needs   ##
    ## - to be in km**2, however this should be checked.  ##
    ########################################################
    directory_str = (
        basedir + "results Utrack Arie Staal/moisture_tracking_intercomparsion/"
    )
    directory = os.fsencode(directory_str)

    dsall = xr.open_mfdataset(
        basedir
        + "results Utrack Arie Staal/moisture_tracking_intercomparsion/*_mixing48h_dt025h_100p.nc",
        combine="nested",
        concat_dim="time",
    )
    a_gridcell_new, l_ew_gridcell, l_mid_gridcell = get_grid_info_new(
        dsall.lat.values, dsall.lon.values
    )  # Calcluate grid cell area

    n_gridcells_pakistan = (71 - 67) / 0.25 * (30 - 24) / 0.25

    srcs_utrack_e1 = dsall["moisture_source"].sum("time") * n_gridcells_pakistan
    dsall = xr.open_mfdataset(
        basedir
        + "results Utrack Arie Staal/moisture_tracking_intercomparsion/*_mixing24h_dt025h_100p.nc",
        combine="nested",
        concat_dim="time",
    )
    srcs_utrack_e2 = dsall["moisture_source"].sum("time") * n_gridcells_pakistan
    dsall = xr.open_mfdataset(
        basedir
        + "results Utrack Arie Staal/moisture_tracking_intercomparsion/*_mixing12h_dt025h_100p.nc",
        combine="nested",
        concat_dim="time",
    )
    srcs_utrack_e3 = dsall["moisture_source"].sum("time") * n_gridcells_pakistan
    dsall = xr.open_mfdataset(
        basedir
        + "results Utrack Arie Staal/moisture_tracking_intercomparsion/*_mixing24h_dt025h_1000p.nc",
        combine="nested",
        concat_dim="time",
    )
    srcs_utrack_e4 = dsall["moisture_source"].sum("time") * n_gridcells_pakistan
    dsall = xr.open_mfdataset(
        basedir
        + "results Utrack Arie Staal/moisture_tracking_intercomparsion/*_mixing24h_dt010h_100p.nc",
        combine="nested",
        concat_dim="time",
    )
    srcs_utrack_e5 = dsall["moisture_source"].sum("time") * n_gridcells_pakistan

    ########################################################
    ## HAMSTER (Ghent)                                    ##
    ########################################################

    # E1: Sodemann #
    temp = xr.open_dataset(
        basedir
        + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens1_sod08/bias_corrected_20220810120000_sod08.nc"
    )
    srcs_ghent_e1 = temp["E2P_BC"].mean(
        "time"
    )  # Mean over time to remove time dimension

    for date in range(11, 25):
        temp = xr.open_dataset(
            basedir
            + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens1_sod08/bias_corrected_202208"
            + str(date)
            + "120000_sod08.nc"
        )
        srcs_ghent_e1 += temp["E2P_BC"].mean("time")

    # E2: FAS19 (Fremme + Sodemann, 2019) #
    temp = xr.open_dataset(
        basedir
        + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens2_fas19/bias_corrected_20220810120000_fas19.nc"
    )
    srcs_ghent_e2 = temp["E2P_BC"].mean("time")

    for date in range(11, 25):
        temp = xr.open_dataset(
            basedir
            + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens2_fas19/bias_corrected_202208"
            + str(date)
            + "120000_fas19.nc"
        )
        srcs_ghent_e2 += temp["E2P_BC"].mean("time")

    # E3: FAS19 (Fremme + Sodemann, 2019) #
    temp = xr.open_dataset(
        basedir
        + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens3_rh20/bias_corrected_20220810120000_rh20.nc"
    )
    srcs_ghent_e3 = temp["E2P_BC"].mean("time")

    for date in range(11, 25):
        temp = xr.open_dataset(
            basedir
            + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens3_rh20/bias_corrected_202208"
            + str(date)
            + "120000_rh20.nc"
        )
        srcs_ghent_e3 += temp["E2P_BC"].mean("time")

    # E4:
    temp = xr.open_dataset(
        basedir
        + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens4_allabl/bias_corrected_20220810120000_allabl.nc"
    )
    srcs_ghent_e4 = temp["E2P_BC"].mean("time")

    for date in range(11, 25):
        temp = xr.open_dataset(
            basedir
            + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens4_allabl/bias_corrected_202208"
            + str(date)
            + "120000_allabl.nc"
        )
        srcs_ghent_e4 += temp["E2P_BC"].mean("time")

    # E5:
    temp = xr.open_dataset(
        basedir
        + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens5_rhadap80/bias_corrected_20220810120000_rhadap80.nc"
    )
    srcs_ghent_e5 = temp["E2P_BC"].mean("time")

    for date in range(11, 25):
        temp = xr.open_dataset(
            basedir
            + "results UGhent HAMSTER/Pakistan simulations/bias_corrected/ens5_rhadap80/bias_corrected_202208"
            + str(date)
            + "120000_rhadap80.nc"
        )
        srcs_ghent_e5 += temp["E2P_BC"].mean("time")

    ########################################################
    ## TRACMASS Dipanjan Dey                              ##
    ## Units in mm/day, so multiplied with # of           ##
    ## event days                                         ##
    ########################################################
    nrdays = 15
    ds_TRACMASS = xr.open_dataset(
        basedir + "results TRACMASS Dipanjan Dey/TRACMASS_diagnostics.nc"
    )  # Evaporative sources (and preicp?) mm/day
    ds_pr_TRACMASS = xr.open_dataset(
        basedir + "results TRACMASS Dipanjan Dey/PR_ERA5_TRACMASS.nc"
    )  # Precip ERA5 and TRACMASS Comparison
    # convert to -180 to 180 lon
    ds_TRACMASS.coords["lon"] = (ds_TRACMASS.coords["lon"] + 180) % 360 - 180
    ds_TRACMASS = ds_TRACMASS.sortby(ds_TRACMASS.lon)

    srcs_TRACMASS = (
        ds_TRACMASS["E_TRACMASS"] * nrdays
    )  # Units of data is mm/day but we want mm over whole time period

    ########################################################
    ## FLEXPART-Watersip TatFanCheng                      ##
    ########################################################
    filename = (
        basedir
        + "results FLEXPART_WaterSip_TatFanCheng/WaterSip_Cb_20220810-20220824_Pakistan_box.nc"
    )
    ds_flexpart_tatfancheng = xr.open_dataset(filename)
    # convert to -180 to 180 lon
    ds_flexpart_tatfancheng.coords["lon"] = (
        ds_flexpart_tatfancheng.coords["lon"] + 180
    ) % 360 - 180
    ds_flexpart_tatfancheng = ds_flexpart_tatfancheng.sortby(ds_TRACMASS.lon)
    srcs_flexpart_tatfancheng = ds_flexpart_tatfancheng.sum("time")["Cb"]

    ########################################################
    ## Flexpart Ru Xu                                     ##
    ########################################################
    ds_flexpart_xu = xr.open_dataset(
        basedir + "results Ru_Xu_FLEXPART/e_daily.nc"
    ).rename(latitude="lat", longitude="lon")
    srcs_flexpart_xu = ds_flexpart_xu["variable"].sum("time")

    ########################################################
    ## Lagranto CHc                                       ##
    ########################################################
    ds_lagranto_CHc = xr.open_dataset(
        basedir + "results CHc LAGRANTO/Pakistan_2022_CHc_eventtotal_ens1.nc"
    ).rename(dimx_N="lon", dimy_N="lat")
    srcs_lagranto_CHc = ds_lagranto_CHc["N"].sum("time").squeeze()
    srcs_lagranto_CHc = srcs_lagranto_CHc.assign_coords(
        lat=srcs_TRACMASS.lat[::-1], lon=srcs_TRACMASS.lon
    )

    ########################################################
    ## FLEXPART UniVie                                    ##
    ########################################################
    ds_flexpart_univie = xr.open_dataset(
        basedir + "results univie FLEXPART/pakistan_univie.nc"
    )
    srcs_flexpart_univie = (
        ds_flexpart_univie["moisture_uptakes_bl"]
        + ds_flexpart_univie["moisture_uptakes_ft"]
    ).sum("time")

    ########################################################
    ## 2LDRM                                              ##
    ########################################################
    ds_2ldrm = xr.open_dataset(
        basedir + "results 2LDRM/2LDRM_Pakistan_case_gl.nc"
    ).rename(latitude="lat", longitude="lon")
    srcs_2ldrm = ds_2ldrm["moisture_source"].sum("time").T

    ########################################################
    ## FLEXPART UiB                                       ##
    ########################################################
    ds_flexpart_uib = xr.open_dataset(
        basedir
        + "results UiB FLEXPART WaterSip/Pakistan_2022_UiB_Sodemann_grid_EN1_regridded.nc"
    )
    srcs_flexpart_uib = ds_flexpart_uib["moisture_uptakes"]

    # Combine cumulative moisture sources for all models in one netcdf
    datasets = {
        "2ldrm": srcs_2ldrm,
        "flexpart_uib": srcs_flexpart_uib,
        "flexpart_univie": srcs_flexpart_univie,
        "lagranto_CHc": srcs_lagranto_CHc,
        "flexpart_xu": srcs_flexpart_xu,
        "flexpart_tatfancheng": srcs_flexpart_tatfancheng,
        "TRACMASS": srcs_TRACMASS,
        "ghent_e5": srcs_ghent_e5,
        "ghent_e4": srcs_ghent_e4,
        "ghent_e3": srcs_ghent_e3,
        "ghent_e2": srcs_ghent_e2,
        "ghent_e1": srcs_ghent_e1,
        "utrack_e5": srcs_utrack_e5,
        "utrack_e4": srcs_utrack_e4,
        "utrack_e3": srcs_utrack_e3,
        "utrack_e2": srcs_utrack_e2,
        "utrack_e1": srcs_utrack_e1,
        "Vigo_e2_Sodemann": srcs_Vigo_e2_Sodemann,
        "Vigo_e1_Stohl": srcs_Vigo_e1_Stohl,
        "wam2layers_new": srcs_wam2layers_new,
    }

    all_maps = xr.Dataset(datasets)

    return all_maps