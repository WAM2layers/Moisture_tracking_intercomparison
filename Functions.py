import numpy as np


def grid_cell_area(latitude, longitude):
    """Return grid cell area"""
    erad = 6.371e6  # [m] Earth radius

    dx = np.abs(longitude[1] - longitude[0])  # [degrees] grid cell size

    # new area size calculation:
    lat_n_bound = np.minimum(90.0, latitude + 0.5 * dx)
    lat_s_bound = np.maximum(-90.0, latitude - 0.5 * dx)

    a_gridcell = np.zeros([len(latitude), 1])
    a_gridcell[:, 0] = (
        (np.pi / 180.0)
        * erad**2
        * abs(np.sin(lat_s_bound * np.pi / 180.0) - np.sin(lat_n_bound * np.pi / 180.0))
        * dx
    )

    return a_gridcell


def calc_fractional_sources(
    sources, precipitation=None, lon_name="lon", lat_name="lat"
):
    """Returns fractional moisture sources

    Standard assumes that the names of longitude and latitude
    dimension are "lon" and "lat" respectively. If this is not the case
    lon_name and lat_name should be used.

    Input:
    sources: (xarray) array containing the moisture sources
    lon_name: Name of longitude dimension
    lat_name: Name of latitude dimension
    precipitation: Optional array containing the precipitation of the model

    If no precipitation array is provided, thre fractional sources will be caluclated
    with respect to the total moisture sources

    Otherwise it will calculate fractional sources with respect to the total precipitation
    of the event*
    *DOES NOT WORK YET!"""

    if precipitation == None:
        a_gridcell_new, l_ew_gridcell, l_mid_gridcell = get_grid_info_new(
            sources[lat_name], sources[lon_name]
        )
        sources_frac = (
            sources
            * a_gridcell_new
            / ((sources * a_gridcell_new).sum((lon_name, lat_name)))
            * 100.0
        )
    else:
        """ Should be adapted so that the area-weighted precipitation can be used to bias correct
        moisture sources """

    return sources_frac


def calc_regional_sources(
    sources, regions, weights=None, lon_name="lon", lat_name="lat"
):
    """Returns regional moisture sources

    Standard assumes that the names of longitude and latitude
    dimension are "lon" and "lat" respectively. If this is not the case
    lon_name and lat_name should be used.

    Weights can be re-used for different ensemble memebers per model"""

    if weights is not None:
        sources_per_region = sources.weighted(weights).sum(dim=(lat_name, lon_name))
    else:
        mask_3D = regions.mask_3D(sources, lon_name=lon_name, lat_name=lat_name)
        weights = np.cos(np.deg2rad(sources[lat_name]))
        weights = (mask_3D * weights).fillna(0)
        sources_per_region = sources.weighted(weights).sum(dim=(lat_name, lon_name))

    return sources_per_region, weights
