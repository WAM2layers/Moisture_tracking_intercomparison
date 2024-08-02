import numpy as np
import cartopy.crs as crs
import cartopy
from cmocean import cm
import matplotlib.pyplot as plt

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
        # Use equal weights for each region as we are already using area-weighted fractional sources.
        # Effectively, the weighting is only used to group grid points of the same region together.
        weights = mask_3D.fillna(0)
        sources_per_region = sources.weighted(weights).sum(dim=(lat_name, lon_name))

    return sources_per_region, weights

# Get grid info
def get_grid_info_new(latitude, longitude):
    """Return grid cell area and lenght sides."""
    dg = 111089.56  # [m] length of 1 degree latitude
    erad = 6.371e6  # [m] Earth radius

    gridcell = np.abs(longitude[1] - longitude[0])  # [degrees] grid cell size

    # new area size calculation:
    lat_n_bound = np.minimum(90.0, latitude + 0.5 * gridcell)
    lat_s_bound = np.maximum(-90.0, latitude - 0.5 * gridcell)

    a_gridcell = np.zeros([len(latitude), 1])
    a_gridcell[:, 0] = (
        (np.pi / 180.0)
        * erad ** 2
        * abs(np.sin(lat_s_bound * np.pi / 180.0) - np.sin(lat_n_bound * np.pi / 180.0))
        * gridcell
    )

    l_ew_gridcell = gridcell * dg  # [m] length eastern/western boundary of a cell
    l_n_gridcell = (
        dg * gridcell * np.cos((latitude + gridcell / 2) * np.pi / 180)
    )  # [m] length northern boundary of a cell
    l_s_gridcell = (
        dg * gridcell * np.cos((latitude - gridcell / 2) * np.pi / 180)
    )  # [m] length southern boundary of a cell
    l_mid_gridcell = 0.5 * (l_n_gridcell + l_s_gridcell)
    return a_gridcell, l_ew_gridcell, l_mid_gridcell

def plot_single_map(ax, data, mask, title, map_lons_extend, map_lats_extend, vmax=5, xlabel="Longitude", ylabel="Latitude"):

    data.plot(ax=ax,vmin=0,vmax=vmax,robust=False,cmap=cm.rain,
              cbar_kwargs=dict(fraction=0.05, shrink=0.5,label=None),)
    ax.set_title(title, loc="left")
    ax.contour(mask['lon'].values, mask['lat'].values, np.squeeze(mask['mask'].values),colors=["r"])
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=.2)

    ax.set_xticks(np.arange(-180, 181, 20))
    ax.set_yticks(np.arange(-90, 91, 10))
    ax.set_xlim(map_lons_extend[0], map_lons_extend[1])
    ax.set_ylim(map_lats_extend[0], map_lats_extend[1])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def plotting_sources(ds_data, mask, ens_names, casename, figwidth=24, figheight=14, vmax=5, central_longitude=0, figrows=5, figcols=5, map_lons_extend=None, map_lats_extend=None, fname="fig"):

    # set map extent if it is not set
    if map_lons_extend is None:
        if casename == "Pakistan":
            map_lons_extend = [0, 140]
        elif casename == "Scotland":
            map_lons_extend = [-85, 40]
        elif casename == "Australia":
            map_lons_extend = [120, 240]

    if map_lats_extend is None:
        if casename == "Pakistan":
            map_lats_extend = [-40,50]
        elif casename == "Scotland":
            map_lats_extend = [10, 80]
        elif casename == "Australia":
            map_lats_extend = [-65, 10]

    my_projection = crs.PlateCarree(central_longitude=central_longitude)

    # Make figure
    fig, axs = plt.subplots(figrows, figcols, figsize=(figwidth, figheight),subplot_kw={'projection': crs.PlateCarree()},sharey=False)

    for n,ens in enumerate(ens_names):

        print("------  Plotting", ens)

        i = n//figcols
        j = n%figcols

        # Dismiss label of y-axis, except for left most column
        if(j > 0):
            ylabel=""
        else:
            ylabel="Latitude"

        # Dismiss label of x-axis except for bottom row
        if(i < figrows-1):
            xlabel=""
        else:
            xlabel="Longitude"

        plot_single_map(axs[i,j], ds_data[ens], mask, ens, map_lons_extend, map_lats_extend, vmax=vmax, xlabel=xlabel, ylabel=ylabel)

    # Remove axes from empty subplots
    for n in range(len(ens_names), figrows*figcols):
        i = n//figcols
        j = n%figcols
        plt.sca(axs[i,j])
        plt.axis("off")

    fig.savefig(fname, dpi=300,  bbox_inches="tight")
