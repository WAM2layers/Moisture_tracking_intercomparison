import numpy as np

## WAM2layers functions needed to load data and convert units

#### Definitions regions/etc. #### 
def load_tagging_region(config, t=None):
    tagging_region = xr.open_dataarray(config.tagging_region)

    if config.tracking_domain is not None:
        tagging_region = select_subdomain(tagging_region, config.tracking_domain)

    if t is not None:
        return tagging_region.sel(time=t, method="nearest")
    return tagging_region

#region =  tagging_region.sel(time=t, method="nearest")
def get_grid_info(ds):
    """Return grid cell area and lenght of the sides."""
    dg = 111089.56  # [m] length of 1 degree latitude
    erad = 6.371e6  # [m] Earth radius

    latitude = ds.latitude.values
    longitude = ds.longitude.values
    grid_spacing = np.abs(longitude[1] - longitude[0])  # [degrees]

    # Calculate area TODO check this calculation!
    lat_n = np.minimum(90.0, latitude + 0.5 * grid_spacing)
    lat_s = np.maximum(-90.0, latitude - 0.5 * grid_spacing)

    a = (
        np.pi
        / 180.0
        * erad**2
        * grid_spacing
        * abs(np.sin(lat_s * np.pi / 180.0) - np.sin(lat_n * np.pi / 180.0))
    )

    # Calculate faces
    ly = grid_spacing * dg  # [m] length eastern/western boundary of a cell
    lx_n_gridcell = ly * np.cos((latitude + grid_spacing / 2) * np.pi / 180)
    lx_s_gridcell = ly * np.cos((latitude - grid_spacing / 2) * np.pi / 180)
    lx = 0.5 * (lx_n_gridcell + lx_s_gridcell)
    return a, ly, lx


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

def calc_fractional_sources(sources, precipitation=None, lon_name="lon",lat_name="lat"):
    """ Returns fractional moisture sources 
    
        Standard assumes that the names of longitude and latitude 
        dimension are "lon" and "lat" respectively. If this is not the case
        lon_name and lat_name should be used."""
    
    if(precipitation == None):
        a_gridcell_new, l_ew_gridcell, l_mid_gridcell = get_grid_info_new(sources[lat_name],
                                                                         sources[lon_name])
        sources_frac = sources*a_gridcell_new/((sources*a_gridcell_new).sum((lon_name,lat_name)))*100.0
    else:
        """ Should be adapted so that the area-weighted precipitation can be used to bias correct
        moisture sources """
        
    return sources_frac