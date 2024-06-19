import numpy as np
import cartopy.crs as crs
import cartopy
from cmocean import cm
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter
import matplotlib.ticker as mticker
import xarray as xr


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
        * erad**2
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


def plot_single_map(
    ax,
    data,
    mask,
    title,
    map_lons_extend,
    map_lats_extend,
    vmax=5,
    xlabel="Longitude",
    ylabel="Latitude",
):
    data.plot(
        ax=ax,
        vmin=0,
        vmax=vmax,
        robust=False,
        cmap=cm.rain,
        cbar_kwargs=dict(fraction=0.05, shrink=0.5, label=None),
    )
    ax.set_title(title, loc="left")
    ax.contour(
        mask["lon"].values,
        mask["lat"].values,
        np.squeeze(mask["mask"].values),
        colors=["r"],
    )
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, linestyle="-", linewidth=0.2)

    ax.set_xticks(np.arange(-180, 181, 20))
    ax.set_yticks(np.arange(-90, 91, 10))
    ax.set_xlim(map_lons_extend[0], map_lons_extend[1])
    ax.set_ylim(map_lats_extend[0], map_lats_extend[1])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plotting_sources(
    ds_data,
    mask,
    ens_names,
    casename,
    figwidth=24,
    figheight=14,
    vmax=5,
    central_longitude=0,
    figrows=5,
    figcols=5,
    map_lons_extend=None,
    map_lats_extend=None,
    fname="fig",
):
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
            map_lats_extend = [-40, 50]
        elif casename == "Scotland":
            map_lats_extend = [10, 80]
        elif casename == "Australia":
            map_lats_extend = [-65, 10]

    my_projection = crs.PlateCarree(central_longitude=central_longitude)

    # Make figure
    fig, axs = plt.subplots(
        figrows,
        figcols,
        figsize=(figwidth, figheight),
        subplot_kw={"projection": crs.PlateCarree()},
        sharey=False,
    )

    for n, ens in enumerate(ens_names):
        print("------  Plotting", ens)

        i = n // figcols
        j = n % figcols

        # Dismiss label of y-axis, except for left most column
        if j > 0:
            ylabel = ""
        else:
            ylabel = "Latitude"

        # Dismiss label of x-axis except for bottom row
        if i < figrows - 1:
            xlabel = ""
        else:
            xlabel = "Longitude"

        plot_single_map(
            axs[i, j],
            ds_data[ens],
            mask,
            ens,
            map_lons_extend,
            map_lats_extend,
            vmax=vmax,
            xlabel=xlabel,
            ylabel=ylabel,
        )

    # Remove axes from empty subplots
    for n in range(len(ens_names), figrows * figcols):
        i = n // figcols
        j = n % figcols
        plt.sca(axs[i, j])
        plt.axis("off")

    fig.savefig(fname, dpi=300, bbox_inches="tight")


def plotting_sources_cases(
    ds_data,
    mask,
    ens_names,
    figwidth=24,
    figheight=14,
    vmax=5,
    central_longitude=0,
    figrows=5,
    figcols=5,
    map_lons_extend=[-85, 40],
    map_lats_extend=[10, 80],
    fname="fig",
):
    my_projection = crs.PlateCarree(central_longitude=central_longitude)
    rows = figrows
    cols = figcols
    # Make figure
    fig, axs = plt.subplots(
        rows,
        cols,
        figsize=(figwidth, figheight),
        subplot_kw={"projection": my_projection},
    )
    i = 0
    j = 0
    for iens, ens in enumerate(ens_names):
        print("------  Plotting", ens)

        ds_data[ens].plot(
            ax=axs[i, j],
            vmin=0,
            vmax=vmax,
            robust=False,
            cmap=cm.rain,
            transform=crs.PlateCarree(),
            extend="max",
            cbar_kwargs=dict(fraction=0.05, shrink=0.5, label=None),
        )
        axs[i, j].set_title(ens, loc="left")

        if len(mask["mask"].values.shape) > 2:
            maskvals = mask["mask"].values[0, :]
        else:
            maskvals = mask["mask"].values

        axs[i, j].contour(
            mask["lon"].values,
            mask["lat"].values,
            maskvals,
            colors=["r"],
            transform=crs.PlateCarree(),
        )
        axs[i, j].add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        axs[i, j].add_feature(cartopy.feature.BORDERS, linestyle="-", linewidth=0.2)

        if central_longitude == 0:
            lons_ticks = np.arange(-180, 181, 20)
        else:
            lons_ticks = np.arange(0, 361, 20)

        axs[i, j].set_xticks(lons_ticks, crs=crs.PlateCarree())
        axs[i, j].set_yticks(np.arange(-90, 91, 10), crs=crs.PlateCarree())

        lon_formatter = LongitudeFormatter(direction_label=False, degree_symbol="")
        axs[i, j].xaxis.set_major_formatter(lon_formatter)
        fig.axes[iens].set_extent(
            [
                map_lons_extend[0],
                map_lons_extend[1],
                map_lats_extend[0],
                map_lats_extend[1],
            ],
            crs.PlateCarree(),
        )

        # Dismiss label of y-axis, except for left most column

        if j > 0:
            axs[i, j].set_ylabel("")
        else:
            axs[i, j].set_ylabel("Latitude")

        if i == rows - 1:
            axs[i, j].set_xlabel("Longitude")
        else:
            axs[i, j].set_xlabel("")

        if j < cols - 1:
            i = i
            j = j + 1
        else:
            i = i + 1
            j = 0

    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close()


def plotting_sources_one_case(
    ds_data,
    mask,
    ens_names,
    figwidth=24,
    figheight=14,
    vmax=5,
    central_longitude=0,
    figrows=1,
    figcols=1,
    map_lons_extend=[-85, 40],
    map_lats_extend=[10, 80],
    fname="fig",
):
    my_projection = crs.PlateCarree(central_longitude=central_longitude)
    rows = figrows
    cols = figcols
    # Make figure

    for iens, ens in enumerate(ens_names):
        fig, axs = plt.subplots(
            rows,
            cols,
            figsize=(figwidth, figheight),
            subplot_kw={"projection": my_projection},
        )

        print("------  Plotting", ens)

        ds_data[ens].plot(
            ax=axs,
            vmin=0,
            vmax=vmax,
            robust=False,
            cmap=cm.rain,
            transform=crs.PlateCarree(),
            extend="max",
            cbar_kwargs=dict(fraction=0.05, shrink=0.5, label=None),
        )
        axs.set_title(ens, loc="left")

        if len(mask["mask"].values.shape) > 2:
            maskvals = mask["mask"].values[0, :]
        else:
            maskvals = mask["mask"].values

        axs.contour(
            mask["lon"].values,
            mask["lat"].values,
            maskvals,
            colors=["r"],
            transform=crs.PlateCarree(),
        )
        axs.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        axs.add_feature(cartopy.feature.BORDERS, linestyle="-", linewidth=0.2)

        if central_longitude == 0:
            lons_ticks = np.arange(-180, 181, 20)
        else:
            lons_ticks = np.arange(0, 361, 20)

        axs.set_xticks(lons_ticks, crs=crs.PlateCarree())
        axs.set_yticks(np.arange(-90, 91, 10), crs=crs.PlateCarree())

        lon_formatter = LongitudeFormatter(direction_label=False, degree_symbol="")
        axs.xaxis.set_major_formatter(lon_formatter)
        fig.axes[0].set_extent(
            [
                map_lons_extend[0],
                map_lons_extend[1],
                map_lats_extend[0],
                map_lats_extend[1],
            ],
            crs.PlateCarree(),
        )

        # Dismiss label of y-axis, except for left most column

        axs.set_ylabel("Latitude")
        axs.set_xlabel("Longitude")

        fig.savefig(fname + "_" + ens + ".png", dpi=300, bbox_inches="tight")
        plt.close()


def plotting_global(means, cases, cmaps, vmaxs, masks, plot_vline, proj, fname):
    prj = crs.PlateCarree()
    ## Robinson Projection
    if proj == 1:
        prjr = crs.PlateCarree(
            central_longitude=70
        )  # crs.Robinson(central_longitude=70)
    elif proj == 2:
        prjr = crs.Robinson(central_longitude=70)
    else:
        print("Only proj=1 and proj=2 are allowed")
        quit()
    fig, ax = plt.subplots(figsize=(14, 10), subplot_kw={"projection": prjr})

    ax1 = fig.add_subplot(131)
    ax1.set_axis_off()
    ax2 = fig.add_subplot(132)
    ax2.set_axis_off()
    ax3 = fig.add_subplot(133)
    ax3.set_axis_off()

    axs = [ax1, ax2, ax3]

    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)

    for i, case in enumerate(cases):
        print("-----> Plotting", case)

        mask = xr.open_dataset(case + "/" + masks[case])

        values = means[case].values
        values[values <= 0.25] = np.nan

        means[case][:, :] = values

        vars()["p_" + case] = means[case].plot.contourf(
            ax=ax,
            levels=np.arange(0, vmaxs[i] + int(vmaxs[i] / 5), int(vmaxs[i] / 5)),
            extend="max",
            cmap=cmaps[i],
            transform=prj,
            add_colorbar=False,
        )

        fig.colorbar(
            vars()["p_" + case],
            ax=axs[i],
            pad=0.05,
            orientation="horizontal",
            label="mm",
        )

        # axs[i].set_title(case+" case")

        if len(mask["mask"].values.shape) > 2:
            maskvals = mask["mask"].values[0, :]
        else:
            maskvals = mask["mask"].values

        ax.contour(
            mask["lon"].values,
            mask["lat"].values,
            maskvals,
            colors=["r"],
            transform=prj,
        )

    gl = ax.gridlines(draw_labels=True, linewidth=0)
    gl.top_labels = False
    gl.right_labels = False

    if proj == 1:
        if plot_vline:
            ax.axvline(-50, color="gray", linewidth=1)
            ax.axvline(50, color="gray", linewidth=1)

        ax.set_extent([-80, 230, -60, 80], crs=prj)
        glons = [
            -80,
            -60,
            -40,
            -20,
            0,
            20,
            40,
            60,
            80,
            100,
            120,
            140,
            160,
            180,
            -160,
            -140,
        ]
        gl.xlocator = mticker.FixedLocator(glons)

    if proj == 2:
        plt.subplots_adjust(bottom=-0.01)

    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close()
