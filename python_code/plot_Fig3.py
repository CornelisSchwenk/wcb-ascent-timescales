import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
import glob
import cartopy.crs as crs
import datetime as dt
import sys
from matplotlib.colors import BoundaryNorm
from matplotlib.collections import LineCollection

plt.style.use("ggplot")
#mpl.use('TkAgg')

path = "../nest_data/"
icon_path = "path/to/icon/output/"

#   These functions will only work if the Eulerian output is available!
#   If you dont have the Eulerian data, rewrite the functions (they are
#   not that long) so that you don't use the constfl2 and constfl3. The
#   only thing missing in the final plot will then be the domain borders. 
#   So make function load_data only return pres, lon and lat and make the 
#   function plot_Fig3 only take pres, lon and lat. 

def load_data():
    files = glob.glob(path + "WCB_tau_vars_normed.nc")
    data0 = Dataset(files[0]).variables
    pres = data0["p"][:] / 100
    z = data0["alt"][:]
    lon = data0["lon"][:] * (180 / np.pi)
    lat = data0["lat"][:] * (180 / np.pi)
    time = data0["rtime"][:]

    constfl2 = icon_path + "ICONCONST_DOM02_DOM02_ML_0001.nc"
    constfl3 = icon_path + "ICONCONST_DOM03_DOM03_ML_0001.nc"

    lon2 = Dataset(constfl2).variables["lon"][:]
    lat2 = Dataset(constfl2).variables["lat"][:]

    lon3 = Dataset(constfl3).variables["lon"][:]
    lat3 = Dataset(constfl3).variables["lat"][:]

    lon[lon > 100] = lon[lon > 100] - 360

    return pres, lon, lat, lon2, lat2, lon3, lat3

def plot_Fig3(pres, lon, lat, lon2, lat2, lon3, lat3):

    fig = plt.figure(figsize=(16, 10))
    ax1 = plt.axes(projection=crs.PlateCarree(central_longitude=0))
    ax1.set_extent([-75, 30, 25, 85], crs=crs.PlateCarree())
    gl = ax1.gridlines(crs=crs.PlateCarree(), draw_labels=True, linewidth=1,
                       color='grey', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.right_labels = False

    cmap = plt.get_cmap('viridis_r').copy()
    cmap.set_over('k')
    cmap.set_under([1, 0.7, 0])
    levels = np.arange(150, 1100, 50)
    norm = BoundaryNorm(levels, cmap.N, clip=False)

    for x in range(0, len(lon[0, ::10])): #plots only every tenth trajectory
        lon_plt = lon[:, x]
        lat_plt = lat[:, x]
        z_plt = pres[:, x]

        lon_plt = lon_plt[~np.isnan(lon_plt)]
        lat_plt = lat_plt[~np.isnan(lat_plt)]
        z_plt = z_plt[~np.isnan(z_plt)]

        points = np.array([lon_plt, lat_plt]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm, zorder=1)
        lc.set_array(z_plt[np.isfinite(z_plt)])
        lc.set_linewidth(1)
        ax1.add_collection(lc)

    # Plot ICON NESTS domain borders. 
    # Comment this part out if Eulerian data is not available
    # -----------
    ax1.plot([min(lon2), min(lon2)], [min(lat2), max(lat2)], color="r", lw=2)
    ax1.plot([max(lon2), max(lon2)], [min(lat2), max(lat2)], color="r", lw=2)
    ax1.plot([min(lon2), max(lon2)], [min(lat2), min(lat2)], color="r", lw=2)
    ax1.plot([min(lon2), max(lon2)], [max(lat2), max(lat2)], color="r", lw=2)

    ax1.plot([min(lon3), min(lon3)], [min(lat3), max(lat3)], color="k", lw=2)
    ax1.plot([max(lon3), max(lon3)], [min(lat3), max(lat3)], color="k", lw=2)
    ax1.plot([min(lon3), max(lon3)], [min(lat3), min(lat3)], color="k", lw=2)
    ax1.plot([min(lon3), max(lon3)], [max(lat3), max(lat3)], color="k", lw=2)

    ax1.plot([min(lon2), min(lon2)], [min(lat2), 50], color="black", lw=2, ls="--")
    ax1.plot([-20, -20], [min(lat2), 50], color="black", lw=2, ls="--")
    ax1.plot([min(lon2), -20], [min(lat2), min(lat2)], color="black", lw=2, ls="--")
    ax1.plot([min(lon2), -20], [50, 50], color="black", lw=2, ls="--")
    # -----------

    ax1.coastlines(resolution='auto', color='k', zorder=4)

    cb = fig.colorbar(lc, extend='both', shrink=0.8)
    cb.set_label(r'p [hPa]')
    cb.ax.invert_yaxis()

    save_title = "../plots/Fig3.png"
    fig.savefig(save_title, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    pres, lon, lat, lon2, lat2, lon3, lat3 = load_data()
    plot_spaghetti(pres, lon, lat, lon2, lat2, lon3, lat3)


