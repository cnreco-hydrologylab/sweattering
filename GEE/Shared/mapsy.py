#this is basically the code that is now found in maps.py

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import scipy
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np

def draw_map(lat_list, lon_list, value_list, z_min = -999, z_max = -999, n_levels = -999, map_type='satellite', mesh_type='contourf', zoom = 15, save = False, ax = None, trajectory = False, log = None, fig = None):

    if n_levels == -999:
        n_levels = 10

    if z_min == -999:
        z_min = min(value_list)

    if z_max == -999:
        z_max = max(value_list)

    center_lon = np.asarray(lon_list).mean()
    center_lat = np.asarray(lat_list).mean()
    width = 1.1*max([max(lon_list)-min(lon_list), max(lat_list)-min(lat_list)])
    lon_min = center_lon-0.5*width
    lon_max = center_lon+0.5*width
    lat_min = center_lat-0.5*width
    lat_max = center_lat+0.5*width

    if map_type=='satellite':
        url = 'http://mt0.google.com/vt/lyrs=s&hl=en&x={x}&y={y}&z={z}'
    elif map_type=='map':
        url = 'http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z} '

    stamen_terrain = cimgt.GoogleTiles(url = url)#cimgt.Stamen('terrain-background')

    # Create a GeoAxes in the tile's projection.
    if ax == None:
        if fig == None:
            fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot(1, 1, 1, projection=stamen_terrain.crs)
        ax.set_aspect('equal')

    # gl = ax.gridlines(draw_labels=True, alpha=0.3, ls = ':', color = 'k')
    # gl.xlabels_top = gl.ylabels_right = False
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    # Limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add the Stamen data at zoom level 8.
    ax.add_image(stamen_terrain, zoom)


    lon_arr          = np.linspace(min(lon_list), max(lon_list), 500)
    lat_arr          = np.linspace(min(lat_list), max(lat_list), 500)
    lon_mesh, lat_mesh = np.meshgrid(lon_arr, lat_arr)
    MMGC_mesh = griddata((np.asarray(lon_list), np.asarray(lat_list)), np.asarray(value_list), (lon_mesh, lat_mesh), method='nearest')
    sigma = [5, 5]
    #sigma = [10, 10]
    #sigma = [0, 0]
#     MMGC_mesh = scipy.ndimage.filters.gaussian_filter(MMGC_mesh, sigma, mode='constant')
# more filters on scipy.ndimage site

    x, y = lon_mesh, lat_mesh

    if mesh_type=='contourf':
        ax.contourf(x, y, MMGC_mesh, alpha=0.8, cmap='Spectral_r', linestyles='None', vmin=z_min, vmax=z_max, levels = n_levels, transform=ccrs.PlateCarree())
        mesh = ax.contourf([x.min(), x.max()], [y.min(), y.max()], [[z_min, z_max],[z_min, z_max]], alpha=1, cmap='Spectral_r', linestyles='None', vmin=z_min, vmax=z_max, levels = n_levels, zorder = -5, transform=ccrs.PlateCarree())
    elif mesh_type=='mesh':
        ax.pcolormesh(x, y, MMGC_mesh, shading='nearest',cmap='Spectral_r', vmin=z_min, vmax=z_max, alpha=0.5, transform=ccrs.PlateCarree())
        mesh = ax.pcolormesh(x, y, MMGC_mesh, shading='nearest', cmap='Spectral_r', vmin=z_min, vmax=z_max, alpha=1, zorder = -5, transform=ccrs.PlateCarree())

    if trajectory:
        try:
            ax.plot(log['Lng'], log['Lat'], transform=ccrs.PlateCarree(), c  = 'k')
        except:
            raise Exception('Log not provided!')

    #plt.scatter(x_orig, y_orig, s=5*data[' Absolute Altitude (m)']+1, c='k')
    #fig.colorbar(contourf, boundaries=np.linspace(vmin, vmax, n_levels))
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False
    # gl.xlines = True
    # gl.xlocator = mticker.FixedLocator([120, 140, 160, 180, -160, -140, -120])
    # gl.ylocator = mticker.FixedLocator([0, 20, 40, 60])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'rotation': 25, 'ha':'right'}
    gl.ylabel_style = {'rotation': 25, 'ha':'right'}
    #gl.xlabel_style = {'rotation': 45}
    #gl.ylabel_style = {'rotation': 45}
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    #ax.scatter(lon_list, lat_list, c = value_list, transform=ccrs.PlateCarree(), zorder = 20, s = 100, cmap='Spectral_r', vmin=vmin, vmax=vmax)

    fig.colorbar(mesh)
    fig.tight_layout()
    if save:
        now = dt.datetime.now()
        current_time = now.strftime("%d_%m_%y_%H_%M_%S")
        fileNAME = f'Out_maps/map_{current_time}.png'
        fig.savefig(fileNAME, bbox_inches = 'tight',pad_inches = 0.3)
        return fileNAME

def draw_scatter(lat_list, lon_list, map_type='satellite', zoom = 15, save = False, **kwargs):

    center_lon = np.asarray(lon_list).mean()
    center_lat = np.asarray(lat_list).mean()
    width = 1.1*max([max(lon_list)-min(lon_list), max(lat_list)-min(lat_list)])
    lon_min = center_lon-0.5*width
    lon_max = center_lon+0.5*width
    lat_min = center_lat-0.5*width
    lat_max = center_lat+0.5*width

    fig = plt.figure(figsize=(12,8))
    if map_type=='satellite':
        url = 'http://mt0.google.com/vt/lyrs=s&hl=en&x={x}&y={y}&z={z}'
    elif map_type=='map':
        url = 'http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z} '

    stamen_terrain = cimgt.GoogleTiles(url = url)#cimgt.Stamen('terrain-background')

    # Create a GeoAxes in the tile's projection.
    ax = fig.add_subplot(1, 1, 1, projection=stamen_terrain.crs)
    ax.set_aspect('equal')

    # gl = ax.gridlines(draw_labels=True, alpha=0.3, ls = ':', color = 'k')
    # gl.xlabels_top = gl.ylabels_right = False
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    # Limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add the Stamen data at zoom level 8.
    ax.add_image(stamen_terrain, zoom)


    x, y = lon_list, lat_list
    mesh = ax.plot(x, y, transform=ccrs.PlateCarree(), **kwargs)

    #plt.scatter(x_orig, y_orig, s=5*data[' Absolute Altitude (m)']+1, c='k')
    #fig.colorbar(contourf, boundaries=np.linspace(vmin, vmax, n_levels))
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False
    # gl.xlines = True
    # gl.xlocator = mticker.FixedLocator([120, 140, 160, 180, -160, -140, -120])
    # gl.ylocator = mticker.FixedLocator([0, 20, 40, 60])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'rotation': 25, 'ha':'right'}
    gl.ylabel_style = {'rotation': 25, 'ha':'right'}
    #gl.xlabel_style = {'rotation': 45}
    #gl.ylabel_style = {'rotation': 45}
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    #ax.scatter(lon_list, lat_list, c = value_list, transform=ccrs.PlateCarree(), zorder = 20, s = 100, cmap='Spectral_r', vmin=vmin, vmax=vmax)

    fig.tight_layout()
    if save:
        now = dt.datetime.now()
        current_time = now.strftime("%d_%m_%y_%H_%M_%S")
        fileNAME = f'Out_maps/map_{current_time}.png'
        fig.savefig(fileNAME, bbox_inches = 'tight',pad_inches = 0.3)
        return fileNAME