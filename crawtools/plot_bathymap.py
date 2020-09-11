#!/usr/bin/env python3
"""
Plot a bathymetry map like in GMT
"""
import math
from datetime import datetime

import numpy as np
from scipy.io import netcdf
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import gradient, sin, cos, arctan2, arctan, cm, pi, hypot

# CONSTANTS
text_bbox = dict(boxstyle="round", pad=0.1, edgecolor='none',
                 facecolor='white', alpha=0.5)


def plot_basemap(grid_x, grid_y, bounds, bathy_map, contour=1000, pastel=False, 
                 zoom_shade=True, debug=True):
    """
    Set up map, plot coastlines (plus bathy if requested) and stations

    :param grid_x: grid spacing in x direction
    :param grid_y: grid spacing in y direction
    :param bounds [list]: [minlon,maxlon,minlat,maxlat]
    :param bathy_map [boolean]: True: use etopo bathy map
                                False: no bathy map
                                [str]:       filename of netcdf bathymetry file
    :param pastel: lighten colors to make pastel
    :param contour: contour spacing (m)
    :param zoom_shade: calculate shading using only plotted part of bathy map
    :returns: map object
    """

    if debug:
        print("In plot_basemap()")
    lon_min, lon_max, lat_min, lat_max = (bounds[0], bounds[1], bounds[2],
                                          bounds[3])
    # Plot coastlines
    if debug:
        print("Plotting coastlines")
    m = _get_basemap_coastlines(lon_min, lat_min, lon_max, lat_max)
    if debug:
        print("Plotting parallels and meridians")
    m.drawparallels(np.arange(math.floor(lat_min / grid_y) * grid_y,
                              math.ceil(lat_max / grid_y) * grid_y, grid_y),
                    labels=[1, 0, 0, 0])  # draw parallels
    m.drawmeridians(np.arange(math.floor(lon_min / grid_x) * grid_x,
                              math.ceil(lon_max / grid_x) * grid_x, grid_x),
                    labels=[0, 0, 0, 1])  # draw meridians

    if not bathy_map:
        pass
    elif isinstance(bathy_map, str):
        # Read in bathy data and make shaded version
        if debug:
            print(f'Plotting bathymetry map "{bathy_map}"')
        f = netcdf.netcdf_file(bathy_map, 'r')
        if 'x' in f.variables:
            if debug:
                print(f'x = {f.variables["x"]}')
                print(f'y = {f.variables["x"]}')
            lon = f.variables['x'][:].copy()
            lat = f.variables['y'][:].copy()
            z = f.variables['z'][:].copy()
        elif 'x_range' in f.variables:
            if debug:
                print(f'x_range = {f.variables["x_range"].data}')
                print(f'y_range = {f.variables["y_range"].data}')
                print(f'spacing = {f.variables["spacing"].data}')
            x_spacing, y_spacing = f.variables['spacing'].data
            lon = np.arange(f.variables['x_range'].data[0],
                            f.variables['x_range'].data[1] + x_spacing/2,
                            x_spacing)
            lat = np.arange(f.variables['y_range'][0],
                            f.variables['y_range'][1]+y_spacing/2,
                            y_spacing)
            # z = np.flipud(np.reshape(f.variables['z'].data, ((len(lat),len(lon))))
            z = np.flipud(np.array(f.variables['z'].data).reshape(len(lat),len(lon)))
        [xx, yy] = np.meshgrid(lon, lat)
        f.close()
        # Get indices corresponding to map ranges (necessary because extent
        # in imshow leaves an offset)
        ixmin = np.flatnonzero((lon >= lon_min)).min()
        ixmax = np.flatnonzero((lon <= lon_max)).max()
        iymin = np.flatnonzero((lat >= lat_min)).min()
        iymax = np.flatnonzero((lat <= lat_max)).max()
        if debug:
            print('z={}x{}, ixmin,max = {}, {}, iymin,max = {}, {}'.format(
                len(z), len(z[0]), ixmin, ixmax, iymin, iymax), flush=True)
        if pastel:
            im_mult, im_offset = 0.5, 0.5
        else:
            im_mult, im_offset = 1.0, 0.0
        if zoom_shade:
            z_shade = _set_shade(np.nan_to_num(z[iymin:iymax,ixmin:ixmax]),
                                 cmap=cm.jet, scale=1.0, azdeg=90)
            m.imshow(im_mult * z_shade + im_offset)
            print(xx[ixmin:ixmax].shape)
            m.contour(xx[iymin:iymax,ixmin:ixmax],
                      yy[iymin:iymax,ixmin:ixmax],
                      z[iymin:iymax,ixmin:ixmax],
                      [-5000, -4000, -3000, -2000, -1000, -1],
                      latlon=True, colors='k', linestyles='solid', linewidths=1)
        else:
            z_shade = _set_shade(np.nan_to_num(z),
                                 cmap=cm.jet, scale=1.0, azdeg=90)
            m.imshow(im_mult * z_shade[iymin:iymax, ixmin:ixmax] + im_offset)
            m.contour(xx, yy, z, [-5000, -4000, -3000, -2000, -1000, -1],
                      latlon=True, colors='k', linestyles='solid', linewidths=1)
        m.fillcontinents()
    else:
        if debug:
            print(f'Plotting etopo bathymetry map')
        m.etopo()
    m.drawcoastlines(linewidth=3.0, color="black")
    # m.drawcountries()

    return m


def plot_station(m, lon, lat, name='', sym='o', color='blue', ms=10, fs = 7,
                 text_offset = 0.01, text_color='black'):
    """
    plot a station
    
    :param m: map object
    :param lon: station longitude
    :param lat: station latitude
    :param name: station name
    :param sym: station symbol
    :param ms: station markersize
    :param fs: fontsize
    :param text_offset: offset of text from station (degrees)
    :param text_color: text color
    :returns: map object
    """
    m.plot(lon, lat, sym, color=color, markersize=ms, latlon=True)
    if name:
        x, y = m(lon + text_offset, lat)
        plt.text(x, y, name, fontsize=fs, va='center', color=text_color)
    return m


def _get_basemap_coastlines(lon_min, lat_min, lon_max, lat_max):
    """
    Get Basemap object at best available resolution
    """
    try:
        m = Basemap(lon_min, lat_min, lon_max, lat_max, projection="merc",
                    resolution='h')
    except:
        print('  mpl_toolkit Basemap high resolution coastline not found, '
              'trying intermediate')
        try:
            m = Basemap(lon_min, lat_min, lon_max, lat_max, projection="merc",
                        resolution='i')
        except:
            print('  mpl_toolkit Basemap intermediate resolution coastline '
                  'not found, using low')
            m = Basemap(lon_min, lat_min, lon_max, lat_max, projection="merc",
                        resolution='l')
    return m


def save_map(filename, title, fontsize=12, show=False, debug=False):
    """
    Saves map to a file

    :param filename:
    :param title: plot title
    :param timedelta_total: total time spent (datetime.timedelta)
    :returns: figure object, base_name of output file
    """
    if debug:
        print("In close_map()")
    base_name = filename
    plt.title(title, fontsize=fontsize)
    fig1 = plt.gcf()    # Needed to save after "show" (which creates new fig)
    if show:
        plt.show()
    fig1.savefig(filename)
    
def _set_shade(a, intensity=None, cmap=cm.jet, scale=10.0, azdeg=165.0,
              altdeg=45.0):
    '''
    sets shading for data array based on intensity layer or data value

    inputs:
        a - a 2-d array or masked array
        intensity - a 2-d array of same size as a (no check on that)
                    representing the intensity layer. if none is given
                    the data itself is used after getting the hillshade values
                    see hillshade for more details.
        cmap - a colormap (e.g matplotlib.colors.LinearSegmentedColormap
              instance)
        scale,azdeg,altdeg - parameters for hilshade function see there for
              more details
    output:
        rgb - an rgb set of the Pegtop soft light composition of the data and
                intensity can be used as input for imshow()
    based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight
    '''
    if intensity is None:
        # hilshading the data
        intensity = _hillshade(a, scale=scale, azdeg=azdeg, altdeg=altdeg)
    else:
        # or normalize the intensity
        intensity = (intensity - intensity.min())\
                    / (intensity.max() - intensity.min())
    # get rgb of normalized data based on cmap
    rgb = cmap((a - a.min()) / float(a.max() - a.min()))[:, :, :3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2 * d * rgb + (rgb**2) * (1 - 2*d)
    return rgb


def _hillshade(data, scale=10.0, azdeg=165.0, altdeg=45.0):
    '''
    Convert data to hillshade based on matplotlib.colors.LightSource class.

    input:
         data - a 2-d array of data
         scale - scaling value of the data. higher number = lower gradient
         azdeg - where the light comes from: 0 south ; 90 east ; 180 north ;
                      270 west
         altdeg - where the light comes from: 0 horizon ; 90 zenith
    output: a 2-d array of normalized hilshade
    '''
    # convert alt, az to radians
    az = azdeg * pi / 180.0
    alt = altdeg * pi / 180.0
    # gradient in x and y directions
    dx, dy = gradient(data / float(scale))
    slope = 0.5 * pi - arctan(hypot(dx, dy))
    aspect = arctan2(dx, dy)
    intensity = sin(alt) * sin(slope)\
        + cos(alt) * cos(slope) * cos(-az - aspect - 0.5 * pi)
    intensity = (intensity - intensity.min())\
        / (intensity.max() - intensity.min())
    return intensity

