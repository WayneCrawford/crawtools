#!/usr/bin/env python3
"""
Class to plot GMT-like bathymetry maps
"""
import numpy as np
from scipy.io import netcdf
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.feature import GSHHSFeature
from pylab import gradient, sin, cos, arctan2, arctan, cm, pi, hypot


class BathyMap():
    """
    Create a bathymetry map figure and axis
    """
    def __init__(self, map_extent, bathy_file,
                 grid_x=None, grid_y=None, intens_file=None):
        """
        Set up a bathymetric map

        :param map_extent [list]: [minlon,maxlon,minlat,maxlat]
        :param bathy_file [str]: name of netcdf bathymetry file
        :param intens_file: name of netcdf intensity file
        :kind intens_file: str, optional
        :param grid_x: x grid spacing for the plot axis
        :param grid_y: y grid spacing for the plot axis
        """
        self.map_extent = map_extent
        # self.lon, self.lat, self.z = self._setup_bathy_vars()
        self.lon, self.lat, self.z = self._read_bathy_map(bathy_file)
        self.intens = None
        if intens_file:
            lon, lat, intens = self._read_bathy_map(intens_file)
            if lon==self.lon and lat==self.lat:
                self.intens = intens
                
        # set up figure and axes
        self.fig = plt.figure()
        self.ax = plt.axes(projection=ccrs.Mercator())
        self.ax.set_extent(self.map_extent)
        if grid_x:
            min_x = grid_x * np.floor(min(self.lon) / grid_x)
            grid_x = np.arange(min_x, max(self.lon), grid_x)
        if grid_y:
            min_y = grid_y * np.floor(min(self.lat) / grid_y)
            grid_y = np.arange(min_y, max(self.lat), grid_y)
        gl = self.ax.gridlines(xlocs=grid_x, ylocs=grid_y, draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

    def _read_bathy_map(self, fname):
        """
        Read netcdf grid file
        
        :param fname: filename
        """
        f = netcdf.netcdf_file(fname, 'r')
        if 'x' in f.variables:
            lon = f.variables['x'][:].copy()
            lat = f.variables['y'][:].copy()
            z = f.variables['z'][:].copy()
        elif 'x_range' in f.variables:
            x_spacing, y_spacing = f.variables['spacing'].data
            lon = np.arange(f.variables['x_range'].data[0],
                            f.variables['x_range'].data[1] + x_spacing/2,
                            x_spacing)
            lat = np.arange(f.variables['y_range'][0],
                            f.variables['y_range'][1]+y_spacing/2,
                            y_spacing)
            z = np.flipud(np.array(f.variables['z'].data).reshape(len(lat),
                                                                  len(lon)))
        f.close()
        # Cut map down to desired range
        ixmin = np.flatnonzero((lon >= self.map_extent[0])).min()
        ixmax = np.flatnonzero((lon <= self.map_extent[1])).max()
        iymin = np.flatnonzero((lat >= self.map_extent[2])).min()
        iymax = np.flatnonzero((lat <= self.map_extent[3])).max()
        return lon[ixmin:ixmax], lat[iymin:iymax], z[iymin:iymax, ixmin:ixmax]

    # def _setup_bathy_vars(self):
    #     """
    #     Sets up bathymetric map variables
    #     """
    #     if not self.bathy_file:
    #         return
    #     f = netcdf.netcdf_file(self.bathy_map, 'r')
    #     if 'x' in f.variables:
    #         lon = f.variables['x'][:].copy()
    #         lat = f.variables['y'][:].copy()
    #         z = f.variables['z'][:].copy()
    #     elif 'x_range' in f.variables:
    #         x_spacing, y_spacing = f.variables['spacing'].data
    #         lon = np.arange(f.variables['x_range'].data[0],
    #                         f.variables['x_range'].data[1] + x_spacing/2,
    #                         x_spacing)
    #         lat = np.arange(f.variables['y_range'][0],
    #                         f.variables['y_range'][1]+y_spacing/2,
    #                         y_spacing)
    #         z = np.flipud(np.array(f.variables['z'].data).reshape(len(lat),
    #                                                               len(lon)))
    #     f.close()
    #     # Cut map down to desired range
    #     ixmin = np.flatnonzero((lon >= self.map_extent[0])).min()
    #     ixmax = np.flatnonzero((lon <= self.map_extent[1])).max()
    #     iymin = np.flatnonzero((lat >= self.map_extent[2])).min()
    #     iymax = np.flatnonzero((lat <= self.map_extent[3])).max()
    #     return lon[ixmin:ixmax], lat[iymin:iymax], z[iymin:iymax, ixmin:ixmax]

    def plot_image(self, pastel=False):
        """
        Plot the bathymetric image

        :param pastel: lighten colors to make pastel
        """
        if self.intens:
            z_shade = _set_shade(self.z, self.intens)
        else:
            z_shade = _set_shade(self.z)
        if pastel:
            z_shade = 0.5 * z_shade + 0.5
        self.ax.imshow(z_shade, origin='lower', extent=self.map_extent,
                       transform=ccrs.PlateCarree())

    def plot_contours(self, levels=500, linewidth=1, color='k'):
        """
        Plot the bathymetric contours

        :param levels: list of contours, or contour interval (m)
        :param linewidth: contour linewidth (1)
        :param colors: contour line color ('k')
        """
        if not isinstance(levels, list):
            interval = levels
            min_level = interval * np.floor(np.amin(self.z)/interval)
            levels = np.arange(min_level, np.amax(self.z), interval)
        plt.contour(self.lon, self.lat, self.z,
                    levels, colors=color,
                    linestyles='solid', linewidths=linewidth,
                    transform=ccrs.PlateCarree())

    def show(self):
        """
        Show the plot on the screen
        """
        plt.show()

    def plot_coastlines(self, resolution):
        """
        Plot coastlines

        :param resolution: what resolution coastlines to include.  Must
            correspond to a NaturalEarth resolution ('10m', '50m', '110m')
            or a GSHSS resolution ('auto', 'low', 'high', 'full'...)
        """
        NE_coast_resolutions = ['10m', '50m', '110m']
        GSHSS_coast_resolutions = ['auto', 'coarse', 'low',
                                   'intermediate', 'high',
                                   'full']
        if resolution in GSHSS_coast_resolutions:
            coast = GSHHSFeature(scale=resolution)
            self.ax.add_feature(coast)
        elif resolution in NE_coast_resolutions:
            coast = NaturalEarthFeature(scale=resolution)
            self.ax.add_feature(coast)
            # self.ax.coastlines(resolution=resolution)
        else:
            print(f'Invalid coastline resolution: "{resolution}"')

    def plot_station(self, lon, lat, name='', sym='o', color='deepskyblue',
                     mec='k', ms=8, fs=7, text_offset=0.01,
                     text_color='black', text_box=True,
                     **kwargs):
        """
        plot a station

        :param m: map object
        :param lon: station longitude
        :param lat: station latitude
        :param name: station name
        :param sym: station symbol
        :param ms: station markersize
        :param mec: marker edge color
        :param fs: fontsize
        :param text_offset: offset of text from station (degrees)
        :param text_color: text color
        :param text_box: Put a box around the text
        :param kwargs: keyword arguments to pass on to the ax.plot()
        """
        self.ax.plot(lon, lat, sym, color=color, ms=ms, mec=mec,
                     transform=ccrs.PlateCarree(),
                     **kwargs)
        if name:
            bbox = None
            if text_box:
                bbox = dict(boxstyle="round", pad=0.1, edgecolor='none',
                            facecolor='white', alpha=0.5)
            self.ax.text(lon + text_offset, lat, name,
                         fontsize=fs, va='center', bbox=bbox,
                         color=text_color, transform=ccrs.PlateCarree())

    def save_map(self, filename, title, fontsize=12):
        """
        Saves map to a file

        :param filename:
        :param title: plot title
        :param timedelta_total: total time spent (datetime.timedelta)
        :returns: figure object, base_name of output file
        """
        self.ax.set_title(title, fontsize=fontsize)
        self.fig.savefig(filename)


def _set_shade(a, intensity=None, cmap=cm.jet, **kwargs):
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
        intensity = _hillshade(a, **kwargs)
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
