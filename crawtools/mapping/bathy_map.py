#!/usr/bin/env python3
"""
Class to plot GMT-like bathymetry maps
"""
import warnings

import numpy as np
from scipy.io import netcdf
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.feature import GSHHSFeature
from pylab import gradient, arctan, cm, hypot
# import xarray as xr


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
        if not isinstance(fname, str):
            if fname is True:
                print('Sorry, ETOPO isn"t working')
            return self.map_extent[:2], self.map_extent[2:4], None
            
            # etopo = xr.open_dataset('etopo1.nc')
            # # Slice desired coordinates.
            # etopo = etopo.sel(x=slice(self.map_extent[0],self.map_extent[1]),
            #                   y=slice(self.map_extent[2],self.map_extent[3]))  
            # # Return sliced variables.
            # return etopo.x, etopo.y, etopo.z
        # Open ETOPO1 file and slice desired area.

        try:
            f = netcdf.netcdf_file(fname, 'r')
            if 'x' in f.variables:
                lon = f.variables['x'][:].copy()
                lat = f.variables['y'][:].copy()
                z = f.variables['z'][:].copy()
            elif 'x_range' in f.variables:
                x_spacing, y_spacing = f.variables['spacing'].data
                lon = np.arange(f.variables['x_range'].data[0],
                                f.variables['x_range'].data[1] + x_spacing / 2,
                                x_spacing)
                lat = np.arange(f.variables['y_range'][0],
                                f.variables['y_range'][1] + y_spacing / 2,
                                y_spacing)
                z = np.flipud(np.array(f.variables['z'].data).reshape(len(lat),
                                                                      len(lon)))
            f.close()
        except Exception as e:
            print(e)
            return None, None, None

        # Cut map down to desired range
        lon_step= lon[1] - lon[0]
        lat_step= lat[1] - lat[0]
        if self.map_extent[0] < min(lon):
            warnings.warn(f'Requested min lon < map minlon: increasing to {min(lon)}')
            self.map_extent[0] = min(lon) + lon_step/10
        if self.map_extent[1] > max(lon):
            warnings.warn(f'Requested max lon > map maxlon: decreasing to {max(lon)}')
            self.map_extent[1] = max(lon) - lon_step/10
        if self.map_extent[2] < min(lat):
            warnings.warn(f'Requested min lat < map minlat: increasing to {min(lat)}')
            self.map_extent[2] = min(lat) + lat_step/10
        if self.map_extent[3] > max(lat):
            warnings.warn(f'Requested max lat > map maxlat: decreasing to {max(lat)}')
            self.map_extent[3] = max(lat) - lat_step/10
        ixmin = np.flatnonzero((lon >= self.map_extent[0])).min()
        ixmax = np.flatnonzero((lon <= self.map_extent[1])).max()
        iymin = np.flatnonzero((lat >= self.map_extent[2])).min()
        iymax = np.flatnonzero((lat <= self.map_extent[3])).max()
        return lon[ixmin:ixmax], lat[iymin:iymax], z[iymin:iymax, ixmin:ixmax]

    def plot_image(self, cmap=cm.jet, pastel=False, shaded=False, intens_file=None, **kwargs):
        """
        Plot the bathymetric image

        :param cmap: colormap (e.g matplotlib.colors.LinearSegmentedColormap
              instance)
        :param pastel: lighten colors to make pastel
        :param shaded: [bool] apply shading
        :param intens_file: file containing intensities to plot (generally for shading)
        :param **kwargs: arguments to pass to hillshade (scale, azdeg, altdeg)
        """
        intens = None
        if intens_file is not None:
            lon, lat, intens = self._read_bathy_map(intens_file)
            if not (len(lon) == len(self.lon)) or not (len(lat) == len(self.lat)):
                warnings.warn('intensity lon/lat values different length than in bathy map')
            elif not all(abs(lon - self.lon) < ((lon[1] - lon[0]) / 10)) \
                    or not all(abs(lat - self.lat) < ((lat[1] - lat[0]) / 10)):
                warnings.warn('intensity lon/lat values do not match bathy map')
        if intens is not None:
            z_shade = _set_shade(self.z, intensity=intens, cmap=cmap, **kwargs)
        else:
            z_shade = _set_shade(self.z, cmap=cmap, shaded=shaded, **kwargs)
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
            min_level = interval * np.floor(np.nanmin(self.z)/interval)
            levels = np.arange(min_level, np.nanmax(self.z), interval)
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


    def plot(self, lon, lat, sym='o', color='deepskyblue',
                     mec='k', ms=8, **kwargs):
        """
        plot a point

        :param lon:  longitude
        :param lat:  latitude
        :param sym:  symbol
        :param color: marker face color
        :param kwargs: keyword arguments to pass on to pyplot.plot()
        """
        self.ax.plot(lon, lat, sym, color=color,
                     transform=ccrs.PlateCarree(),
                     **kwargs)

    def arrow(self, x, y, dx, dy, **kwargs):
        """
        plot a arrow

        :param x:  longitude
        :param y:  latitude
        :param dx:  arrow length
        :param dy:  arrow length
        :param kwargs: keyword arguments to pass on to pyplot.arrow()
        """
        self.ax.arrow(x, y, dx, dy, transform=ccrs.PlateCarree(), **kwargs)

    def text(self, lon, lat, name, text_box=True,
                   ha='left', va='center', color='black', **kwargs):
        """
        plot a station label

        :param m: map object
        :param lon: longitude
        :param lat: latitude
        :param name: station name
        :param text_box: Put a semi-opaque box under the text
        :param ha: horizontal text alignment ('left', 'center', 'right')
        :param va: vertical text alignment ('center', 'middle', 'bottom')
        :param color: text color
        :param kwargs: keyword arguments to pass on to pyplot.text()
        """
        if text_box:
            kwargs['bbox'] = dict(boxstyle="round", pad=0.1, ec='none', fc='white', alpha=0.5)
        self.ax.text(lon, lat, name, va=va, ha=ha,
                     color=color, transform=ccrs.PlateCarree(), **kwargs)

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


def _set_shade(a, intensity=None, cmap=cm.jet, shaded=True, **kwargs):
    '''
    sets shading for data array based on intensity layer or data value

    inputs:
        a - a 2-d array or masked array
        intensity - a 2-d array of same size as a, representing the intensity
                    layer. if none is given the data itself is used after
                    getting the hillshade values see hillshade for more details.
        cmap - a colormap (e.g matplotlib.colors.LinearSegmentedColormap
              instance)
        shaded: boolean calculate shading
        scale,azdeg,altdeg - parameters for hilshade function see there for
              more details
    output:
        rgb - an rgb set of the Pegtop soft light composition of the data and
                intensity can be used as input for imshow()
    based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight
    '''
    if intensity is None:
        if shaded is True:
            # hilshade the data
            intensity = _hillshade(a, **kwargs)
        else:
            intensity = np.ones(a.shape)
    else:
        # or normalize the intensity
        intensity -= np.nanmin(intensity)
        intensity /= np.nanmax(intensity)
    # get rgb of normalized data based on cmap
    rgb = cmap((a - np.nanmin(a)) / float(np.nanmax(a) - np.nanmin(a)))[:, :, :3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2 * d * rgb + (rgb**2) * (1 - 2*d)
    return rgb


def _hillshade(data, scale=1.0, azdeg=165.0, altdeg=45.0):
    '''
    Convert data to hillshade based on matplotlib.colors.LightSource class.

    input:
         data - a 2-d array of data
         scale - scaling value of the data. higher number = lower gradient
         azdeg - where the light comes from: 0 south ; 90 east ; 180 north ;
                      270 west
         altdeg - where the light comes from: 0 horizon ; 90 zenith
    output: a 2-d shading array (normalized 0 to 1)
    '''
    # gradient in x and y directions
    dx, dy = gradient(data / float(scale))
    # slope = np.pi/2. - arctan(hypot(dx, dy))
    slope = np.pi/2. - np.arctan(np.sqrt(dx*dx + dy*dy))
    aspect = np.arctan2(dx, dy)

    # convert alt, az to radians
    azimuth = np.radians(azdeg)
    altitude = np.radians(altdeg)
    intensity = np.sin(altitude) * np.sin(slope)\
        + np.cos(altitude) * np.cos(slope)\
        * np.cos((-azimuth - np.pi/2.) - aspect)
    # Scale from 0 to 1
    intensity -= np.nanmin(intensity)
    intensity /= np.nanmax(intensity)
    return intensity
