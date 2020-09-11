import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.feature import GSHHSFeature
import matplotlib.pylab as plt
import numpy as np

NE_coast_resolutions = ['10m', '50m', '110m']
GSHSS_coast_resolutions = ['auto', 'coarse', 'low', 'intermediate', 'high', 'full']

class MapCrossSection():
    """
    Create axes map with cross-sections to right (east) and below (south)

    ax_map is a Mercator projection, when you plot to it you must
        specify transform=cartopy.crs.Geodetic()
    """
    def __init__(self, figsize=None, fill_figure=True, pad_left=0.6,
                 pad_right=0.6, pad_top=0.4, pad_bottom=0.4,
                 map_extent=None, depth_extent=None, extent_buffer=0.2,
                 latitudes=None, longitudes=None, depths=None, gridlines=True,
                 coastlines=None):
        """
        :param figsize: width, height in inches
        :kind figsize: (float, float), optional
        :param fill_figure: Expand map bounds to fill the figure
        :param pad_left: inches to left of axes for labels/text/emptiness
        :param pad_right: inches to right of axes for labels/text/emptiness
        :param pad_top: inches above axes for labels/text/emptiness
        :param pad_bottom: inches below axes for labels/text/emptiness
        :param map_extent: bounds of map (x0, x1, y0, y1) or None
        :param depth_extent: depth bounds (z0, z1) or None
        :param latitudes: list of "event" latitudes, used to calculate map
            extents if map_extents == None
        :param longitudes: list of "event" longitudes, used to calculate map
            extents if map_extents == None
        :param longitudes: list of "event" depths, used to calculate depth
            extents if depth_extents == None
        :param extent_buffer: extents will exceed event limits by this
            fraction of the event range
        :param gridlines: Put gridlines on map view (boolean)
        :param coastlines: what resolution coastlines to include.  Must
            correspond to a NaturalEarth resolution ('10m', '50m', '110m')
            or a GSHSS resolution ('auto', 'low', 'high', 'full'...)
        """
        self.km_per_deg_lat = 1.852 * 60.
        self.map_extent = self._get_map_extent(map_extent, latitudes,
                                               longitudes, extent_buffer)
        self.depth_extent = self._get_depth_extent(depth_extent, depths,
                                                   extent_buffer)
        self.pad_left = pad_left
        self.pad_right = pad_right
        self.pad_top = pad_top
        self.pad_bottom = pad_bottom
        
        self.coastlines = None
        if coastlines is not None:
            if ((coastlines in NE_coast_resolutions)
                or (coastlines in GSHSS_coast_resolutions)):
                    self.coastlines=coastlines
            else:
                print(f'Invalid coastlines string: "{coastlines}"')

        fig = plt.figure(figsize)
        rect_map, rect_xz, rect_zy = self._make_rects(fig, fill_figure)

        # Setup Lat-Z axis (lower left)
        self.ax_xz = plt.axes(rect_xz,
                              xlim=self.map_extent[:2],
                              ylim=self.depth_extent[::-1])
        self.ax_xz.tick_params(direction='in')

        # Setup Z-lon axis (upper right)
        self.ax_zy = plt.axes(rect_zy,
                              xlim=self.depth_extent,
                              ylim=self.map_extent[2:])
        self.ax_zy.tick_params(direction='in')
        self.ax_zy.xaxis.set_ticks_position('top')
        self.ax_zy.yaxis.set_ticks_position('right')

        # Setup map axis (upper left)
        self.ax_map = plt.axes(rect_map, projection=ccrs.Mercator())
        self.ax_map.set_extent(self.map_extent)
        if self.coastlines in GSHSS_coast_resolutions:
            coast = GSHHSFeature(scale=self.coastlines)
            feature = self.ax_map.add_feature(coast)
        elif self.coastlines in NE_coast_resolutions:
            self.ax_map.coastlines(resolution='10m')
        gl = self.ax_map.gridlines(draw_labels=True)
        gl.right_labels = False
        gl.bottom_labels = False
        if not gridlines:
            gl.xlines = False
            gl.ylines = False

    @property
    def lat_min(self):
        return self.map_extent[2]

    @property
    def lat_max(self):
        return self.map_extent[3]

    @property
    def lon_min(self):
        return self.map_extent[0]

    @property
    def lon_max(self):
        return self.map_extent[1]

    @property
    def depth_min(self):
        return self.depth_extent[0]

    @property
    def depth_max(self):
        return self.depth_extent[1]

    @property
    def km_per_deg_lon(self):
        return self.km_per_deg_lat\
            * np.cos(np.radians((self.lat_min + self.lat_max) / 2.))

    @property
    def lonrange(self):
        return self.lon_max - self.lon_min

    @property
    def latrange(self):
        return self.lat_max - self.lat_min

    @property
    def depthrange(self):
        return self.depth_extent[1] - self.depth_extent[0]

    @property
    def lonrange_km(self):
        return self.lonrange * self.km_per_deg_lon

    @property
    def latrange_km(self):
        return self.latrange * self.km_per_deg_lat

    @property
    def xrange_km(self):
        return self.lonrange_km + self.depthrange

    @property
    def yrange_km(self):
        return self.latrange_km + self.depthrange

    @property
    def plot_aspect(self):
        return self.xrange_km / self.yrange_km

    def _get_map_extent(self, map_extent, latitudes, longitudes,
                        extent_buffer):
        if map_extent:
            assert(len(map_extent) == 4),\
                f'len(map_extent) != 4 ({map_extent})'
            # print(map_extents)
            assert map_extent[1] > map_extent[0],\
                f'map_extent: x0({map_extent[0]:g}) >= x1({map_extent[1]:g})'
            assert map_extent[3] > map_extent[2], 'map_extent: y0 >= y1'
            lon_min, lon_max, lat_min, lat_max = map_extent
        else:
            assert latitudes is not None and longitudes is not None,\
                'need either extent or latitudes and longitudes'
            lon_min, lon_max = min(longitudes), max(longitudes)
            lat_min, lat_max = min(latitudes), max(latitudes)
            lon_buffer = (lon_max - lon_min) * extent_buffer
            lat_buffer = (lat_max - lat_min) * extent_buffer
            lon_min -= lon_buffer
            lon_max += lon_buffer
            lat_min -= lat_buffer
            lat_max += lat_buffer
            map_extent = (lon_min, lon_max, lat_min, lat_max)
        return map_extent

    def _get_depth_extent(self, depth_extent, depths, extent_buffer):
        if depth_extent:
            assert depth_extent[1] > depth_extent[0], 'depth_extent: z0 >= z1'
        else:
            assert depths is not None, 'need either depth_extents or depths'
            depth_min, depth_max = min(depths), max(depths)
            depth_buffer = (depth_max - depth_min) * extent_buffer
            depth_min -= depth_buffer
            depth_max += depth_buffer
            depth_extent = (depth_min, depth_max)
        return depth_extent

    def _adjustminmax(self, fig_aspect):
        if fig_aspect > self.plot_aspect:   # expand x_range
            add_kms = ((self.yrange_km * fig_aspect) - self.xrange_km) / 2.
            add_deg = add_kms / self.km_per_deg_lon
            self.map_extent[0] -= add_deg
            self.map_extent[1] += add_deg
        else:   # expand y_range
            add_kms = ((self.xrange_km / fig_aspect) - self.yrange_km) / 2.
            add_deg = add_kms / self.km_per_deg_lat
            self.map_extent[2] -= add_deg
            self.map_extent[3] += add_deg
        return

    def _make_figure(self, figsize):
        fig = plt.figure(figsize)
        fig_inches = fig.get_size_inches()
        fig_xrange = fig_inches[0] - self.pad_left - self.pad_right
        fig_yrange = fig_inches[1] - self.pad_top - self.pad_bottom
        return fig, fig_aspect

    def _make_rects(self, fig, fill_figure):
        """
        Make rectangles for the three axes
        
        :param fig: figure to make rects in
        :param fill_figure: if true, extend map bounds to fill figure
        """
        fig_inches_x, fig_inches_y = fig.get_size_inches()
        fig_xrange = fig_inches_x - self.pad_left - self.pad_right
        fig_yrange = fig_inches_y - self.pad_top - self.pad_bottom
        fig_aspect = fig_xrange / fig_yrange

        if fill_figure:
            # Expand map extent to fill figure
            self._adjustminmax(fig_aspect)
        
        inch_per_km_x = fig_xrange / self.xrange_km
        inch_per_km_y = fig_yrange / self.yrange_km
        if self.plot_aspect > fig_aspect:
            inch_per_km_y *= fig_aspect / self.plot_aspect
        elif self.plot_aspect < fig_aspect:
            inch_per_km_x *= self.plot_aspect / fig_aspect

        # make axes rects
        left_width = self.lonrange_km * inch_per_km_x / fig_inches_x
        right_width = self.depthrange * inch_per_km_x / fig_inches_x
        lower_height = self.depthrange * inch_per_km_y / fig_inches_y
        upper_height = self.latrange_km * inch_per_km_y / fig_inches_y
        
        left_left = self.pad_left / fig_inches_x
        right_left = left_left + left_width
        lower_bottom = self.pad_bottom / fig_inches_y
        upper_bottom = lower_bottom + lower_height

        rect_map = [left_left, upper_bottom, left_width, upper_height]
        rect_xz = [left_left, lower_bottom, left_width, lower_height]
        rect_zy = [right_left, upper_bottom, right_width, upper_height]
        return rect_map, rect_xz, rect_zy


if __name__ == "__main__":
    pt = (45.5, -13.0, 35)
    for map_extent in [[45.0, 45.5, -13.2, -12.4],
                       [45.0, 46.0, -13.2, -12.4],
                       [44.5, 46.5, -13.2, -12.4]]:
        for fill_fig in [False, True]:
            mcs = MapCrossSection(map_extent=map_extent, depth_extent=[0, 50],
                                  fill_figure=fill_fig)
            mcs.ax_map.plot(pt[0], pt[1], '+', transform=ccrs.Geodetic())
            mcs.ax_xz.plot(pt[0], pt[2], '+')
            mcs.ax_zy.plot(pt[2], pt[1], '+')
            plt.show()
