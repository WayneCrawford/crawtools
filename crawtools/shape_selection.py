"""
Classes to select events within different shapes
"""
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event import Catalog
import matplotlib.pyplot as plt
import matplotlib.path as paths
import matplotlib.patches as patches
import numpy as np


class ShapeSelection(object):
    """ Superclass for including or excluding events within the shape"""
    def include(self, catalog):
        """
        Return catalog with events within the shape

        :param catalog: earthquake catalog
        :kind catalog: ~class `obspy.core.event.catalog.Catalog`
        """
        new_catalog = Catalog()
        for event in catalog:
            origin = event.preferred_origin()
            if self.contains(origin.longitude, origin.latitude):
                new_catalog.append(event)
        return new_catalog

    def exclude(self, catalog):
        """
        Return catalog without any events in the shape

        :param catalog: earthquake catalog
        :kind catalog: ~class `obspy.core.event.catalog.Catalog`
        """
        new_catalog = Catalog()
        for event in catalog:
            origin = event.preferred_origin()
            if not self.contains(origin.longitude, origin.latitude):
                new_catalog.append(event)
        return new_catalog


class Circle(ShapeSelection):
    """
    Selects events within a given distance of a given point
    """
    def __init__(self, longitude, latitude, radius):
        """
        :param latitude: latitude at the center of the circle
        :param longitude: longitude at the center of the circle
        :param radius: radius of the circle in degrees latitude
        :kind latitude, longitude, radius_km: num

        >>> Circle(-12.3, 45.6, 10)
        Circle(-12.3, 45.6, 10)
        """
        self.longitude = longitude
        self.latitude = latitude
        self.radius = radius

    def __repr__(self):
        return "Circle({:g}, {:g}, {:g})".format(
            self.longitude, self.latitude, self.radius)

    def contains(self, event_longitude, event_latitude):
        """
        Returns True if event is inside or on the Circle, False if not

        :param event_latitude: latitude of the event
        :param event_longitude: longitude of the event
        :type event_latitude, event_longitude: numeric

        >>> Circle(-12.3, 45.6, 10)._iswithin(-12.35, 45.65)
        True
        >>> Circle(-12.3, 45.6, 10)._iswithin(-12.4, 45.7)
        False
        """
        dist = gps2dist_azimuth(event_latitude, event_longitude,
                                self.latitude, self.longitude)[0]
        return dist/(1852*60) <= self.radius

    def plot(self, ax=None, **kwargs):
        if not ax:
            fig, ax = plt.subplots()
        degrees_long_per_lat = 1./np.cos(np.radians(self.latitude))
        thetas = np.radians(np.arange(0, 361))
        ax.plot(self.longitude
                + self.radius*np.cos(thetas) * degrees_long_per_lat,
                self.latitude + self.radius*np.sin(thetas))

    def add_patch(self, ax=None, **kwargs):
        patch = patches.Circle((self.longitude, self.latitude),
                               self.radius, **kwargs)
        if not ax:
            fig, ax = plt.subplots()
        ax.add_patch(patch)


class Polygon(ShapeSelection):
    """
    Selection of events within or outside of a polygon
    """
    def __init__(self, polygon):
        """
        :param polygon: list of [longitudes, latitudes]

        No need to close the polygon
        """
        self.polygon = polygon

    def contains(self, event_longitude, event_latitude):
        """
        Returns True if event is inside or on the Polygon, False if not

        :param event_longitude: longitude of the event
        :param event_latitude: latitude of the event
        :type event_latitude, event_longitude: numeric

        >>> Polygon([[12, 45], [12, 47], [14, 47], [14, 45]])._iswithin(13, 46)
        True
        >>> Polygon([[12, 45], [12, 47], [14, 47], [14, 45]])._iswithin(15, 46)
        False
        """
        path = paths.Path(self.polygon)
        return path.contains_point((event_longitude, event_latitude))

    def plot(self, ax=None, **kwargs):
        """Plots the shape's outline"""
        if not ax:
            fig, ax = plt.subplots()
        p = self.polygon
        if not p[0] == p[-1]:
            p.append(p[0])
        ax.plot([x[0] for x in p], [x[1] for x in p], **kwargs)

    def add_patch(self, ax=None, **kwargs):
        """Adds a patch corresponding to the shape"""
        path = paths.Path(self.polygon)
        if not ax:
            fig, ax = plt.subplots()
        patch = patches.PathPatch(path, lw=2, **kwargs)
        ax.add_patch(patch)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
