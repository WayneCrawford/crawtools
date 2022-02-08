"""
Spectral Density Functions
"""
import pickle
from dataclasses import dataclass

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from obspy.core.inventory import Inventory
from obstools.atacr import DayNoise, StaNoise
from .utils import _prol1pi, _prol4pi
from .transfer_function import TransferFunction

np.seterr(all='ignore')
# np.set_printoptions(threshold=sys.maxsize)


class XFCleaner():
    """
    Class for calculating and applying TransferFunction-based data cleaning
    """
    def __init__(sdf, remove_list, noise_channel, n_to_reject=3):
        """
        Args:
            sdf (:class: ~SpectralDensity): one-sided spectral density functs
            remove_list (list): list of channels to remove, in order
            noise_channel (str): which channel the noise is on.  Choices are:
                "input", "output", "equal", "unknown", "model"
            n_to_reject (int): Number of neighboring frequencies for which
                the coherence must be above the 95% signifcance level in order
                to use for cleaning (0 means use all frequencies)
        Returns:
            clean_list (list): list of transfer functions to apply, in order
        """
        sdf = sdf.copy()
        self.xf_list = []
        ignore_chans = []
        for in_chan in remove_list:
            ignore_chans += in_chan
            out_chans = [x in sdf.channels if x not in ignore_chans]
            xf = TransferFunc(sdf, in_chan, None, noise_channel, n_to_reject)
            sdf = self._removeXFMultiplied(sdf, xf)
            self.xf_list += xf

    def apply(self, sdf, clean_list):
        """
        Args:
            sdf (:class: ~SpectralDensity): one-sided spectral density functs
            clean_list (list): list of transfer functions to apply, in order
        Return:
            new SpectralDensity object
        """
        newsdf = sdf.copy()
        for xf in self.xf_list:
            newsdf = self._removeXFMultiplied(newsdf, xf)
        return newsdf        

    def apply_stream(self, stream, transfer_function_list):
        """
        Args:
            stream (Stream): list of channels to remove, in order
            noise_channel (str): which channel the noise is on.  Choices are:
                "input", "output", "equal", "unknown", "model"
            n_to_reject (int): Number of neighboring frequencies for which
                the coherence must be above the 95% signifcance level in order
                to use for cleaning (0 means use all frequencies)
        """
        print('apply_stream() is not yet implemented')

    @staticmethod
    def _removeXFMultiplied(sdf, xfs):
        """
        Clean a channel or channels by removing another channel multiplied
        by the given transfer function(s)
        """
        in_chan = list(xfs.xfs.coords['input'].values)[0]:
        in_auto = sdf.autospect(in_chan)
        for out_chan in list(xfs.xfs.coords['output'].values):
            out_auto_ptr = sdf.autospect_ptr(out_chan)
            cross_ptr = sdf.crossspect_ptr(in_chan, out_chan)
            # Would prefer not to have to calculate this next
            crossT_ptr = sdf.crossspect_ptr(out_chan, in_chan)
            xf = xfs.xf(in_chan)
            out_auto_ptr -= in_auto * np.abs(xf)**2  # B&P equation 6.35
            cross_ptr -= out_auto * xf               # B&P equation 6.36
            crossT_ptr = np.conj(out_auto * xf)
