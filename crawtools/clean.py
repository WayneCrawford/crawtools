# Copyright 2019 Pascal Audet & Helen Janiszewski
#
# This file is part of OBStools.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import pickle

import numpy as np
from matplotlib import pyplot as plt

from obstools.atacr import utils  # , plotting
from .day_noise import DayNoise
from .sta_noise import StaNoise

np.seterr(all='ignore')
# np.set_printoptions(threshold=sys.maxsize)


class Clean(object):
    """
    Calculate and apply transfer-functions to remove the effect of one
    channel on another

    Attributes
    ----------
    seq : list
        ordered list of transfer functions to apply to remove noise
    """
    def __init__(self, spect_coher, remove_order):
        """
        Calculates the sequence of transfer functions needed to remove noise
        
        Args:
            spect_coher: :class:`.SpectCoher`
                Spectra and coherencies
            remove_order: list of 'str'
                channels to "clean" from other channels, from first to last
        """
        # Argument checks
        if not isinstance(remove_order, 'list'):
            raise ValueError('remove_order is not a list')
        if not isinstance(spect_coher, SpectCoher):
            raise ValueError('spect_coher is not a SpectCoher object')
        for r in remove_order:
            if not isinstance(r, str):
                raise ValueError('remove_order contains non-string elements')
            if not r in spec_coher.drive.coords:
                raise ValueError(f'"{r}" is not a channel in spec_coher')
        
        # construct the transfer function sequence
        self.seq = []
        for r in remove_order:
            xf = TransferFunc(spect_coher, drive=r)
            spect_coher = self._clean_spect([xf])
            self.seq.extend(xf)
            
    def _clean_spect(self, spect_coher):
        """
        Return spectra cleaned using self sequence
        """
        return self._clean_spect(self.seq)
    
    def _clean_spect(self, spect_coher, xf_seq):
        """
        Return spectra cleaned using the given sequence of transfer functions

        Arguments:
        ----------
        spect_coher : :class:`.SpectCoher`
        xf_seq: list
            Sequence of TransferFunction objects to use to clean data.  Each
            TransferFunction can only have one drive channel
        """
        s_c = spect_coher.copy()
        if len(xf.drive) < 1:
            raise ValueError('transfer function has more than one drive channel')
        for xf in xf_seq:
            for r = xf.resp:
                d = xf.drive[0]
                s_c.spectra.sel(drive=r, resp=r) = spect_coher(drive=r, resp=r)\
                    - xf.sel(d = drive, r=resp) * spect_coher(drive=d, resp=d)
                # How to recalculate cross-spectra?

    def plot(self, key, coher_too=True):
        """
        """
        pass

    def plot_one(self, key, subkey, fig=None, fig_grid=(1, 1),
                 plot_spot=(0, 0), xlabel=True, ylabel=True,
                 coher_too=True):
        """
        """
        pass

    def save(self, filename):
        """
        """
        pass