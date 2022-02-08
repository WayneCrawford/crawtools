# Copyright 2021 Wayne Crawford
#
# This file is based on OBStools:TFNoise.
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
import xarray as xr
from matplotlib import pyplot as plt

# from obstools.atacr import utils
from .spectral_density import SpectralDensity

np.seterr(all='ignore')
# np.set_printoptions(threshold=sys.maxsize)


class TransferFunction(object):
    """
    Calculates and contains transfer functions for a given input channel.

    Attributes:
        xf (XArray): Transfer functions
        xferrs (XArray): Transfer function errors
        noise_chan: channel on which the noise is found
        n_winds: number of spectral windows used
    """
    def __init__(self, spect_density, in_chan, out_chans=None,
                 noise_chan='output', n_to_reject=3):
        """
        Args:
            spect_coher (:class:`.SpectralDensity`): Cross-spectral density
                matrix
            in_chan (str): input channel
            out_chans (list of str): output channels  (None => all))
            noise_chan (str): 'input', 'output', 'equal', 'unknown'
            n_to_reject (int): number of neighboring frequencies for which the
                coherence must be above the 95% significance level in order
                to use for cleaning (0 means use all frequencies)
        """
        if out_chans is None:
            out_chans = list(spect_density.sdf.coords['output'].values)
        self._check_init_arg_errors(spect_density, in_chan, out_chans)

        f = spect_density.sdf.coords['f'].values
        self.noise_chan = noise_chan
        self.n_winds = spect_density.n_winds
        shape = (1, len(out_chans), len(f))
        self.xfs = xr.DataArray(
            data=np.zeros(shape, dtype='complex'),
            dims=('input', 'output', 'f'),
            coords={'input': [in_chan], 'output': out_chans, 'f': f})
        self.xferrs = self.xfs.copy()
        for out in out_chans:
            xf, xferr = self._calcxf(spect_density, in_chan, out, noise_chan,
                                     n_to_reject)
            self.xfs.loc[dict(input=in_chan, output=out)] = xf
            self.xferrs.loc[dict(input=in_chan, output=out)] = xferr
         
    @property
    def f(self):
        """Return frequencies"""
        return self.sdf.coords['f'].values

    @property
    def xf(self, input_chan):
        """Return one transfer function"""
        return self.xfs.sel(input=input_chan)

    @property
    def xf_ptr(self, input_chan):
        """Return pointer to one transfer function"""
        return self.xfs.loc[dict(input=input_chan)]
        
    @property
    def xferr(self, input_chan):
        """Return one transfer function uncertainty"""
        return self.xferrs.sel(input=input_chan)

    @property
    def xferr_ptr(self, input_chan):
        """Return pointer to one transfer function uncertainty"""
        return self.xferrs.loc[dict(input=input_chan)]
        
    @property
    def coh_signif_95(self):
        return self.coh_signif(prob=0.95)

    def coh_signif(self, prob=0.95):
        """
        Definition: L_1(alpha, q) = sqrt(1-alpha**(1/q))
        
        where alpha = 1-prob and 2(q+1) = nwinds (degree of freedom)
        
        For nwinds >> 1, L1 ~ sqrt(1-alpha**(2/nwinds))
        For a 95% signif level this comes out to
            sqrt(1-.05**(2/nwinds)) for nwinds >> 1.
        I previously used sqrt(2/nwinds) for the 95% signif level (alpha=0.05),
        but L1 is much closer to sqrt(6/nwinds).
        
        Args:
            prob (float): significance level (between 0 and 1)
        """
        assert prob > 0 and prob < 1
        alpha = 1 - prob
        n = self.n_winds/2
        return sqrt(1 - alpha ** (1. / (n-1)))
        
        np.sqrt(2/self.n_winds)

    @staticmethod
    def _check_init_arg_errors(spect_density, in_chan, out_chans):
        if not isinstance(spect_density, SpectralDensity):
            raise TypeError("Error: A TransferFunc object must be initialized "
                            "with a SpectralDensity object")
        if not isinstance(out_chans, list):
            raise TypeError("Error: out_chans is not a list")
        if not isinstance(in_chan, str):
            raise TypeError("Error: in_chan is not a str")

    def _calcxf(self, spect_density, input, output, noise_chan="output",
                n_to_reject=1):
        """
        Calculate transfer function between a given input and output channel

        Args:
            spect_density(:class: ~SpectralDensity): cross-spectral density
                matrix
            input (str): input channel name
            output (str): output channel name
            noise_channel (str): which channel has the noise
            n_to_reject (int): only use values for which more than this
                many consecutive coherences are above the 95% significance level
                (0 = use all)
        """
        Gxx = spect_density.sdf.sel(input=input, output=input)
        Gxy = spect_density.sdf.sel(input=input, output=output)
        coh = spect_density.coherence(input, output)
        # Shouldn't need abs() here, coh is real positive
        coh_mag_sq = abs(coh) * abs(coh)
        H = Gxy / Gxx  # B&P Equation 6.69
        errbase = np.sqrt((np.ones(coh.shape) - coh_mag_sq)
                          / (2 * self.n_winds * coh_mag_sq))
        if n_to_reject > 0:
            goods = coh > self.coh_signif_95
            # for n == 1, should do nothing, for n == 2, shift once, etc
            for n in range(n_to_reject - 1):
                goods = np.logical_and(goods, np.concatenate((True, goods[1:])))   
            H[~goods] = 0
        if noise_chan == 'output':
            xf = H * coh
            xferr = np.abs(xf) * errbase
        elif noise_chan == 'input':
            xf = H / coh
            xferr = np.abs(xf) * errbase
        elif noise_chan == 'equal':
            xf = H
            xferr = np.abs(xf) * errbase
        elif noise_chan == 'unknown':
            xf = H
            # VERY ad-hoc error guesstimate
            maxerr = np.abs(coh**(-1)) + errbase
            minerr = np.abs(coh) - errbase
            xferr = np.abs(xf * (maxerr-minerr)/2)
        else:
            raise ValueError(f'unknown noise channel: "{noise_chan}"')
        return xf, xferr

    def plot(self):
        """
        Plot transfer functions

        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        inputs = list(self.xfs.coords['input'].values)
        outputs = list(self.xfs.coords['output'].values)
        rows = len(inputs)
        cols = len(outputs)
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        for input, i in zip(inputs, range(len(inputs))):
            for output, j in zip(outputs, range(len(outputs))):
                axa, axp = self.plot_one(input, output, fig, (rows, cols),
                                         (i, j), ylabel=j == 0,
                                         xlabel=i == rows-1)
            ax_array[i, j] = (axa, axp)
        return ax_array

    def plot_one(self, in_chan, out_chan, fig=None, fig_grid=(1, 1),
                 plot_spot=(0, 0), xlabel=True, ylabel=True):
        """
        Plot one transfer function

         Arguments:
            in_chan (str): input channel
            out_chan (str): output channel
            fig (:class: ~matplotlib.figure.Figure): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            subplot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            xlabel (bool): put an xlabel on this subplot
            ylabel (bool): put a y label on this subplot

         Returns:
            tuple:
                transfer function amplitude plot
                transfer function phase plot
            """
        print(f'{self.xfs=}')
        xf = self.xfs.sel(input=in_chan, output=out_chan).copy()
        print(f'{xf=}')
        xferr = self.xferrs.sel(input=in_chan, output=out_chan).copy()
        f = self.xfs.coords['f'].values
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        fig.suptitle("Transfer Functions")
        ax_a = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                (3*plot_spot[0]+0, plot_spot[1]+0),
                                rowspan=2)
        xf[xf == 0] = None
        ax_a.vlines(f, np.abs(xf-xferr), np.abs(xf+xferr),
                    label=f"'{out_chan}' / '{in_chan}'")
        ax_a.set_xscale('log')
        ax_a.set_yscale('log')
        # ax_a.loglog(f, np.abs(xf), label=f"'{out_chan}' / '{in_chan}'")
        ax_a.set_xlim(f[1], f[-1])

        # legend_1 = ax_a.legend()
        if ylabel:
            ax_a.set_ylabel('TF')
        else:
            ax_a.set_yticklabels([])
        ax_a.set_xticklabels([])
        # Plot phase
        ax_p = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                (3*plot_spot[0]+2, plot_spot[1]+0))
        ax_p.semilogx(f, np.degrees(np.angle(xf)), marker='.', linestyle='')
        ax_p.set_ylim(-180, 180)
        ax_p.set_xlim(f[1], f[-1])
        ax_p.set_yticks((-180, 0, 180))
        if ylabel:
            ax_p.set_ylabel('Phase')
        else:
            ax_p.set_yticklabels([])
        if xlabel:
            ax_p.set_xlabel('Frequency (Hz)')
        return ax_a, ax_p

    def save(self, filename):
        """
        Method to save the object to file using `~Pickle`.

        Args:
            filename (str): File name

        Examples
        --------

        Run demo through all methods

        >>> from obstools.atacr import DayNoise, StaNoise, TFNoise
        >>> daynoise = DayNoise('demo')
        Uploading demo data - March 04, 2012, station 7D.M08A
        >>> daynoise.QC_daily_spectra()
        >>> daynoise.average_daily_spectra()
        >>> tfnoise_day = TFNoise(daynoise)
        >>> tfnoise_day.transfer_func()
        >>> stanoise = StaNoise('demo')
        Uploading demo data - March 01 to 04, 2012, station 7D.M08A
        >>> stanoise.QC_sta_spectra()
        >>> stanoise.average_sta_spectra()
        >>> tfnoise_sta = TFNoise(stanoise)
        >>> tfnoise_sta.transfer_func()

        Save object

        >>> tfnoise_day.save('tf_daynoise_demo.pkl')
        >>> tfnoise_sta.save('tf_stanoise_demo.pkl')

        Check that everything has been saved

        >>> import glob
        >>> glob.glob("./tf_daynoise_demo.pkl")
        ['./tf_daynoise_demo.pkl']
        >>> glob.glob("./tf_stanoise_demo.pkl")
        ['./tf_stanoise_demo.pkl']

        """

        if not self.transfunc:
            print("Warning: saving before having calculated the transfer "
                  "functions")

        # Remove traces to save disk space
        file = open(filename, 'wb')
        pickle.dump(self, file)
        file.close()
