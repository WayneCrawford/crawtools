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

np.seterr(all='ignore')
# np.set_printoptions(threshold=sys.maxsize)


@dataclass
class SpectralDensity():
    """
    Contains the spectral density matrix for different signals.
    """
    sdf: xr.DataArray    # one-sided spectral density functions
    windowtype: str      # window type
    starttimes: list     # list of starttimes (UTCDateTimes)
    inv: Inventory       # obspy inventory containing instrument responses
    n_winds: int = None  # Number of windows used to calculate spec densities

    @property
    def channels(self):
        assert (list(self.sdf.coords['input'].values)
                == list(self.sdf.coords['output'].values))
        return list(self.sdf.coords['input'].values)

    @property
    def f(self):
        """Return frequencies"""
        return self.sdf.coords['f'].values

    @property
    def autospect(self, channel):
        """Return an auto-spectral_density"""
        return self.sdf.sel(input=channel, output=channel)

    @property
    def crossspect(self, in_chan, out_chan):
        """Return a cross-spectral_density"""
        return self.sdf.sel(input=in_chan, output=out_chan)

    @classmethod
    def from_ATACR(cls, objnoise, horizontal_format='separate'):
        """
        Initiate class from ATACR DayNoise or StaNoise class

        Args:
            objnoise (:class:`.DayNoise` or :class:`.StaNoise`):
                noise spectra and frequencies
            horizontal_format (str): which type of horizontal channels to use:
                'aligned': one horizontal channel, aligned with highest noise
                'separate': two orthogonal horizontal channels
        """
        if (not objnoise and not isinstance(objnoise, DayNoise) and
                not isinstance(objnoise, StaNoise)):
            raise TypeError("Error: A TFNoise object must be initialized with"
                            " only one of type DayNoise or StaNoise object")

        if not objnoise.av:
            raise(Exception("Error: Noise object has not been processed (QC "
                            "and averaging) - aborting"))

        if horizontal_format == 'separate':
            chans = ['1', '2', 'Z', 'P']
        elif horizontal_format == 'aligned':
            chans = ['L', 'Z', 'P']
        else:
            raise ValueError('horizontal_format not "separate" or "aligned"')
        shape = (len(chans), len(chans), len(objnoise.f))
        sdf = xr.DataArray(data=np.zeros(shape, dtype='complex'),
                           dims=('input', 'output', 'f'),
                           coords={"input": chans,  "output": chans,
                                       "f": objnoise.f})
        if horizontal_format == 'separate':
            sdf.loc[dict(input='1', output='1')] = objnoise.power.c11
            sdf.loc[dict(input='1', output='2')] = objnoise.cross.c12
            sdf.loc[dict(input='1', output='Z')] = objnoise.cross.c1Z
            sdf.loc[dict(input='1', output='P')] = objnoise.cross.c1P
            sdf.loc[dict(input='2', output='1')] = np.conj(objnoise.cross.c12)
            sdf.loc[dict(input='2', output='2')] = objnoise.power.c22
            sdf.loc[dict(input='2', output='Z')] = objnoise.cross.c2Z
            sdf.loc[dict(input='2', output='P')] = objnoise.cross.c2P
            sdf.loc[dict(input='Z', output='1')] = np.conj(objnoise.cross.c1Z)
            sdf.loc[dict(input='Z', output='2')] = np.conj(objnoise.cross.c2Z)
            sdf.loc[dict(input='Z', output='Z')] = objnoise.power.cZZ
            sdf.loc[dict(input='Z', output='P')] = objnoise.cross.cZP
            sdf.loc[dict(input='P', output='1')] = np.conj(objnoise.cross.c1P)
            sdf.loc[dict(input='P', output='2')] = np.conj(objnoise.cross.c2P)
            sdf.loc[dict(input='P', output='Z')] = np.conj(objnoise.cross.cZP)
            sdf.loc[dict(input='P', output='P')] = objnoise.power.cPP
        elif horizontal_format == 'aligned':
            sdf.loc[dict(input='L', output='L')] = objnoise.rotation.cHH
            sdf.loc[dict(input='L', output='Z')] = objnoise.rotation.cHZ
            sdf.loc[dict(input='L', output='P')] = objnoise.rotation.cHP
            sdf.loc[dict(input='Z', output='L')] = np.conj(objnoise.rotation.cHZ)
            sdf.loc[dict(input='Z', output='Z')] = objnoise.power.cZZ
            sdf.loc[dict(input='Z', output='P')] = objnoise.cross.cZP
            sdf.loc[dict(input='P', output='L')] = np.conj(objnoise.rotation.cHP)
            sdf.loc[dict(input='P', output='Z')] = np.conj(objnoise.cross.cZP)
            sdf.loc[dict(input='P', output='P')] = objnoise.power.cPP
        if hasattr(objnoise, 'nwins'):
            n_winds = np.sum(objnoise.nwins)
        else:
            n_winds = np.sum(objnoise.goodwins)
        return cls(sdf, "hanning", None, None, n_winds)

    @classmethod
    def from_stream(cls, stream, window_s=2000, windowtype='hanning',
                    inv=None):
        """
        Calculate sdfs from the provided stream

        Should add a window selection algorithm, for now just steps by
        the window length
df
        Args:
            stream (~class `obspy.core.stream.Stream`): data
            window_s (float): desired window length in seconds
            windowtype (str): window type, must be a valid
        """
        stream = _align_traces(stream)

        # Select windows
        ws = window_s*stream[0].sample_rate
        ws = 2**(ws-1).bit_length()
        # window_starts = WindowSelect(stream, ws, windowtype)

        ft = {}
        ids = [tr.id for tr in stream]
        for tr in stream:
            ft[tr.id], f = _calculate_windowed_fft(tr, ws, ws, windowtype)
        sdf = xr.DataArray(data=np.zeros(shape=(len(ids), len(ids), len(f)),
                                         dtype='complex'),
                           dims=('input', 'output', 'f'),
                           coords={"input": ids,  "output": ids, "f": f})
        for inp in ids:
            for outp in ids:
                locid = dict(input=inp, output=outp)
                sdf.loc[locid] = np.abs(np.mean(ft[inp]*np.conj(ft[outp]),
                                                axis=0))[0:len(f)]
        return cls(sdf, windowtype, stream[0].stats.starttime, inv)

    def coherence(self, in_chan, out_chan):
        """
        Return the coherence for the given input and output channels

        This is a real-valued quantity
        Args:
            in_chan (str): input channel.  Must match one of the
                coordinates in sdf
            out_chan (str): output channel.  Must match one of the
                coordinates in sdf

        Bendat & Piersol (1986), eq 6.27
        """
        if in_chan not in self.sdf.input:
            raise ValueError('"in_chan" not in spectral density matrix')
        if out_chan not in self.sdf.output:
            raise ValueError('"out_chan" not in spectral density matrix')
        return (abs(self.sdf.sel(input=in_chan, output=out_chan))**2
                / (self.sdf.sel(input=in_chan, output=in_chan)
                   * self.sdf.sel(input=out_chan, output=out_chan)))

    def plot(self, x=None, overlap=False):
        self.plot_autospect()

    def plot_autospect(self, x=None, overlap=False, show_coher=False):
        """
        Plot autospectra

        Args:
            x (list of str): limit to the listed channels
            overlap (bool): put all spect on one axis
        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        if x is None:
            x = list(self.sdf.coords['input'].values)
        else:
            for key in x:
                if key not in self.transfunc:
                    ValueError('key "{key}" not in self.transfunc')
        n_subkeys = len(x)
        if n_subkeys == 1:
            rows, cols = 1, 1
        elif n_subkeys == 2:
            rows, cols = 1, 2
        elif n_subkeys <= 4:
            rows, cols = 2, 2
        elif n_subkeys <= 6:
            rows, cols = 2, 3
        else:
            rows, cols = 3, 3
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        fig.suptitle('Auto-spectra')
        for key, i in zip(x, range(len(x))):
            i_row = int(i/cols)
            i_col = i - cols*i_row
            ylabel = 'PSD' if i_col == 0 else None
            axa, axp = self.plot_one_spect(key, key, fig, (rows, cols),
                                           (i_row, i_col), ylabel=ylabel,
                                           show_xlabel=i_row == rows-1,
                                           show_coher=show_coher)
            ax_array[i_row, i_col] = (axa, axp)
        return ax_array

    def plot_xspect(self, x=None, overlap=False, show_coher=False):
        """
        Plot cross (and auto) spectra

        Args:
            x (list of str): limit to the listed channels
            overlap (bool): put all spect on one axis
        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        if x is None:
            x = list(self.sdf.coords['input'].values)
        else:
            for key in x:
                if key not in self.transfunc:
                    ValueError('key "{key}" not in self.transfunc')
        n_subkeys = len(x)
        rows, cols = n_subkeys, n_subkeys
        ax_array = np.ndarray((rows, cols), dtype=tuple)
        fig, axs = plt.subplots(rows, cols, sharex=True)
        fig.suptitle('Cross-spectra')
        for drive, i in zip(x, range(len(x))):
            for resp, j in zip(x, range(len(x))):
                title = resp if i == 0 else None
                ylabel = drive if j == 0 else None
                axa, axp = self.plot_one_spect(drive, resp, fig, (rows, cols),
                                               (i, j), ylabel=ylabel,
                                               show_xlabel=i == rows-1,
                                               show_legend=False,
                                               title=title,
                                               show_coher=show_coher)
                ax_array[i, j] = (axa, axp)
        return ax_array

    def plot_one_spect(self, key, subkey, fig=None, fig_grid=(1, 1),
                       plot_spot=(0, 0), show_xlabel=True, ylabel=None,
                       show_legend=True, title=None, show_coher=False):
        """
        Plot one spectral density

        Args:
            key (str): input (driving) channel
            subkey (str): output (response) channel
            fig (:class: ~matplotlib.figure.Figure): figure to plot on, if
                None this method will plot on the current figure or create
                a new figure.
            fig_grid (tuple): this plot sits in a grid of this many
                              (rows, columns)
            subplot_spot (tuple): put this plot at this (row,column) of
                                  the figure grid
            show_xlabel (bool): put an xlabel on this subplot
            ylabel (str): y label on this subplot
            title (str): title to put on this subplot
            show_legend (bool): put a legend on this subplot
            show_coher (bool): draw coherency on the same plot

        Returns:
            tuple:
                amplitude plot axis
                phase plot axis
        """
        sdf = self.sdf.sel(input=key, output=subkey)
        f = self.sdf.coords['f'].values
        if fig is None:
            fig = plt.gcf()
        # Plot amplitude
        # print(f'{subkey=}, {plot_spot=}')
        ax_a = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                (3*plot_spot[0]+0, plot_spot[1]+0),
                                rowspan=2)
        if show_coher:
            ax2 = ax_a.twinx()
            ax2.semilogx(f, np.abs(self.coherence(key, subkey)),
                         color='red', linewidth=0.5, alpha=0.8)
            ax2.axhline(np.sqrt(2/self.n_winds), color='red', linewidth=0.5,
                        alpha=0.8, ls='--')
            ax2.set_ylim(0, 1)
            if plot_spot[1] == fig_grid[1]-1:  # Rightmost column
                ax2.set_ylabel('Coher', color='red')
            else:
                ax2.set_yticklabels([])
        sdf[sdf == 0] = None
        # print(f'{self.f[0:2]=},{self.f[-2:-1]=}')
        ax_a.loglog(f, np.abs(sdf), label=subkey)
        ax_a.set_xlim(f[1], f[-1])

        if show_legend:
            legend_1 = ax_a.legend()
            if show_coher:
                legend_1.remove()
                ax2.add_artist(legend_1)
        if ylabel:
            ax_a.set_ylabel(ylabel)
        else:
            ax_a.set_yticklabels([])
        ax_a.set_xticklabels([])
        if title:
            ax_a.set_title(title)
        # Plot phase
        ax_p = plt.subplot2grid((3*fig_grid[0], 1*fig_grid[1]),
                                (3*plot_spot[0]+2, plot_spot[1]+0))
        ax_p.semilogx(f, np.degrees(np.angle(sdf)), marker='.', linestyle='')
        ax_p.set_ylim(-180, 180)
        ax_p.set_xlim(f[1], f[-1])
        ax_p.set_yticks((-180, 0, 180))
        if ylabel:
            # ax_p.set_ylabel('Phase')
            pass
        else:
            ax_p.set_yticklabels([])
        if show_xlabel:
            ax_p.set_xlabel('Frequency (Hz)')
        return ax_a, ax_p

    def plot_cohers(self, x=None, y=None, overlap=False):
        """Plot coherences

        Args:
            x (list of str): limit to the listed x coordinates
            x (list of str): limit to the listed y coordinates
            overlap (bool): put all coher on one plot

        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        print('plot_cohers() is not yet implemented')

    def plot_one_coher(self, x=None, y=None, overlap=False):
        """Plot coherences

        Args:
            x (list of str): limit to the listed x coordinates
            x (list of str): limit to the listed y coordinates
            overlap (bool): put all coher on one plot

        Returns:
            (numpy.ndarray): array of axis pairs (amplitude, phase)
        """
        print('plot_one_coher() is not yet implemented')

    def save(self, filename):
        """
        Method to save the object to file using `~Pickle`.

        Parameters
        ----------
        filename : str
            File name

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


def _align_traces(stream):
    """Trim stream so that all traces are aligned and same length"""
    # Verify that all traces have the same sample rate
    first_start = last_start = stream[0].starttime
    first_end = last_end = stream[0].endtime
    sample_rate = stream[0].stats.sample_rate
    for tr in stream[1:]:
        if tr.starttime > last_start:
            last_start = tr.starttime
        elif tr.starttime < first_start:
            first_start = tr.starttime
        if tr.endtime < first_end:
            first_end = tr.endtime
        elif tr.endtime > last_end:
            last_end = tr.endtime
        if not tr.stats.sample_rate == sample_rate:
            raise ValueError("not all traces have same sample rate")

    if last_start >= first_end:
        raise ValueError("There are non-overlapping traces")
    if last_start - first_start > 1/sample_rate:
        print("Cutting up to {last_start-first_start}s from trace starts")
    if last_end - first_end > 1/sample_rate:
        print("Cutting up to {last_start-first_start}s from trace ends")
    stream.trim(last_start, first_end)
    min_len = min([tr.stats.npts for tr in stream])
    max_len = max([tr.stats.npts for tr in stream])
    if not max_len == min_len:
        for tr in stream:
            tr.data = tr.data[:min_len]
    return stream


# ## COPIED FROM ATACR, but other tapers added/allowed
def _calculate_windowed_fft(trace, ws, ss=None, win_taper='hanning'):
    """
    Calculates windowed Fourier transform

    Args:
        trace (:class:`~obspy.core.Trace`): Input trace data
        ws (int): Window size, in number of samples
        ss (int): Step size, or number of samples until next window
        win_taper (str): taper to apply to data ['hanning', 'prol4pi',
            'prol1pi']

    Returns:
        ft (:class:`~numpy.ndarray`): Fourier transform of trace
        f (:class:`~numpy.ndarray`): Frequency axis in Hz
    """
    n2 = _npow2(ws)
    f = trace.stats.sampling_rate/2. * np.linspace(0., 1., int(n2/2) + 1)
    # Extract sliding windows
    tr, nd = _sliding_window(trace.data, ws, ss, win_taper)
    # Fourier transform
    ft = np.fft.fft(tr, n=n2)
    return ft, f


def _sliding_window(a, ws, ss=None, win_taper='hanning'):
    """
    Split a data array into overlapping, tapered sub-windows

    Args:
        a (:class:`~numpy.ndarray`): 1D array of data to split
        ws (int): Window size in samples
        ss (int): Step size in samples. If not provided, window and step size
            are equal.
        win_taper (str): taper to apply to data ['hanning', 'prol4pi',
            'prol1pi', 'bartlett', 'blackman', 'hamming']

    Returns:
        out (:class:`~numpy.ndarray`): 1D array of windowed data
        nd (int): Number of windows
    """
    if ss is None:
        # no step size was provided. Return non-overlapping windows
        ss = ws
    # Calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ss
    nd = (valid) // ss
    out = np.ndarray((nd, ws), dtype=a.dtype)
    if win_taper in ['hanning', 'hamming', 'blackman', 'bartlett']:
        eval(f'taper = np.{win_taper}(ws)')
    elif win_taper == 'prol1pi':
        taper = _prol1pi(ws)
    elif win_taper == 'prol4pi':
        taper = _prol4pi(ws)
    else:
        raise ValueError(f'Unknown taper type "{win_taper}"')
    if nd == 0:
        out = a * taper
    for i in range(nd):
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        out[i] = a[start: stop] * taper
    return out, nd


def _npow2(x):
    return 1 if x == 0 else 2**(x-1).bit_length()
