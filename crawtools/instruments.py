"""
Classes used to plot instrument dynamic ranges
"""
from matplotlib import pyplot as plt
import numpy as np
from obspy.core.inventory.response import (Response, InstrumentSensitivity,
                                           PolesZerosResponseStage as PZStage,
                                           CoefficientsTypeResponseStage as CoeffStage)
from tiskitpy import Peterson_noise_model, Pressure_noise_model


# The following are useless but needed to work with a CoefficientsTypeResponseStage
coeff_args = {'decimation_input_sample_rate': 125,
              'decimation_factor': 1,
              'decimation_offset': 0,
              'decimation_delay': 0,
              'decimation_correction': 0}

class Sensor():
    def __init__(self, gain, zeros=[], poles=[], input_units='m/s',
                 self_noise=None, clip_level=None, gain_freq=10., norm_freq=1.,
                 comments=''):
        """
        Args:
            gain (float): Gain (V/{input_unit}})
            zeros (list): complex zeros (rad/s)
            poles (list): complex poles (rad/s)
            input_units (str): 'm/s' or 'Pa'
            self_noise (list): [[freqs],[vals]], vals in dB ref 1 (m/s^2 or Pa)^2/Hz
            clip_level (list): [[freqs],[vals]], vals in m/s^2 or Pa
            gain_freq (float): gain reference frequency (Hz)
            norm_freq (float): normalization frequency (Hz)
            comments (str): notes on the values provided
        """
        assert input_units in ['m/s', 'Pa']
        self.input_units = input_units
        # Use 'm/s' as input units because obspy doesn't know Pa
        self.response = PZStage(1, gain, gain_freq, 'm/s', 'V',
                                'LAPLACE (RADIANS/SECOND)', norm_freq,
                                zeros, poles)
        for testee in (self_noise, clip_level):
            if testee is not None:
                for x in testee:
                    assert len(x) == 2
        self.self_noise = self_noise
        self.clip_level = clip_level
        self.comments = comments

    def __str__(self):
        return ("Sensor: \n"
                f"    response={str(self.response)}"
        )

class Preamplifier():
    # No min or max values specified
    def __init__(self, gain, zeros=[], poles=[]):
        """
        Args:
            gain (float): Gain (V/V)
            zeros (list): complex zeros (rad/s)
            poles (list): complex poles (rad/s)
        """
        self.response = PZStage(2, gain, 1, 'V', 'V',
                                'LAPLACE (RADIANS/SECOND)', 0., zeros, poles)
    def __str__(self):
        return ("Preamplifier: \n"
                f"    response={str(self.response)}"
        )


class Logger():
    def __init__(self, fs_counts, fs_volts, dynamic_range_dB):
        """
        Args:
            fs_counts (float): Full-scale amplitude (total counts/2 for a signed A/D)
            fs_volts (float): Volts at the full-scale amplitude
            dynamic_range_dB (float): dynamic range in decibels
        """
        gain = fs_counts / fs_volts
        self.response = CoeffStage(3, gain, 1, 'V', 'count', 'DIGITAL',
                                   numerator=[], denominator=[], **coeff_args)
        # min and max values in counts
        self.max_counts = fs_counts
        self.min_counts = fs_counts / np.power(10., dynamic_range_dB/20)

    def __str__(self):
        return ("Logger: \n"
                f"    response={str(self.response)}"
                f"    max_counts={self.max_counts}\n"
                f"    min_counts={self.min_counts}"
        )

class Instrument():
    def __init__(self, sensor, preamplifier, logger, resp_freq=10.):
        self.sensor = sensor
        self.preamplifier = preamplifier
        self.logger = logger
        self.response = self._calculate_response(resp_freq)

    def __str__(self):
        return ("Instrument: \n"
                f"    sensor={self.sensor})\n"
                f"    preamplifier={self.preamplifier}\n"
                f"    logger={self.logger}\n"
                f"    resp={self.response}"
        )
    def _calculate_response(self, resp_freq):
        resp = Response(response_stages=[self.sensor.response,
                                         self.preamplifier.response,
                                         self.logger.response],
                        instrument_sensitivity = InstrumentSensitivity(
                            1., resp_freq,
                            self.sensor.response.input_units,
                            self.logger.response.output_units))
        resp.recalculate_overall_sensitivity(resp_freq)
        return resp

    @property
    def input_units(self):
        return self.sensor.input_units

    def plot(self, min_freq=0.001):
        self.response.plot(min_freq=min_freq)

    def max_values(self, frequencies):
        """
        In sensor.input_units
        """
        return self.logger.max_counts / self.evalresp(frequencies)

    def min_values(self, frequencies):
        """
        In 1/self.evalresp.output_units
        """
        # return self.logger.min_counts / self.evalresp(frequencies)
        return  (self.logger.min_counts/self.logger.max_counts) * self.max_values(frequencies)

    def max_dBs(self, frequencies, divisor=3.46):
        # divisor = sqrt(12) is for sine wave inputs(?)
        return 20*np.log10(np.abs(self.max_values(frequencies)/divisor))

    def min_dBs(self, frequencies, divisor=3.46):
        # divisor = sqrt(12) is for sine wave inputs(?)
        return 20*np.log10(np.abs(self.min_values(frequencies)/divisor))

    def evalresp(self, frequencies):
        """ in counts/pa or counts/m/s^2 """
        if self.sensor.input_units.lower() == 'pa':
            return self.response.get_evalresp_response_for_frequencies(frequencies, 'DEF')
        else:
            return self.response.get_evalresp_response_for_frequencies(frequencies, 'ACC')

    def plot_minmax(self, frequencies, color='b', plot_sense='f', show=True, ax=None,
                    label='range', ls = '-', alpha=0.05):
        """
        Plot the minimum and maximum PSD levels for an Instrument.
        
        This cannot be accurate because the PSD level depends on the variance,
        not the maximum value.  For example, a full-scale sine wave will have
        much less variance than randomly selected values in the full-scale range.
        
        If the sensor has a:
            - clip level, it will be plotted as a thick solid line
            - self-noise level, it will be plotted as a thick dashed line.
        
        Args:
            frequencies (:class: `numpy.array`): frequencies (Hz) at which to plot
            color (str): color of the line
            plot_sense (str): If starts with 'f', plot by frequencies
                              If starts with 'T', plot by periods
            show (bool): show the plot now
            ax (:class: `matplotlib.pyplot.Axis`): axis on which to plot
            label (str): label for this line
            ls (str): linestyle
            alpha (float): Transparency value for fill between the levels
        Returns:
            (:class: `matplotlib.pyplot.Axis`): axis plotted on
        """
        if ax is None:
            fig, ax = plt.subplots()
        x_ax_sense = plot_sense[0].lower()
        if x_ax_sense == 'f':
            x_values = frequencies
            x_label = 'Frequency (Hz)'
        elif x_ax_sense == 'p':
            x_values = np.power(frequencies, -1.)
            x_label = 'Period (s)'
        else:
            raise ValueError(f"{plot_sense=} doesn't start with an 'f' or a 'p'")

        ax.semilogx(x_values,
                    self.min_dBs(frequencies),
                    color=color, label=label, ls=ls)
        ax.semilogx(x_values, 
                    self.max_dBs(frequencies),
                    color=color, ls=ls)
        ax.fill_between(x_values,
                        self.min_dBs(frequencies),
                        self.max_dBs(frequencies),
                        alpha=alpha, color=color)
        if self.sensor.clip_level is not None:
            arr = self.sensor.clip_level
            if x_ax_sense == 'f':
                x_vals = [x[0] for x in arr]
            else:
                x_vals = [1. / x[0] for x in arr]
            ax.semilogx(x_vals,
                        20*np.log10(np.abs([x[1] for x in arr])),
                        color=color, ls='-', lw=3)
        if self.sensor.self_noise is not None:
            arr = self.sensor.self_noise
            if x_ax_sense == 'f':
                x_vals = [x[0] for x in arr]
            else:
                x_vals = [1. / x[0] for x in arr]
            ax.semilogx(x_vals, [x[1] for x in arr],
                        color=color, ls='--', lw=3)
        plt.legend()
        ax.set_xlabel(x_label)
        if self.sensor.input_units.lower() == 'pa':
            PSD_units = 'Pa^2/Hz'
        elif self.sensor.input_units.lower() in ('m', 'm/s', 'm/s^2'):
            PSD_units = '(m^2/s^4)/Hz'
        ax.set_ylabel(f'dB ref 1 {PSD_units}')
        if show:
            plt.show()
        return ax

    @staticmethod
    def plot_several(title, f, inst_list, plot_sense='frequencies',
                     show_bounds=True, outfile=None):
        """
        Plot PSD min_maxs for several Instruments
    
        All insts must be of same type (seismometer or pressure)
    
        Args:
            title (str): string at top of plot
            f (list or np.array): frequencies at which to plot
            inst_list (list): list of Instrument objects to plot
            plot_sense (str): 'frequencies': x-axis is frequencies
                              'periods': x-axis is periods
            show_bounds (bool): Show high/low noise reference levels
            outfile (str or None): write to the given file
        """
        S_low, S_high, model_label = _get_noise_model(f, inst_list)

        inst, color, label = inst_list[0]
        ax = inst.plot_minmax(f, color, show=False, plot_sense=plot_sense,
                              label=label)
        for x in inst_list[1:]:
            inst, color, label = x
            inst.plot_minmax(f, color, ax=ax, show=False, plot_sense=plot_sense,
                             label=label)
        if show_bounds is True:
            if plot_sense[0].lower() == 'p':
                ax.plot(np.power(f, -1.), S_low, 'k--', label=model_label)
                ax.plot(np.power(f, -1.), S_high, 'k--')
            else:
                ax.plot(f, S_low, 'k--', label=model_label)
                ax.plot(f, S_high, 'k--')
        ax.set_title(title)
        plt.legend()
        if outfile is not None:
            plt.savefig(outfile)
        plt.show()


def _get_noise_model(f, inst_list):
    # Figure out if sensors are pressure or velocity
    input_units = inst_list[0][0].input_units.lower()
    # Make sure that all sensors have the same input_units
    if len(inst_list) > 1:
        for inst in inst_list[1:]:
            assert inst[0].input_units.lower() == input_units
            
    if input_units.lower() == 'pa':
        S_low, S_high = Pressure_noise_model(f, as_freqs=True)
        return S_low, S_high, "Pressure Noise Model"
    if input_units.lower() == 'm/s':
        S_low, S_high = Peterson_noise_model(f, as_freqs=True)
        return S_low, S_high, "Peterson Noise Model"
    else:
        raise ValueError(f'Unknown {input_units=}')
    