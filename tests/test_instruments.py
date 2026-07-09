# test_instruments.py

import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from crawtools import instruments


# ----------------------------------------------------------------------
# Sensor
# ----------------------------------------------------------------------

def test_sensor_requires_valid_input_units():
    with pytest.raises(AssertionError):
        instruments.Sensor(
            gain=1000,
            input_units="invalid"
        )


def test_sensor_stores_attributes():
    sensor = instruments.Sensor(
        gain=1000,
        input_units="m/s",
        comments="test sensor"
    )
    assert sensor.input_units == "m/s"
    assert sensor.comments == "test sensor"
    assert sensor.self_noise is None
    assert sensor.clip_level is None


def test_sensor_validates_self_noise_structure():
    with pytest.raises(AssertionError):
        instruments.Sensor(
            gain=1000,
            self_noise=[[1.0]]
        )


def test_sensor_accepts_valid_self_noise_structure():
    sensor = instruments.Sensor(
        gain=1000,
        self_noise=[
            [1.0, -120.0],
            [10.0, -130.0],
        ]
    )

    assert len(sensor.self_noise) == 2


def test_sensor_validates_clip_level_structure():
    with pytest.raises(AssertionError):
        instruments.Sensor(
            gain=1000,
            clip_level=[[1.0]]
        )


# ----------------------------------------------------------------------
# Preamplifier
# ----------------------------------------------------------------------

def test_preamplifier_creates_response():
    preamp = instruments.Preamplifier(gain=20)

    assert preamp.response.stage_sequence_number == 2
    assert preamp.response.input_units == "V"
    assert preamp.response.output_units == "V"


# ----------------------------------------------------------------------
# Logger
# ----------------------------------------------------------------------

def test_logger_gain_and_dynamic_range():
    logger = instruments.Logger(fs_counts=1_000_000,
                                fs_volts=10,
                                dynamic_range_dB=120)
    expected_gain = 1000000 / 10
    assert logger.response.stage_gain == expected_gain

    expected_min = 1_000_000 / (10 ** (120 / 20))
    assert np.isclose(logger.min_counts, expected_min)

    assert logger.max_counts == 1_000_000


# ----------------------------------------------------------------------
# Instrument
# ----------------------------------------------------------------------

@pytest.fixture
def instrument_obj():
    sensor = instruments.Sensor(gain=1000)
    preamp = instruments.Preamplifier(gain=10)
    logger = instruments.Logger(fs_counts=1_000_000,
                                fs_volts=10,
                                dynamic_range_dB=120)
    return instruments.Instrument(sensor, preamp, logger)


def test_instrument_input_units_property(instrument_obj):
    assert instrument_obj.input_units == "m/s"


# def test_instrument_calculates_min_values():
#     frequencies = np.array([1.0, 10.0])
# 
#     inst = MagicMock()
#     inst.logger.min_counts = 100
#     inst.evalresp.return_value = np.array([10.0, 20.0])
# 
#     result = instruments.Instrument.min_values(inst, frequencies)
# 
#     np.testing.assert_allclose(result, [10.0, 5.0])
# 
# 
# def test_instrument_calculates_max_values():
#     frequencies = np.array([1.0, 10.0])
# 
#     inst = MagicMock()
#     inst.logger.max_counts = 1000
#     inst.evalresp.return_value = np.array([10.0, 20.0])
# 
#     result = instruments.Instrument.max_values(inst, frequencies)
# 
#     np.testing.assert_allclose(result, [100.0, 50.0])


# def test_evalresp_uses_acc_for_seismic_sensor(instrument_obj):
#     frequencies = np.array([1.0])
# 
#     instrument_obj.response = MagicMock()
#     instrument_obj.response.get_evalresp_response_for_frequencies.return_value = (
#         np.array([123.0])
#     )
# 
#     instrument_obj.evalresp(frequencies)
# 
#     instrument_obj.response.get_evalresp_response_for_frequencies.assert_called_once_with(
#         frequencies,
#         "ACC",
#     )


def test_evalresp_uses_def_for_pressure_sensor():
    sensor = instruments.Sensor(
        gain=1000,
        input_units="Pa"
    )

    preamp = instruments.Preamplifier(gain=10)

    logger = instruments.Logger(
        fs_counts=1000,
        fs_volts=10,
        dynamic_range_dB=80,
    )

    inst = instruments.Instrument(sensor, preamp, logger)

    inst.response = MagicMock()
    inst.response.get_evalresp_response_for_frequencies.return_value = np.array([1])

    inst.evalresp(np.array([1.0]))

    inst.response.get_evalresp_response_for_frequencies.assert_called_once_with(
        np.array([1.0]),
        "DEF",
    )

def test_instrument_max_level():
    sensor_gain = 3  # m/s / volt
    preamp_gain = 5  # V/V
    logger_volt_range = 7  # V / full range
    logger_count_range = 11  # counts / full range
    logger_dynamic_range_dB = 13
    inst = instruments.Instrument(
        instruments.Sensor(gain=sensor_gain),
        instruments.Preamplifier(gain=preamp_gain),
        instruments.Logger(fs_counts=logger_count_range,
                           fs_volts=logger_volt_range,
                           dynamic_range_dB=logger_dynamic_range_dB))
    f = np.array([1/(2*np.pi)])
    print(inst)
    print(f'{np.abs(inst.max_values(f))=}')
    print(f'{np.abs(inst.evalresp(f))=}')
    # inst.plot(min_freq=0.01)
    assert np.abs(inst.max_values(f)) == np.array([logger_volt_range/(sensor_gain*preamp_gain)])
    
    

# ----------------------------------------------------------------------
# plot_minmax
# ----------------------------------------------------------------------

@patch("matplotlib.pyplot.show")
def test_plot_minmax_returns_axis(mock_show, instrument_obj):
    frequencies = np.array([1.0, 10.0])

    instrument_obj.min_values = MagicMock(
        return_value=np.array([1.0, 1.0])
    )
    instrument_obj.max_values = MagicMock(
        return_value=np.array([10.0, 10.0])
    )

    ax = instrument_obj.plot_minmax(
        frequencies,
        show=False
    )

    assert ax is not None
    mock_show.assert_not_called()


def test_plot_minmax_invalid_plot_sense(instrument_obj):
    frequencies = np.array([1.0, 10.0])

    with pytest.raises(ValueError):
        instrument_obj.plot_minmax(
            frequencies,
            plot_sense="invalid",
            show=False,
        )


# ----------------------------------------------------------------------
# _get_noise_model
# ----------------------------------------------------------------------

# @patch("instruments.Peterson_noise_model")
# def test_get_noise_model_seismic(mock_model):
#     mock_model.return_value = (
#         np.array([-150]),
#         np.array([-100]),
#     )
# 
#     sensor = instruments.Sensor(gain=1000, input_units="m/s")
#     preamp = instruments.Preamplifier(gain=1)
#     logger = instruments.Logger(1000, 10, 80)
# 
#     inst = instruments.Instrument(sensor, preamp, logger)
# 
#     low, high, label = instruments._get_noise_model(
#         np.array([1.0]),
#         [(inst, "b", "test")]
#     )
# 
#     assert label == "Peterson Noise Model"
#     np.testing.assert_array_equal(low, np.array([-150]))
#     np.testing.assert_array_equal(high, np.array([-100]))
# 
# 
# @patch("instruments.Pressure_noise_model")
# def test_get_noise_model_pressure(mock_model):
#     mock_model.return_value = (
#         np.array([-80]),
#         np.array([-50]),
#     )
# 
#     sensor = instruments.Sensor(gain=1000, input_units="Pa")
#     preamp = instruments.Preamplifier(gain=1)
#     logger = instruments.Logger(1000, 10, 80)
# 
#     inst = instruments.Instrument(sensor, preamp, logger)
# 
#     low, high, label = instruments._get_noise_model(
#         np.array([1.0]),
#         [(inst, "b", "test")]
#     )
# 
#     assert label == "Pressure Noise Model"


def test_get_noise_model_rejects_mixed_units():
    sensor1 = instruments.Sensor(gain=1000, input_units="m/s")
    sensor2 = instruments.Sensor(gain=1000, input_units="Pa")

    preamp = instruments.Preamplifier(gain=1)
    logger = instruments.Logger(1000, 10, 80)

    inst1 = instruments.Instrument(sensor1, preamp, logger)
    inst2 = instruments.Instrument(sensor2, preamp, logger)

    with pytest.raises(AssertionError):
        instruments._get_noise_model(
            np.array([1.0]),
            [
                (inst1, "b", "seismic"),
                (inst2, "r", "pressure"),
            ]
        )