"""
Functions to convert rotation to horizotal seismometer signals, and vice versa,
assuming only rotation through the gravitational field

Currently uses small-angle theory, simply multiplying radians by the gravitational attraction.


Functions are:
- `rotation_to_motion` convert a rotation signal to a motion signal
- `motion_to_rotation` inverse of the above
- `convert_rotation`: convert rotation from one units to another
- `convert_motion`: convert motion from one units to another
- `to_rads`: convert rotation units to radians
- `from_rads`: convert rotation units from radians
"""
import numpy as np


def rotation_to_motion(periods, vals, in_units, out_units="m/s^2"):
    """
    Return the horizontal PSD that will result from rotation through grav field

    PSD is assumed to be in units/sqrt(Hz)

    Args:
        periods (np.array): periods in seconds
        vals (np.array or float): rotationalPSDvalue at each period.  If a scalar, the
            same value is applied at every period.
        in_units (str): Reference units for rotation: "rad", "rad.s", "rad/s",
            "nrad", "nrad.s", "nrad/s"
        out_units (str): Reference units for horizontal sensor: 'm', 'm/s' or 'm/s^2'
    """
    vals = _match_dims(vals, periods)
    vals = _to_rads(periods, vals, in_units)
    outp = vals * 9.8
    outp = _from_accel(periods, outp, out_units)
    return outp


def motion_to_rotation(periods, vals, in_units='m/s^2', out_units='rad'):
    """
    Convert horizontal motion to rotational equivalent

    PSD is assumed to be in in_units/sqrt(Hz)

    Args:
        periods (np.array): periods in seconds
        vals (np.array or float):  PSD value in (in_units)/sqrt(Hz) at each period.
            If a scalar, the same value is applied at every period.
        out_units (str): Reference units for rotation: '(n)rad', '(n)rad/s', '(n)rad.s'
    """
    vals = _match_dims(vals, periods)
    vals = _to_accel(periods, vals, in_units)
    outp = vals / 9.8
    outp = _from_rads(periods, outp, out_units)
    return outp
    
def convert_rotation(periods, vals, in_units, out_units):
    """Convert rotation units"""
    vals = _match_dims(vals, periods)
    vals = _to_rads(periods, vals, in_units)
    return _from_rads(periods, vals, out_units)

def convert_motion(periods, vals, in_units, out_units):
    """Convert rotation units"""
    vals = match_dims(vals, periods)
    vals = _to_accel(periods, vals, in_units)
    return _from_accel(periods, vals, out_units)


def _match_dims(vals, array):
    """
    Verify dimensions of vals against array
    
    If vals is a scalar, convert to constant array of same size as array
    """
    if isinstance(vals, (int, float)):
        vals = vals*np.ones(array.shape)
    assert array.shape == vals.shape
    return vals


def _to_rads(periods, vals, in_units):
    # Convert rotation input units to rad
    match in_units:
        case "rad":
            return vals
        case "nrad":
            return vals/1.e9
        case "nrad/s":
            return periods*vals/(2*np.pi*1.e9)
        case "nrad.s":
            return (2*np.pi)*vals/(1.e9*periods)
        case "rad/s":
            return periods*vals/(2*np.pi)
        case "rad.s":
            return (2*np.pi)*vals/periods
        case _:
            raise ValueError(f'{in_units=} not in ("rad", "rad/s", "rad.s", "nrad", "nrad.s", "nrad/s")')

def _from_rads(periods, vals, out_units):
    # Convert rotation input units from rads to another unit
    match out_units:
        case "rad":
            return vals
        case "nrad":
            return vals*1.e9
        case "rad/s":
            return (2*np.pi)*vals/periods
        case "nrad/s":
            return (2*np.pi*1.e9)*vals/periods
        case "rad.s":
            return periods*vals/(2*np.pi)
        case "nrad.s":
            return (1.e9*periods)*vals/(2*np.pi)
        case _:
            raise ValueError(f'{out_units=} not in ("rad", "rad/s", "rad.s", "nrad", "nrad.s", "nrad/s")')

def _to_accel(periods, vals, in_units):
    """ Convert displacement input units to m/s^2"""
    match in_units:
        case "m/s^2":
            return vals
        case "m/s":
            return (2*np.pi) * vals / periods
        case "m":
            return (2*np.pi)**2 * vals / periods
        case _:
            raise ValueError(f'{in_units=} not in ("m", "m/s", "m/s^2")')


def _from_accel(periods, vals, out_units):
    """ Convert m/s^2 to other displacement units"""
    match out_units:
        case "m/s^2":
            return vals
        case "m/s":
            return periods*vals/(2*np.pi)
        case "m":
            return vals*(periods/(2*np.pi))**2
        case _:
            raise ValueError(f'{out_units=} not in ("m", "m/s", "m/s^2")')


