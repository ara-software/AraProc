import logging
import numpy as np
import numpy.typing as npt
import ROOT
from scipy import interpolate

import importlib.resources as pkg_resources
from . import data

def load_arasim_phase_response_as_spline():

    """
    Load the AraSim phase response of the system to use for de-dispersion.

    Returns
    -------
    the_phase_spline : interp1d
        A spline of the unwrapped phases as a function of frequency.
        The frequency units are Ghz.
        And the phase is unwrapped, in units of radians.
    """

    file = pkg_resources.open_text(data, 
                                   "ARA_Electronics_TotalGain_TwoFilters.txt")
    file_content = np.genfromtxt(file, 
                                 delimiter=",", skip_header=3,
                                 names=["freq", "gain", "phase"], 
                                )

    freq_ghz = file_content["freq"]/1.E3 # convert to GHz
    phs_unwrapped = np.unwrap(file_content["phase"]) # unwrapped phase in radians

    the_phase_spline = interpolate.interp1d(
        freq_ghz, phs_unwrapped,
        bounds_error=False,
        fill_value=0
    )
    file.close()

    return the_phase_spline

def dedisperse_wave(
        times, # in nanoseconds,
        volts, # in volts,
        phase_spline # the  phase spline
        ):
    
    """
    Fetch a specific calibrated event

    Parameters
    ----------
    times : np.ndarray(dtype=np.float64)
        A numpy array of floats containing the times for the trace,
        in nanoseconds.
    volts : np.ndarray(dtype=np.float64)
        A numpy array of floats containing the voltages for the trace,
        in volts.
    phase_spline : interp1d
        A spline of the unwrapped phase (in radians) vs frequency (in GHz).
        When the function was first written, it was meant to utilize
        the output of `dedisperse.load_arasim_phase_response_as_spline`.
        So check that function for an example of how to do it.

    Returns
    -------
    dedispersed_wave : np.ndarray(dtype=np.float64)
        The dedispersed wave
        
    """

    if len(times) != len(volts):
        raise Exception("The time and volts arrays are mismatched in length. Abort.")

    

