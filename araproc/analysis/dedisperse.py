import numpy as np
from scipy import interpolate

import importlib.resources as pkg_resources
from . import data

from araproc.framework import waveform_utilities as wu

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

    the_phase_spline = interpolate.Akima1DInterpolator(
        freq_ghz, phs_unwrapped,
        method="makima",
    )
    # turn off extrapolation outside the region of support
    the_phase_spline.extrapolate = False
    file.close()

    return the_phase_spline

def eval_splined_phases(phase_spline, freqs_to_evaluate):
    """"
    Just a little helper function.
    This is necessary because the Akima Interpolator will return NaN 
    when called out of the range of support, but we'd rather it gave zeros.
    """
    these_phases = phase_spline(freqs_to_evaluate)
    these_phases = np.nan_to_num(these_phases) # convert nans to zeros
    return these_phases

def dedisperse_wave(
        times, # in nanoseconds,
        volts, # in volts,
        phase_spline, # the phase spline
        pre_pad_ns=40.0,  # amount of zero padding before dedispersion (ns)
        final_crop_ns=10.0  # amount of padding to keep after dedispersion (ns)
        ):
    
    """
    Fetch a specific calibrated event with configurable padding

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
    pre_pad_ns : float, optional
        Amount of zero padding to add before dedispersion, in ns
    final_crop_ns : float, optional
        Amount of the pre-padding to keep in final result, in ns.
        The default value of 10 ns is selected by looking at the 
        "average" group delay (by eye) within the band.
        

    Returns
    -------
    times : np.ndarray(dtype=np.float64)
        The time array for the dedispersed wave
    dedispersed_wave : np.ndarray(dtype=np.float64)
        The dedispersed wave
        
    """

    if len(times) != len(volts):
        raise Exception("The time and volts arrays are mismatched in length. Abort.")

    # Calculate time step from the input times array
    dt = times[1] - times[0]  # assuming uniform sampling
    
    # Calculate number of padding samples needed
    n_pad_samples = int(pre_pad_ns / dt)
    
    # Create padded time array: from (t_i - pre_pad_ns) to t_f
    t_start_padded = times[0] - pre_pad_ns
    padded_times = np.arange(t_start_padded, times[-1] + dt/2, dt)
    
    # Create padded voltage array with zeros at the beginning
    padded_volts = np.zeros(len(padded_times))
    padded_volts[n_pad_samples:n_pad_samples+len(volts)] = volts
    
    # Get the frequency domain representation of the padded trace
    freqs, spectrum = wu.time2freq(padded_times, padded_volts)

    # interpolate the *unwrapped phases* to the correct frequency base
    phased_interpolated = eval_splined_phases(phase_spline, freqs)
    
    # convert these into a complex number
    phased_rewrapped = np.exp((0 + 1j)*phased_interpolated)
    
    # do complex division to do the dedispersion
    spectrum /= phased_rewrapped

    # back to the time domain
    dedispersed_times, dedispersed_volts = wu.freq2time(padded_times, spectrum)
    
    # Crop the result: keep from (t_i - final_crop_ns) to t_f
    t_start_final = times[0] - final_crop_ns
    
    # Find the indices for cropping
    crop_start_idx = np.argmin(np.abs(dedispersed_times - t_start_final))
    original_end_idx = np.argmin(np.abs(dedispersed_times - times[-1]))
    
    # Crop to final desired range
    crop_times = dedispersed_times[crop_start_idx:original_end_idx+1]
    crop_volts = dedispersed_volts[crop_start_idx:original_end_idx+1]
    
    return crop_times, crop_volts
