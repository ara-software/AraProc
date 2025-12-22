import numpy as np
import ROOT

from araproc.framework import waveform_utilities as wfu

def get_number_saturated(waveform, saturation_threshold=1500):
    """
    Counts number of points in channel's trace which are saturated.

    Parameters
    ----------
    waveform: TGraphs or np.ndarrays
        Waveform TGraphs or np.ndarrays to be averaged.
    saturation_threshold : float
        Voltage beyond which a sample is considered saturated (in mV).

    Returns
    -------
    n_saturated : int
        Number of points which are above saturation threshold.
    """
        
    if(isinstance(waveform, ROOT.TGraph)):
      _, trace = wfu.tgraph_to_arrays(waveform) # in mV
    elif(isinstance(waveform, np.ndarray)):
      trace = np.copy(waveform) # in mV

    abs_trace = np.abs(trace)
    above_threshold = abs_trace > saturation_threshold
    n_saturated = above_threshold.sum()

    return n_saturated

def check_waveform_saturation(wave_bundle, excluded_channels=[]):

    """
    Checks whether any channel's trace is saturated.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to check.
    excluded_channels: list
        List of dictionary keys to exclude from check.

    Returns
    -------
    is_saturated : bool
        Whether a waveform is saturated.

    """
    
    chans = list(wave_bundle.keys())
    
    # set parameters for saturation check
    n_points = 5
    saturation_threshold = 1500 # mV

    # run check
    is_saturated = False
    for chan in chans:
        if chan in excluded_channels:
            continue 

        waveform = wave_bundle[chan]
        n_above_threshold = get_number_saturated(waveform)

        if(n_above_threshold > n_points):
            is_saturated = True
            break 

    return is_saturated

def check_valid_waveforms(wave_bundle, excluded_channels=[]):
    """
    Check the trace values are sensible.
    
    Parameters
    ----------
    wave_bundle : dict
        A dict with three entries:
          "event" : int  
            Event number
          "waveforms" : dict
            A dictionary of the 16 waveforms.
            The key is the RF channel number.
            The value is a TGraph.
            There should be 16 entries, even if you don't intend to use all
            16 traces in your interferometry.
            The exclusions are handled further down under the excluded channels section.
          "trace_type" : string
            Waveform type requested by which_trace
    excluded_channels: list
        List of dictionary keys to exclude from check.

    Returns
    -------
    is_invalid : bool
        Whether the waveform has invalid values. 
    """

    chans = list(wave_bundle.keys())

    # run check   
    is_invalid = False
    for chan in chans:
        if chan in excluded_channels:
            continue

        waveform = wave_bundle[chan]
        time, trace = wfu.tgraph_to_arrays(waveform)

        # check for nans
        if np.any(np.isnan(time)) or np.any(np.isnan(trace)):
            is_invalid = True
            break

        # check for flatlining 
        if trace.std() == 0.0:
            is_invalid = True
            break

    return is_invalid

def is_bad_waveform_quality(wavepacket, excluded_channels=[]):

    """
    Wrapper to run all waveform quality checks.

    Parameters
    ----------
    wavepacket : dict
        A dict with three entries:
          "event" : int  
            Event number
          "waveforms" : dict
            A dictionary of the 16 waveforms.
            The key is the RF channel number.
            The value is a TGraph.
            There should be 16 entries, even if you don't intend to use all
            16 traces in your interferometry.
            The exclusions are handled further down under the excluded channels section.
          "trace_type" : string
            Waveform type requested by which_trace
    excluded_channels: list
        List of dictionary keys to exclude from check.

    Returns
    -------
    is_bad_quality : bool
        Whether the waveform quality is bad.
    """

    is_bad_quality = False

    is_bad_quality = is_bad_quality or check_valid_waveforms(wavepacket["waveforms"], excluded_channels)

    return is_bad_quality


