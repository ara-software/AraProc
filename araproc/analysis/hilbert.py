# this calculates Hilbert envelope SNR for each channel and their average

import numpy as np
from scipy.signal import hilbert

from araproc.analysis import snr
from araproc.framework import waveform_utilities as wfu



def get_hilbert_snr(waveform):

    """
    Calculates Hilbert envelope SNR of a single voltage trace.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.

    Returns
    -------
    hill_snr: float
        The Hilbert envelope SNR of the waveform.
    """

    hill = wfu.get_hilbert_envelope(waveform)
    hill_max_idx = np.argmax(hill)
    hill_max = hill[hill_max_idx]
    hill_rms = snr.get_jackknife_rms(hill)

    if(hill_rms == 0.0):
      return 0

    hill_snr = hill_max/hill_rms

    return hill_snr




def get_avg_hilbert_snr(wave_bundle, excluded_channels = []):

    """
    Calculates channel-wise averaged Hilbert envelope SNR.

    Parameters
    ----------
    wave_bundle: dict of TGraphs
        Dictionary of waveform TGraphs to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    avg_hill_snr: float
        The average Hilbert envelope SNR.
    """

    chans = list(wave_bundle.keys())

    avg_hill_snr = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      hill_snr = get_hilbert_snr(wave_bundle[chan])
      avg_hill_snr.append(hill_snr)

    avg_hill_snr = np.mean(avg_hill_snr)

    return avg_hill_snr


def get_peak_over_avg_power(waveform, fraction=0.5):
    
    """
    Calculates the peak power at maximum of the Hilbert envelope divided
    by the average power in fraction of the full trace around the peak. 

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.

    Returns
    -------
    peak_over_rms: float
        The peak power over rms.
    """

    # calculate the power of the Hilbert envelope
    h = wfu.get_hilbert_envelope(waveform)
    envelope_power = np.square(np.abs(h))
    
    # find the peak and sort the power enevelope
    idx_max = np.argmax(envelope_power) 
    maxPower = envelope_power[idx_max]
    idx = np.arange(0, len(envelope_power), 1)
    closeness = abs(idx - idx_max)
    sorted_idx = np.argsort(closeness)
    
    envelope_power = envelope_power[sorted_idx]

    # get the portion of the trace within fraction of the peak
    fraction_from_peak = np.linspace(0, 1, len(envelope_power))
    mask = fraction_from_peak > fraction
    
    envelope_power = envelope_power[mask]

    # get average within fraction of the peak
    avg_power = np.mean(envelope_power)

    peak_over_rms = maxPower / avg_power

    return peak_over_rms

def get_moments_peak_over_avg_power(wave_bundle, excluded_channels=[]):
    
    """
    Calculates the peak power at maximum of the Hilbert envelope divided
    by the average power in fraction of the full trace around the peak for
    each channel, and then calculates the mean and std of this value over channels. 

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    mean_peak_over_rms: float
        Average peak power over rms.
    std_peak_over_rms: float
        Standard deviation of peak power over rms.
    """
    
    chans = list(wave_bundle.keys())

    poaps = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      waveform = wave_bundle[chan]
      poap = get_peak_over_avg_power(waveform)
      poaps.append(poap)
   
    mean_peak_over_rms = np.mean(poaps)
    std_peak_over_rms = np.std(poaps) 

    return mean_peak_over_rms, std_peak_over_rms

def get_rsd_peak_over_avg_power(wave_bundle, excluded_channels=[]):
    
    """
    Calculates the relative standard deviation of the POAP of the Hilbert enveloped waveform. 

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    rsd_peak_over_rms: float
        Relative standard deviation of peak power over rms.
    """
   
    mean_poap, std_poap = get_moments_peak_over_avg_power(wave_bundle, excluded_channels=excluded_channels)

    if np.isclose(mean_poap, 0):
        return 0. 

    rsd_poap = std_poap / mean_poap

    return rsd_poap



