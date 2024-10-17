# this calculates Hilbert envelope SNR for each channel and their average

import numpy as np
from scipy.signal import hilbert

from araproc.analysis import snr
from araproc.framework import waveform_utilities as wfu



def get_envelope(waveform, nsegs = 8):

    """
    Calculates the Hilbert envelope voltage and rms of a trace.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.
    nsegs: int
        Number of segments to break waveform into.

    Returns
    -------
    hill_max: float
        Hilbert envelope maximum voltage in same units as trace.
    hill_rms: float
        Hilbert envelope minimum rms among all segments in same units as trace..
    """ 

    _, trace = wfu.tgraph_to_arrays(waveform)
   
    hill = np.abs(hilbert(trace))
    hill_max_idx = np.nanargmax(hill)
    hill_max = hill[hill_max_idx]
  
    ''' 
    segments = np.array_split(hill, nsegs)
    rms_segments = [np.sqrt(np.nanmean(segment**2)) for segment in segments]
    sorted_rms_segments = sorted(rms_segments)
    hill_rms = np.mean(sorted_rms_segments[:3])    
    '''

    hill_rms = snr.get_min_segmented_rms(hill)
 
    return hill_max, hill_rms




def get_hill_snr(waveform):

    """
    Calculates Hilbert envelope SNR of a single voltage trace.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.

    Returns
    -------
    hill_snr : float
        The Hilbert envelope SNR of the waveform.
    """

    hill_max, hill_rms = get_envelope(waveform)
   
    if(hill_rms == 0.0):
      return 0

    hill_snr = hill_max/hill_rms

    return hill_snr




def get_avg_hill_snr(wave_bundle, excluded_channels = []):

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
    avg_hill_snr : float
        The average Hilbert envelope SNR.
    """

    chans = list(wave_bundle.keys())

    avg_hill_snr = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      hill_snr = get_hill_snr(wave_bundle[chan])
      avg_hill_snr.append(hill_snr)

    avg_hill_snr = np.mean(avg_hill_snr)

    return avg_hill_snr

