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

