import numpy as np
from scipy.ndimage import uniform_filter1d
from araproc.framework import waveform_utilities as wfu
from araproc.analysis import snr

def get_rpr(waveform):

    """
    Computes the RPR (Root Power ratio) value, similar to SNR, for the given waveform.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform containing time and voltage values.

    Returns
    -------
    rpr_val: float
        The RPR value for the waveform, calculated as the ratio of the maximum voltage 
        (after smoothing) to noise RMS of the waveform.
    """

    # Extract time and voltage arrays from the waveform
    channel_time, channel_wf = wfu.tgraph_to_arrays(waveform)
    wf_len = len(channel_wf)
    
    # Square the waveform data for further processing
    channel_wf = channel_wf ** 2

    # Calculate the smoothing window size based on sampling rate
    dt =  wfu.get_dt_and_sampling_rate(channel_time)[0]
    sum_win = 25  # Smoothing window in ns
    sum_win_idx = int(np.round(sum_win / dt))  # Convert window size to sample points
    channel_wf = np.sqrt(uniform_filter1d(channel_wf, size=sum_win_idx, mode='constant'))

    # Find the maximum value of the smoothed waveform
    max_bin = np.argmax(channel_wf)
    max_val = channel_wf[max_bin]

    # Get noise rms from snr module
    noise_sigma = snr.get_min_segmented_rms(channel_wf)

    # Calculate and return the RPR value
    rpr_val = max_val / noise_sigma

    # Calculate the smoothing window size based on sampling rate
    dt =  wfu.get_dt_and_sampling_rate(channel_time)[0]
    sum_win = 25  # Smoothing window in ns
    sum_win_idx = int(np.round(sum_win / dt))  # Convert window size to sample points
    channel_wf = np.sqrt(uniform_filter1d(channel_wf, size=sum_win_idx, mode='constant'))

    # Find the maximum value of the smoothed waveform
    max_bin = np.argmax(channel_wf)
    max_val = channel_wf[max_bin]

    # Read noise rms from snr module
    noise_sigma = snr.get_min_segmented_rms(channel_wf)
    # Calculate and return the RPR value
    rpr_val = max_val / noise_sigma

    return rpr_val

def get_avg_rpr(wave_bundle, chans=None, excluded_channels=[]):
    """
    Calculates the average RPR across selected channels from a given set of waveforms.

    Parameters
    ----------
    wave_bundle: dict of tuples
        Dictionary containing tuples of (voltage array, time array) for each channel.
    chans: list, optional
        List of channels to calculate the average RPR for. If None, averages over all channels in wave_bundle.
    excluded_channels: list, optional
        List of channels to exclude from the calculation.

    Returns
    -------
    avg_rpr: float
        The average RPR value across the selected channels.
    """
    chans = list(wave_bundle.keys()) if chans is None else chans
    avg_rpr = []

    for chan in chans:
        if chan in excluded_channels:
            continue
        waveform = wave_bundle[chan]  # Unpack voltage and time arrays from the wave_bundle
        rpr = get_rpr(waveform)  # Calculate RPR for the channel
        avg_rpr.append(rpr)

    # Return the average RPR value across all selected channels
    return np.mean(avg_rpr)

