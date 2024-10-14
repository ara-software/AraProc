import numpy as np
from scipy.ndimage import uniform_filter1d
from araproc.framework import waveform_utilities as wfu


def get_max_info(channel_wf, channel_time, use_debug=False):
    """
    Identifies the maximum voltage and its corresponding time.

    Parameters
    ----------
    channel_wf: numpy array
        Array of voltage values (squared waveform data).
    channel_time: numpy array
        Array of time values corresponding to the waveform data.
    use_debug: bool, optional
        If True, prints debug information.

    Returns
    -------
    max_bin: int
        Index of the maximum voltage in the waveform.
    max_val: float
        Maximum voltage value in the waveform.
    max_time: float
        Time corresponding to the maximum voltage value.
    """
    max_bin = np.nanargmax(channel_wf)
    max_val = channel_wf[max_bin]
    max_time = channel_time[max_bin]

    if use_debug:
        print(f"Max value: {max_val}, Max time: {max_time}")

    return max_bin, max_val, max_time

def get_sqrt_volt_sum_wf(channel_wf, dt):
    """
    Smooth the squared voltage waveform using a rolling mean filter.

    Parameters
    ----------
    channel_wf: numpy array
        Array of squared voltage values.
    dt: float
        Step size from time array.

    Returns
    -------
    smoothed_waveform: numpy array
        Smoothed version of the input waveform.
    """
    sum_win = 25
    sum_win_idx = int(np.round(sum_win / dt))

    return np.sqrt(uniform_filter1d(channel_wf, size=sum_win_idx, mode='constant'))

def get_mean_sigma_in_no_max(channel_wf, max_bin, wf_len, use_debug=False):
    """
    Computes the mean and standard deviation of the waveform, excluding the region
    around the maximum voltage.

    Parameters
    ----------
    channel_wf: numpy array
        Smoothed voltage values of the waveform.
    max_bin: int
        Index of the maximum voltage in the waveform.
    wf_len: int
        Total number of samples in the waveform.
    use_debug: bool, optional
        If True, prints debug information.

    Returns
    -------
    noise_mean: float
        Mean voltage excluding the peak window.
    noise_sigma: float
        Standard deviation excluding the peak window.
    """
    bin_8 = wf_len // 8
    pad_no_max = channel_wf.copy()
    front_idx = max(0, int(max_bin - bin_8))

    pad_no_max[front_idx:max_bin + bin_8 + 1] = np.nan
    noise_mean = np.nanmean(pad_no_max)
    noise_sigma = np.nanstd(pad_no_max)

    if use_debug:
        print(f"Mean: {noise_mean}, Sigma: {noise_sigma}")

    return noise_mean, noise_sigma

def get_ch_sliding_v2_snr_uw(channel_wf, channel_time, wf_len, use_debug=False):
    """
    Computes the RPR value, which is similar to SNR, for the given waveform.

    Parameters
    ----------
    channel_wf: numpy array
        Voltage values (squared waveform data).
    channel_time: numpy array
        Time values corresponding to the waveform.
    wf_len: int
        Total number of samples in the waveform.
    use_debug: bool, optional
        If True, prints debug information.

    Returns
    -------
    rpr_val: float
        The RPR value (signal-to-noise-like ratio) for the waveform.
    """
    dt =  wfu.get_dt_and_sampling_rate(channel_time)[0] 
    channel_wf = get_sqrt_volt_sum_wf(channel_wf, dt)
    max_bin, max_val, max_time = get_max_info(channel_wf, channel_time, use_debug)
    noise_mean, noise_sigma = get_mean_sigma_in_no_max(channel_wf, max_bin, wf_len, use_debug)
    rpr_val = (max_val - noise_mean) / noise_sigma

    if noise_sigma <= 0:
        print(f'Negative std detected, which may be unphysical')
        rpr_val = -9999

    if use_debug:
        print(f"RPR: {rpr_val}, Max Val: {max_val}, Mean: {noise_mean}, Sigma: {noise_sigma}")

    return rpr_val

def run_rpr_calculation(waveform, use_debug=False):
    """
    Computes the RPR value for a given waveform by preparing the data and invoking the RPR calculation.

    Parameters
    ----------
    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.
    snr: float, optional
        Placeholder for potential SNR calculation (currently unused).
    use_debug: bool, optional
        If True, prints debug information.

    Returns
    -------
    rpr_val: float
        The calculated RPR value for the waveform.
    """
    channel_time, channel_wf = wfu.tgraph_to_arrays(waveform) 
    wf_len = len(channel_wf)
    channel_wf = channel_wf ** 2
    channel_wf[np.isnan(channel_wf)] = 0
    rpr_val = get_ch_sliding_v2_snr_uw(channel_wf, channel_time, wf_len, use_debug)
    return rpr_val

def get_avg_rpr(wave_bundle, chans=None, individual_antenna=False, use_debug=False):
    """
    Calculates channel-wise averaged RPR for a given set of waveforms.

    Parameters
    ----------
    wave_bundle: dict of tuples
        Dictionary containing tuples of (voltage array, time array) for each channel.
    chans: list, optional
        List of channels to average over. If None, averages over all channels in wave_bundle.
    individual_antenna: bool, optional
        If True, returns the RPR value for each channel individually along with the average.
    use_debug: bool, optional
        If True, prints debug information.

    Returns
    -------
    rpr_collect: list of floats, optional
        List of RPR values for each individual channel (if individual_antenna is True).
    avg_rpr: float
        The average RPR across the selected channels.
    """
    if chans is None:
        chans = list(wave_bundle.keys())

    avg_rpr = []
    if individual_antenna:
        rpr_collect = []

    for chan in chans:
        waveform = wave_bundle[chan]  # Unpack voltage and time arrays from the wave_bundle
        rpr = run_rpr_calculation(waveform, use_debug=use_debug)  # Calculate RPR for the channel
        avg_rpr.append(rpr)
        
        if individual_antenna:
            rpr_collect.append(rpr)
    
    avg_rpr = np.mean(avg_rpr)
    
    if use_debug:
        print('Average RPR: ', avg_rpr)

    if individual_antenna:
        return rpr_collect, avg_rpr
    else:
        return avg_rpr





