# this looks for and identifies the spark events by calculating power ratio
# and comparing with a set threshold 

import numpy as np

from araproc.framework import waveform_utilities as wfu





def get_rms_power(waveform):
  
    """
    Calculates RMS power of a waveform.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.

    Returns
    -------
    power: float
        RMS power of the waveform.
    """

    _, trace = wfu.tgraph_to_arrays(waveform)
    power = np.nanmean(trace**2)

    return power


def threshold(threshold_ratio = 5.0):
    
    """
    Defines the spark event power ratio; fixed '5.0' for all stations, can be overwritten.

    Parameters
    ----------
    threshold_ratio: int
       Spark event power ratio.

    Returns
    -------
    cut_val: float
       Spark event power ratio.
    """
        
    cut_val = threshold_ratio

    return cut_val


def get_power_ratio(wave_bundle, excluded_channels = [], cut_ratio = None):
    
    """
    Calculates power ratio between largest and second largest rms power strings.

    Parameters
    ----------
    wave_bundle: dict of TGraphs
        Dictionary of waveform TGraphs to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.
    cut_ratio: float
        Power ratio threshold.

    Returns
    -------
    power_ratio: float
        Power ratio of largest to second largest rms power strings.
    is_spark_event: bool
        0: is not a spark event, 1: is a spark event.
    """

    if cut_ratio is not None:
       cut_ratio = float(cut_ratio)
    else:
       cut_ratio = threshold()


    num_strings = 4

    chans = list(wave_bundle.keys())
    string_chans = [chans[i::4] for i in range(num_strings)]

    power = []
    for ch in np.array(string_chans).flatten():
       if ch in excluded_channels:
          power.append(np.nan)
       else:
          power.append(get_rms_power(wave_bundle[ch])) 
    
    power = np.asarray(power)
    power_avg = np.nanmean(np.array_split(power, num_strings), axis = 1)
    power_avg_sort = -np.sort(-power_avg)
    power_ratio = power_avg_sort[0]/power_avg_sort[1]
    
    is_spark_event = int(power_ratio > cut_ratio)
    
    return power_ratio, is_spark_event
