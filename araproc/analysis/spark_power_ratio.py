# this calculates power ratio of the largest to second largest mean power strings

import warnings
import numpy as np
from collections import defaultdict

from araproc.framework import constants as const
from araproc.framework import waveform_utilities as wfu


def get_mean_power(waveform):
  
    """
    Calculates mean power (mean of square of voltage) of a waveform.

    Parameters
    ----------
    waveform: TGraph
        A TGraph of the waveform.

    Returns
    -------
    power: float
        Mean power of the waveform.
    """

    _, trace = wfu.tgraph_to_arrays(waveform)
    power = np.mean(trace**2)

    return power


def get_power_ratio(wave_bundle, excluded_channels = []):
    
    """
    Calculates power ratio between largest and second largest mean power strings.

    Parameters
    ----------
    wave_bundle: dict of TGraphs
        Dictionary of waveform TGraphs to be averaged.
    excluded_channels: list
    
    Returns
    -------
    power_ratio: float
        Power ratio of largest to second largest mean power strings.
    """

    chans = sorted(list(wave_bundle.keys()))

    power = defaultdict(list)
    for chan in chans:
       if chan in excluded_channels:
          continue

       thisPower = get_mean_power(wave_bundle[chan])
       strNo = chan % const.num_string
       power[strNo].append(thisPower)
    
    power_avg = []
    for strNo in power:   
       power_avg.append(np.mean(power[strNo]))

    power_avg = np.asarray(power_avg)
    power_avg_sort = -np.sort(-power_avg)
    power_ratio = power_avg_sort[0]/power_avg_sort[1]
    
    return power_ratio
