# this calculates correlation snr of pair of channels (ch) separately 
# for VPols (RF ch 0-7) and HPols (RF ch 8-15)

import scipy
import itertools
import numpy as np
from scipy.signal import correlate as co

from araproc.analysis import snr
from araproc.framework import waveform_utilities as wfu


def channel_pairs():
    
    """
    Makes the lists of rf channel pairs, and, VPol and HPol pairs.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    pairs: list
       List of rf channel pairs.
    paits_V: list
       List of VPol rf channel pairs.
    paits_H: list
       List of HPol rf channel pairs.
    """

    rf_chs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    pairs = list(itertools.combinations(rf_chs, 2))
    pairs_V = list(itertools.combinations(rf_chs[:8], 2))
    pairs_H = list(itertools.combinations(rf_chs[8:], 2))

    return pairs, pairs_V, pairs_H



def pairing(channels):
    
    """
    Makes channel pairs.
    
    Parameters
    ----------
    channels: list
       List of given channels.

    Returns
    -------
    pairs: list
       List of all possible pair combinations of given channels.
    """
 
    pairs = list(itertools.combinations(channels, 2))
    
    return pairs



def get_corr(y1, y2):
    
    """
    Calculates normalized correlation between two waveforms, y1 and y2.
    
    Parameters
    ----------
    y1: list
       List of values of waveform 1.
    y2: list
       List of values of waveform 2.
    
    Returns
    -------
    corr: array
       Array of values of cross correlation function.  
    """

    _, y1 = wfu.tgraph_to_arrays(y1)
    _, y2 = wfu.tgraph_to_arrays(y2)
    
    n = len(y1)
    corr = co(y2, y1, mode = 'same')/np.sqrt(co(y1, y1, mode = 'same')[int(n/2)]*co(y2, y2, mode = 'same')[int(n/2)])    

    return corr



def get_corr_snr(first_ch, second_ch):

    """
    Calculates channel-pair correlation SNR
    
    Parameters
    ----------
    first_ch: array 
       Array of a given channel waveform.
    second_ch: array 
       Array of another given channel waveform.
     
    Returns
    -------
    corr_snr: float
       Channel pair correlation SNR.
    """
    
    corr = get_corr(first_ch, second_ch)
    corr_snr = snr.get_snr(corr)

    return corr_snr


def get_avg_corr_snr(wave_bundle, excluded_channels = []):

    """
    Makes useful dictionaries for channel-pair correlation SNR and their average.

    Parameters
    ----------
    wave_bundle: dict of TGraphs
        Dictionary of waveform TGraphs to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    avg_corr_snr : float
        Average channel-pair correlation SNR.
    avg_vpol_corr_snr : float
        The VPol average channel-pair correlation SNR.
    avg_hpol_corr_snr : float
        The HPol average channel-pair correlation SNR.
    """


    chans = list(wave_bundle.keys())

    AllPol_dict, VPol_dict, HPol_dict = {}, {}, {}
    corr_snr, VPol_corr_snr, HPol_corr_snr = {}, {}, {}

    for chan in chans:

      if(chan in excluded_channels):
        continue

      AllPol_dict.update({chan: wave_bundle[chan]})

      if chan in chans[:8]:
        VPol_dict.update({chan: wave_bundle[chan]})
      if chan in chans[8:]:
        HPol_dict.update({chan: wave_bundle[chan]})
    
    avg_corr_snr, avg_vpol_corr_snr, avg_hpol_corr_snr = [], [], []

    pairs = pairing(list(AllPol_dict.keys()))
    pairs_V = pairing(list(VPol_dict.keys()))
    pairs_H = pairing(list(HPol_dict.keys()))

    for pair in pairs:
        corr_snr = get_corr_snr(AllPol_dict[pair[0]], AllPol_dict[pair[1]])
        avg_corr_snr.append(corr_snr)

    for pair in pairs_V:
        corr_snr = get_corr_snr(VPol_dict[pair[0]], VPol_dict[pair[1]])
        avg_vpol_corr_snr.append(corr_snr)
  
    for pair in pairs_H:
        corr_snr = get_corr_snr(HPol_dict[pair[0]], HPol_dict[pair[1]])
        avg_hpol_corr_snr.append(corr_snr)

    avg_corr_snr = np.mean(avg_corr_snr)
    avg_vpol_corr_snr = np.mean(avg_vpol_corr_snr)
    avg_hpol_corr_snr = np.mean(avg_hpol_corr_snr)

    return avg_corr_snr, avg_vpol_corr_snr, avg_hpol_corr_snr



def calculate_corr_snr_dict(d):
    
    """
    Calculates channel-pair correlation SNR
    
    Parameters
    ----------
    d: dict
       Dictionary of given wave bundles.
     
    Returns
    -------
    corr_snr_dict: dict
       Dictionary of channel pair SNR.
    """
    
    chans = list(d.keys())
    pairs = pairing(chans)

    corr_snr_dict = {}
    for pair in pairs:
        corr = get_corr(d[pair[0]], d[pair[1]])
        corr_snr = snr.get_snr(corr)
        corr_snr_dict.update({pair:corr_snr})
            
    return corr_snr_dict    


def get_corr_snr_dict(wave_bundle, excluded_channels = []):

    """
    Makes useful dictionaries for channel-pair correlation SNR and their average.

    Parameters
    ----------
    wave_bundle: dict of TGraphs
        Dictionary of waveform TGraphs to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    corr_snr : dict
        The channel-pair correlation SNR, their average.
    VPol_corr_snr : dict
        The VPol channel-pair correlation SNR, their average.
    HPol_corr_snr : dict
        The HPol channel-pair correlation SNR, their average.
    """
   
    
    chans = list(wave_bundle.keys())

    AllPol_dict, VPol_dict, HPol_dict = {}, {}, {}
    corr_snr, VPol_corr_snr, HPol_corr_snr = {}, {}, {}

    for chan in chans:

      if(chan in excluded_channels):
        continue

      AllPol_dict.update({chan: wave_bundle[chan]})
      
      if chan in chans[:8]:
        VPol_dict.update({chan: wave_bundle[chan]})
      if chan in chans[8:]:
        HPol_dict.update({chan: wave_bundle[chan]})
      
    corr_snr = calculate_corr_snr(AllPol_dict)  
    VPol_corr_snr = calculate_corr_snr(VPol_dict)
    HPol_corr_snr = calculate_corr_snr(HPol_dict)

    return corr_snr, VPol_corr_snr, HPol_corr_snr
