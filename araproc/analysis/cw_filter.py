import logging
import ROOT
import yaml
import copy
import numpy as np

from araproc.framework import constants as const

def apply_filters_one_channel(cw_filters, waveform_in):
    """
    Apply CW filters to a single waveform.
    This function will not delete `waveform` for you,
    and it will return a new waveform.
    So mind your memory.

    Parameters
    ----------
    cw_filters : dictionary
        A dictionary of sine subtract filters to be applied.
        Key is a string, and is the name of the filter from the config file.
        Value is the sine subtract filter to applied.
    waveform_in : ROOT.TGraph
        A single waveform to be filtered.
    
    Returns
    -------
    latest_waveform : ROOT.TGraph
        The filtered waveform as a TGraph.
    """

    latest_waveform = copy.deepcopy(waveform_in) # start here, with a *copy*
    ROOT.SetOwnership(latest_waveform, True) # take posession

    if len(cw_filters) > 0:
        for filter_i, filter in cw_filters.items():
            local_waveform = filter["filter"].subtractCW(latest_waveform, -1)
            ROOT.SetOwnership(local_waveform, True) # take posession 
            del latest_waveform
            latest_waveform = local_waveform # update the latest waveform
    
    return latest_waveform

def apply_filters(cw_filters, waveform_bundle, cw_ids = None):

    """
    Apply CW filters to a bundle of waveforms.

    Parameters
    ----------
    cw_filters : dictionary of dictionaries
        A dictionary of sine subtract filters available to be applied.
        First key is a string, and is the name of the filter from the config file.
        Second key is the filter information, one of: ["filter", "min_freq", "max_freq", "min_power_ratio"]
        Value is the sine subtract filter to applied.
    waveform_bundle : dictionary
        A dictionary of waveforms to be processed.
        Key is channel id.
        Value is the TGraph to be filtered.
    cw_ids : dictionary
        A dictionary of cw id info.
        If None, all set CW filters allowed to filter. Otherwise, only those 
        covering frequencies in the cw id info will be allowed to filter.

    Returns
    -------
    filtered_waveforms : dictionary
        A dictionary of filtered wavforms.
        The keys channel ids.
        The values are TGraphs after the application of all filters.
    """

    # run sanity check to ensure cw id frequencies are reasonable
    check_cw_ids(cw_ids)

    filtered_waveforms = {}
    active_cw_filters_v = get_active_filters(cw_filters, cw_ids, 0)
    active_cw_filters_h = get_active_filters(cw_filters, cw_ids, 8)
    for ch_id, wave in waveform_bundle.items():

        if ch_id in const.vpol_channel_ids:
            active_cw_filters = active_cw_filters_v
        else:
            active_cw_filters = active_cw_filters_h
        filtered_waveforms[ch_id] = apply_filters_one_channel(active_cw_filters, wave)


    return filtered_waveforms

def get_active_filters(cw_filters, cw_ids, chan):
    """
    Apply CW filters to a bundle of waveforms.

    Parameters
    ----------
    cw_filters : dictionary
        A dictionary of sine subtract filters available to be applied..
        Key is a string, and is the name of the filter from the config file.
        Value is the sine subtract filter to applied.
    cw_ids : dictionary
        A dictionary of cw id info.
        If None, all available CW filters are activated. Otherwise, only those 
        covering frequencies in the cw id info will be activated.
    chan : int
        Channel number which filters will be applied to.   
 
    Returns
    -------
    active_filters : dictionary
        A dictionary of sine subtract filters to be applied for this event and channel.
    """

    if cw_ids is None:
        return cw_filters

    active_filters = {}
    
    if chan in const.vpol_channel_ids:
      pol = "v"
    else: 
      pol = "h"

    scan_directions = ["fwd", "bwd"]
    for direction in scan_directions:

        key = f"badFreqs_{direction}_{pol}"
        badFreqs = np.asarray(cw_ids[key])

        
        isFiltered = np.zeros(len(badFreqs)).astype(bool)
        for filter_i, filter in cw_filters.items():

            fmin = filter["min_freq"]
            fmax = filter["max_freq"]
            
            mask = np.logical_and(badFreqs >= fmin, badFreqs <= fmax) # find frequencies covered by this filter  
            isFiltered = np.logical_or(isFiltered, mask) # track frequencies that have covered by a filter
 
            # if filter covers any flagged frequency, activate it 
            if mask.any(): 
                active_filters[filter_i] = filter

        if not isFiltered.all():
            unFilteredFreqs = badFreqs[~isFiltered]
            raise Exception(f"IDed CW at {unFilteredFreqs} GHz has no corresponding filter! Please add one and rerun.")

    # filters are applied in order they appear in dict
    # there's some advantage to apply them in order of descending 
    # min_power_ratio, so let's quickly enforce that
    active_filters = {k : v for k, v in sorted(active_filters.items(), key=lambda x: x[1]["min_power_ratio"], reverse=True)}

    return active_filters

def check_cw_ids(cw_ids):
    """
    Do a quick sanity check of the frequencies in the passed cw_ids.

    Parameters
    ----------
    cw_ids : dictionary
        A dictionary of cw id info.
    """

    if cw_ids is None:
        return

    for key in cw_ids:
        if 'badFreqs' not in key:
            continue 

        badFreqs = np.asarray(cw_ids[key])
        if np.any(badFreqs < 0.):
            raise ValueError("Negative frequency detected in CW IDs. Abort.")
        if np.any(badFreqs > 2.):
            raise ValueError("Frequency >2 detected in CW IDs. Please ensure units are GHz. Abort.")


