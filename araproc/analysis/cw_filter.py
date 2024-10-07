import logging
import ROOT
import yaml
import copy

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
            local_waveform = filter.subtractCW(latest_waveform, -1)
            ROOT.SetOwnership(local_waveform, True) # take posession 
            del latest_waveform
            latest_waveform = local_waveform # update the latest waveform
    
    return latest_waveform

def apply_filters(cw_filters, waveform_bundle):

    """
    Apply CW filters to a bundle of waveforms.

    Parameters
    ----------
    cw_filters : dictionary
        A dictionary of sine subtract filters to be applied.
        Key is a string, and is the name of the filter from the config file.
        Value is the sine subtract filter to applied.
    waveform_bundle : dictionary
        A dictionary of waveforms to be processed.
        Key is channel id.
        Value is the TGraph to be filtered.

    Returns
    -------
    filtered_waveforms : dictionary
        A dictionary of filtered wavforms.
        The keys channel ids.
        The values are TGraphs after the application of all filters.
    """

    filtered_waveforms = {}
    for ch_id, wave in waveform_bundle.items():
        filtered_waveforms[ch_id] = apply_filters_one_channel(cw_filters, wave)


    return filtered_waveforms
