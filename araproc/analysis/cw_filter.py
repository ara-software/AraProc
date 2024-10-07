import logging
import ROOT
import yaml
import copy

import importlib.resources as pkg_resources
from . import data

def get_filters(station_id):

    """
    Load the yaml file containing the filter settings for this station.

    Parameters
    ----------
    station : int
        The station for which we want to load filter configurations.

    Returns
    -------
    cw_filters : dictionary
        A dictionary of filters.
        The keys are filter names (e.g. "filt1").
        The values are the FFTtools.sineSubtract filter objects, configured
        to do the filtering.
    """

    if station_id not in [1, 2, 3, 4, 5]:
        raise KeyError(f"Station {station_id} is not supported")

    file = pkg_resources.open_text(data, 
                                   "analysis_configs.yaml")
    file_content = yaml.safe_load(file)

    cw_filters = {}

    try:
        this_station_filter_configs = file_content[f"station{station_id}"]
    except:
        logging.error(f"Could not find station {station_id} in the cw filter config file")
        raise
    
    if this_station_filter_configs is not None:
        for filter_name, config_settings in this_station_filter_configs.items():
            the_filter = ROOT.FFTtools.SineSubtract(3,
                                                   config_settings["min_power_ratio"],
                                                   False)
            the_filter.setVerbose(False)
            the_filter.setFreqLimits(config_settings["min_freq"], 
                                    config_settings["max_freq"])
            ROOT.SetOwnership(the_filter, True) # give python full ownership
            cw_filters[filter_name] = the_filter

    file.close()

    return cw_filters

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
        This is normally the output of `cw_filter.get_filters`.
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
        This is normally the output of `cw_filter.get_filters`.
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
