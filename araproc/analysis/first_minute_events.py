import numpy as np


def is_first_minute_event(first_unix, useful_event):

    """
    Identifies events that occur during the first minute of a run.

    Parameters
    ----------
    first_unix: signed int
        UnixTime of the first event calculated from d.get_useful_event(0).unixTime.
    useful_event : object
        The useful event object containing event information.

    Returns
    -------
    is_first_min_evt: bool
        True: if an event is from first minute, False: if it isn't
    """

    unix_cut = first_unix + 60
    unix_time = useful_event.unixTime

    is_first_min_evt = unix_time < unix_cut

    return is_first_min_evt



def first_minute_events(d):
    
    """
    Collect event number in the first minute of a run.

    Parameters
    ----------
    d: araproc.framework.dataset.AnalysisDataset object 
        The dataset object which is called during run processing. It takes following inputs
            station_id: int,
            path_to_data_file: string,
            path_to_pedestal_file (optional): string,
            is_simulation: bool
    
    Returns
    -------
    first_minute_events: list
        List of event number in the first minute of a run.
    """

    num_evts = d.num_events
    first_unix = d.get_useful_event(0).unixTime

    first_minute_events = []
    for e in range(num_evts):

        useful_event = d.get_useful_event(e)
        
        is_first_min_evt = is_first_minute_event(first_unix, useful_event)

        if is_first_min_evt == True:
           first_minute_events.append(useful_event.eventNumber)
        else:
           break
 
    return first_minute_events



