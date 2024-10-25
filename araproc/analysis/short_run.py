import numpy as np


def is_short_run(d, event_cut = 1000, duration_cut = 3.0):

    """
    Identifies if a run is short.

    Parameters
    ----------
    d: araproc.framework.dataset.AnalysisDataset object 
        The dataset object which is called during run processing. It takes following inputs
            station_id: int,
            path_to_data_file: string,
            path_to_pedestal_file (optional): string,
            is_simulation: bool

    event_cut: int
        Threshold on the number of events for short run (here for 10% data).
    duration_cut: float
        Threshold on the run duration (in minutes) for short run (here for 10% data).

    Returns
    -------
    is_short_r: bool
        True: if a run is short, False: if it isn't.
    """

    num_evts = d.num_events
    first_unix = d.get_useful_event(0).unixTime
    last_unix = d.get_useful_event(num_evts - 1).unixTime

    total_duration = float(np.abs(last_unix - first_unix)/60.0)
    time_flag = total_duration < duration_cut

    total_event = num_evts
    event_flag = total_event < event_cut

    is_short_r = np.logical_or(time_flag, event_flag)

    return is_short_r
