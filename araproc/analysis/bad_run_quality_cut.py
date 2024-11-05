from araproc.framework import ara_known_issues
import numpy as np


def get_known_bad_run_events(dataset):

        bad_surface_run = ara_known_issues.get_bad_surface_run(dataset.station_id)
        bad_run = ara_known_issues.get_bad_run(dataset.station_id)
        L0_to_L1_Processing = ara_known_issues.get_L0_to_L1_Processing_run(dataset.station_id)
        ARARunLogDataBase = ara_known_issues.get_ARARunLogDataBase(dataset.station_id)
        software_dominant_run = ara_known_issues.get_software_dominant_run(dataset.station_id)
        ob_bad_run = ara_known_issues.get_obviously_bad_run(dataset.station_id)
        bad_runs = np.concatenate((bad_surface_run, bad_run, L0_to_L1_Processing, ARARunLogDataBase, software_dominant_run, ob_bad_run), axis = None, dtype = int)
        bad_runs = np.unique(bad_runs).astype(int)
        del bad_surface_run, bad_run, L0_to_L1_Processing, ARARunLogDataBase, software_dominant_run

        run_flag = int(dataset.run_number in bad_runs) #ASG should I return "True" instead?

        known_bad_run_events = np.full((dataset.num_events), run_flag, dtype = int)
        del bad_runs, run_flag

        return known_bad_run_events

def get_known_bad_unix_time_events(useful_event, station_id, add_unchecked_unix_time = False):

        bad_unix_event = 0

        if ara_known_issues.get_bad_unixtime(useful_event.unixTime, station_id):
            bad_unix_event = 1

        if add_unchecked_unix_time:
            if ara_known_issue.get_unchecked_unixtime(useful_event.unixTime, station_id):
                bad_unix_event = 1

        return bad_unix_event