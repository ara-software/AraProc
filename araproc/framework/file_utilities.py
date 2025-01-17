import os
import ctypes
import ROOT
import datetime

def file_is_safe(file_path):
    
    is_safe = False

    if not isinstance(file_path, str):
        logging.error(f"Path to file must be a string. Is {type(file_path)}")
    if not os.path.exists(file_path):
        logging.error(f"File ({file_path}) not found")
    if not os.path.isfile(file_path):
        logging.error(f"{file_path} is not a file")
    else:
        is_safe = True
    
    return is_safe

def get_cvmfs_ped_file_name(station_id, run_number):

    if station_id not in const.valid_station_ids:
        raise Exception(f"Station id {station_id} is not supported")

    if not np.isfinite(run_number) or run_number < 0:
        raise Exception(f"Run number {run_number} is not supported")
    
    cvmfs_top_dir = "/cvmfs/icecube.osgstorage.org/icecube/PUBLIC/groups/arasoft/pedestals_v2"
    start = start = run_number - (run_number % 1000)
    stop = start + 999
    file=f"station_{station_id}/{start:07d}-{stop:07d}/station_{station_id}_run_{run_number:07d}.gz"
    return os.path.join(cvmfs_top_dir, file)

def guess_location_of_hk_files(dataset):

    """
    A utility function to guess the location of the housekeeping file
    associated with a run.
    You need to be at the WIPAC datawarehosue for this to work.
    """

    # get the year from the date of the first event in the file
    t_start = dataset.get_raw_event(0).unixTime
    timestamp = datetime.datetime.fromtimestamp(t_start)
    year = timestamp.strftime('%Y')


    # always look for the 100pct file
    file_path = os.path.join("/data/exp/ARA/", str(year), 
                                f"L1/100pct/ARA0{dataset.station_id}")
    
    sensor_hk_file = None
    event_hk_file = None

    for root, dirs, files in os.walk(file_path):
        for file in files:
            if str(dataset.run_number) in file:
                if "sensorHk" in file:
                    sensor_hk_file = os.path.join(root, file)
                if "eventHk" in file:
                    event_hk_file = os.path.join(root, file)
    
    return sensor_hk_file, event_hk_file
