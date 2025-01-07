import os
import ctypes
import ROOT

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