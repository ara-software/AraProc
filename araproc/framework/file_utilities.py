import os
import ctypes
import ROOT
import datetime
import glob
import numpy as np
from araproc.framework import constants as const
import logging

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

def guess_location_of_hk_files(station_id, run_number):

    """
    A utility function to guess the location of the housekeeping file
    associated with a run.
    You need to be at the WIPAC datawarehosue for this to work.
    """

    inpath = f"/data/exp/ARA/*/L1/100pct/ARA0{station_id}/*/run{run_number:06d}/"

    # first the sensor hk file
    sensor_hk_file = f"sensorHk{run_number:06d}.root"
    sensor_hk_full_path = inpath + sensor_hk_file
    files = glob.glob(sensor_hk_full_path)
    if not files:
        raise Exception(f"Requested run file not found: {fullpath}. Abort.")
    if len(files) > 1:
        raise Exception(f"Requested run found in more than one directory!\n{files}")
    
    sensor_hk_file = files[0]

    files = None # reset

    # now the event hk file
    event_hk_file = f"eventHk{run_number:06d}.root"
    event_hk_full_path = inpath + event_hk_file
    files = glob.glob(event_hk_full_path)
    if not files:
        raise Exception(f"Requested run file not found: {fullpath}. Abort.")
    if len(files) > 1:
        raise Exception(f"Requested run found in more than one directory!\n{files}")
    
    event_hk_file = files[0]

    return sensor_hk_file, event_hk_file

def guess_location_of_daq_config_file(station_id, run_number):


    """
    A utility function to guess the location of the daq config file.
    You need to be at the WIPAC datawarehosue for this to work.
    """

    inpath = f"/data/exp/ARA/*/L1/100pct/ARA0{station_id}/*/run{run_number:06d}/"

    config_file = f"configFile.run{run_number:06d}.dat"
    config_file_full_path = inpath + config_file
    files = glob.glob(config_file_full_path)
    if not files:
        raise Exception(f"Requested run file not found: {config_file_full_path}. Abort.")
    if len(files) > 1:
        raise Exception(f"Requested run found in more than one directory!\n{files}")
    
    config_file = files[0]

    return config_file
