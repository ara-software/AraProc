import array
import copy
import ctypes
import logging
import numpy as np
import ROOT

from araproc.framework import file_utilities as futil

class SensorHkWrapper:

    """
    A class for representing an ARA sensor housekeeping dataset.

    This wraps around ara sensor housekeeping

    ...

    Attributes
    ----------
    path_to_hk_file : str
       the full path to the sensor hk root file (event data)
    root_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened sensor hk file
    sensor_tree : ROOT TTree
        the sensorHkTree
    num_data_points: int
        number of data points recorded in the sensor hk file
    """

    def __init__(self,
                 path_to_hk_file : str = None,
                 ):
        
        self.path_to_hk_file = None
        self.root_tfile = None
        self.sensor_tree = None
        self.num_entries = None

        if futil.file_is_safe(path_to_hk_file):
            self.path_to_hk_file = path_to_hk_file
        else:
            raise Exception(f"{path_to_hk_file} has a problem!")
        
        self.__open_tfile_load_ttree()
        self.num_data_points = self.sensor_tree.GetEntries()
    
    def __open_tfile_load_ttree(self):

        # open the TFile
        try:
            self.root_tfile = ROOT.TFile(self.path_to_hk_file, "READ")
            logging.debug(f"Successfully opened {self.path_to_hk_file}")
        except:
            logging.critical(f"Opening {self.path_to_hk_file} failed")
            raise

        # set the TTree
        try:
            self.sensor_tree = self.root_tfile.Get("sensorHkTree")
            logging.debug("Successfully got sensorHkTree")
        except:
            logging.critical("Loading the sensorHkTree failed")
            self.root_tfile.Close() # close the file
            raise

        # load up the sensor hk data point
        self.sensor_hk_ptr = ROOT.AtriSensorHkData()
        
        # for whatever reason, pyroot is happier when C++ remains in control
        # when this is True, or not called at all, garbage collection triggers a segfault
        # so let's put this here I guess (this will probably come back to bite me...)
        ROOT.SetOwnership(self.sensor_hk_ptr, False) 
        
        try:
            self.sensor_tree.SetBranchAddress("sensorHk",ROOT.AddressOf(self.sensor_hk_ptr))
            logging.debug("Successfully assigned AtriSensorHkData branch")
        except:
            logging.critical("Assigning the sensor_hk_ptr in the sensorHkTree failed")
            self.root_tfile.Close() # close the file
            raise
    
    def get_data_point(self, point_idx):
        if point_idx is None:
            raise KeyError(f"Requested data point index {point_idx} is invalid")
        if point_idx >= self.num_data_points:
            raise KeyError(f"Requested data point index {point_idx} exceeds number of events in the run ({self.num_data_points})")
        if point_idx <0:
            raise KeyError(f"Requested event index {point_idx} is invalid (negative)")

        try:
            self.sensor_tree.GetEntry(point_idx)
            logging.debug(f"Called root get entry {point_idx}")
        except:
            logging.critical(f"Getting entry {point_idx} failed.")
            raise
        
        this_sensor_hk_ptr = None
        try:
            this_sensor_hk_ptr = copy.deepcopy(self.sensor_hk_ptr)
            logging.debug(f"Copied hk data point {point_idx}")
        except:
            logging.critical(f"Copying data point index {point_idx} failed.")
            raise 
        ROOT.SetOwnership(this_sensor_hk_ptr, False) ## again, c++ needs to retain control of this for... reasons
        return this_sensor_hk_ptr

class EventHkWrapper:

    """
    A class for representing an ARA event housekeeping dataset.

    This wraps around ara event housekeeping

    ...

    Attributes
    ----------
    path_to_hk_file : str
       the full path to the event hk root file (event data)
    root_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened event hk file
    hk_tree : ROOT TTree
        the eventHkTree
    num_data_points: int
        number of data points recorded in the event hk file
    """

    def __init__(self,
                 path_to_hk_file : str = None,
                 ):
        
        self.path_to_hk_file = None
        self.root_tfile = None
        self.hk_tree = None
        self.num_entries = None

        if futil.file_is_safe(path_to_hk_file):
            self.path_to_hk_file = path_to_hk_file
        else:
            raise Exception(f"{path_to_hk_file} has a problem!")
        
        self.__open_tfile_load_ttree()
        self.num_data_points = self.hk_tree.GetEntries()
    
    def __open_tfile_load_ttree(self):

        # open the TFile
        try:
            self.root_tfile = ROOT.TFile(self.path_to_hk_file, "READ")
            logging.debug(f"Successfully opened {self.path_to_hk_file}")
        except:
            logging.critical(f"Opening {self.path_to_hk_file} failed")
            raise

        # set the TTree
        try:
            self.hk_tree = self.root_tfile.Get("eventHkTree")
            logging.debug("Successfully got eventHkTree")
        except:
            logging.critical("Loading the eventHkTree failed")
            self.root_tfile.Close() # close the file
            raise

        # load up the event hk data point
        self.event_hk_ptr = ROOT.AtriEventHkData()
        
        # for whatever reason, pyroot is happier when C++ remains in control
        # when this is True, or not called at all, garbage collection triggers a segfault
        # so let's put this here I guess (this will probably come back to bite me...)
        ROOT.SetOwnership(self.event_hk_ptr, False) 
        
        try:
            self.hk_tree.SetBranchAddress("eventHk",ROOT.AddressOf(self.event_hk_ptr))
            logging.debug("Successfully assigned AtriEventHkData branch")
        except:
            logging.critical("Assigning the event_hk_ptr in the eventHkTree failed")
            self.root_tfile.Close() # close the file
            raise
    
    def get_data_point(self, point_idx):
        if point_idx is None:
            raise KeyError(f"Requested data point index {point_idx} is invalid")
        if point_idx >= self.num_data_points:
            raise KeyError(f"Requested data point index {point_idx} exceeds number of events in the run ({self.num_data_points})")
        if point_idx <0:
            raise KeyError(f"Requested event index {point_idx} is invalid (negative)")

        try:
            self.hk_tree.GetEntry(point_idx)
            logging.debug(f"Called root get entry {point_idx}")
        except:
            logging.critical(f"Getting entry {point_idx} failed.")
            raise
        
        this_event_hk_ptr = None
        try:
            this_event_hk_ptr = copy.deepcopy(self.event_hk_ptr)
            logging.debug(f"Copied hk data point {point_idx}")
        except:
            logging.critical(f"Copying data point index {point_idx} failed.")
            raise 
        ROOT.SetOwnership(this_event_hk_ptr, False) ## again, c++ needs to retain control of this for... reasons
        return this_event_hk_ptr

class ConfigFileWrapper:

    """
    A class for representing an ARA data taking run config file

    This wraps around ara data taking config file

    ...

    Attributes
    ----------
    path_to_config_file : str
       the full path to the config file
    """

    def __init__(self,
                 path_to_config_file : str = None,
                 ):
        
        self.path_to_config_file = None

        if futil.file_is_safe(path_to_config_file):
            self.path_to_config_file = path_to_config_file
        else:
            raise Exception(f"{path_to_config_file} has a problem!")
        
        self.cal_pulser_info = self.parse_cal_pulser_info()
    
    def parse_cal_pulser_info(self, end_key = ';', num_vals = 1):

        with open(self.path_to_config_file, 'r') as c_file:
             config_file_read = c_file.read()

        calpulser_key = ['antennaIceA#I1=', 'antennaIceB#I1=', 'opIceA#I1=', 'opIceB#I1=', 'attIceA#I1=', 'attIceB#I1=']

        cal_pulser_info = {}m

        for key in calpulser_key:

            # check whether there is a same config but commented out one
            # if there is, delete it before begin the search
            old_key = '//'+key
            new_config_file_read = config_file_read.replace(old_key, '')

            # find the index of the key    
            key_idx = new_config_file_read.find(key)

            # check whether key is in the txt_read or not
            if key_idx != -1:
               key_idx += len(key)
               # find the end_key index after key_idx
               end_key_idx = new_config_file_read.find(end_key, key_idx)

               # if there are multiple values for same key
               if num_vals != 1:
                  val = np.asarray(new_config_file_read[key_idx:end_key_idx].split(","), dtype = int)
               else:
                  val = int(new_config_file_read[key_idx:end_key_idx])

            # it there is no key in the txt_read, output numpy nan
            else:
               val = np.full((num_vals), np.nan, dtype = float)

            cal_pulser_info[key.split("#",1)[0]] = val # when assigning the key, drop everything after the #

        return cal_pulser_info
    