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
    num_entries: int
        number of entries in the sensor hk file
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
        # self.num_events = self.sensor_tree.GetEntries()

    # def __del__(self):
    #     # destructor
    #     print("Call destructor")
    #     # del self.sensor_hk_ptr
    
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

        # load up the ARA event
        self.sensor_hk_ptr = ROOT.AtriSensorHkData()
        try:
            self.sensor_tree.SetBranchAddress("sensorHk",ROOT.AddressOf(self.sensor_hk_ptr))
            logging.debug("Successfully assigned AtriSensorHkData branch")
        except:
            logging.critical("Assigning the sensor_hk_ptr in the sensorHkTree failed")
            self.root_tfile.Close() # close the file
            raise

