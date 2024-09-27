import ctypes
import logging
import os
import ROOT

class AraDataset:

    """
    A class for representing an ARA dataset

    ...

    Attributes
    ----------
    path_to_root_file : str
       the full path to the data root file (event data)
    path_to_pedestal_file : str
        the full path to the pedestal file to be used to calibrate this data file
    root_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened data event file
    event_tree : ROOT TTree
        the sound that the animal makes
    run_number: int
        ARA run number for this dataset
        This will be inferred from the data itself
    station_id: 
        station id from 1->5, 100 (only A1-A5 supported for now)
    num_events: int
        number of data events in the data ROOT file
    calibrator: AraEventCalibrator
        the AraEventCalibrator that carries the appropriate pedestals for this run

    """

    def __init__(self, 
                 station_id : str,
                 path_to_data_file : str,
                 path_to_pedestal_file : str,
                 ):
          
        self.path_to_root_file = None
        self.path_to_pedestal_file = None
        self.root_tfile = None
        self.event_tree = None
        self.run_number = None
        self.station_id = None
        self.num_events = None
        self.calibrator = None

        # sanitize station id
        if not isinstance(station_id, int):
            raise TypeError("station_id must be a string")
        if station_id not in [1, 2, 3, 4, 5, 100]:
            raise ValueError(f"The requested station id ({station_id}) is not supported")
        self.station_id = station_id
        
        # check if they gave us strings
        if not isinstance(path_to_data_file, str):
            raise TypeError("Path to root data file must be a string")
        if not isinstance(path_to_pedestal_file, str):
            raise TypeError("Path to pedestal file must be a string")

        # check if the file exists
        if not os.path.exists(path_to_data_file):
            raise FileNotFoundError(f"Data root File ({path_to_data_file}) not found")
        if not os.path.exists(path_to_pedestal_file):
            raise FileNotFoundError(f"Pedestal File ({path_to_pedestal_file}) not found")

        # check they gave us files and not directories
        if not os.path.isfile(path_to_data_file):
            raise ValueError(f"{path_to_data_file} looks like a directory, not a file")
        # check they gave us files and not directories
        if not os.path.isfile(path_to_pedestal_file):
            raise ValueError(f"{path_to_pedestal_file} looks like a directory, not a file")
        
        # now that we're sure these are clean, we can assign the variables
        self.path_to_data_file = path_to_data_file
        self.path_to_pedestal_file = path_to_pedestal_file
        
        # open the file, establish the tree, and its properties
        self.open_tfile_load_ttree()
        self.num_events = self.event_tree.GetEntries()

        self.establish_run_number() # set the run number
        self.load_pedestal() # load the custom pedestals

    def open_tfile_load_ttree(self):

        # open the TFile
        try:
            self.root_tfile = ROOT.TFile(self.path_to_data_file, "READ")
            logging.debug(f"Successfully opened {self.path_to_data_file}")
        except:
            logging.critical(f"Opening {self.path_to_data_file} failed")
            raise

        # set the TTree
        try:
            self.event_tree = self.root_tfile.Get("eventTree")
            logging.debug(f"Successfully got eventTree")
        except:
            logging.critical("Loading the eventTree failed")
            self.root_tfile.Close() # close the file
            raise

        # load up the ARA event
        self.raw_event_ptr = ROOT.RawAtriStationEvent()
        try:
            self.event_tree.SetBranchAddress("event",ROOT.AddressOf(self.raw_event_ptr))
            logging.debug(f"Successfully assigned RawAtriStationEvent branch")
        except:
            logging.critical("Assigning the rawEventPtr in the eventTree failed")
            self.root_tfile.Close() # close the file
            raise

    def establish_run_number(self):

        run_ptr = ctypes.c_int()
        try:
            self.event_tree.SetBranchAddress("run",run_ptr)
            logging.debug("Setting the run branch address worked")
        except:
            logging.critical("Could not assign run branch address")
            raise

        try:
            self.event_tree.GetEntry(0)
            logging.debug("Got entry zero for purposes of establishing run number")
        except:
            logging.critical("Could not get entry zero for purposes of establishing run number")
            raise
    
        self.run_number = run_ptr.value

    def load_pedestal(self):
        try:
            self.calibrator = ROOT.AraEventCalibrator.Instance()
            logging.debug("Instantiating the AraEventCalibrator was successful")
        except:
            logging.critical("Instantiating the AraEventCalibrator failed")
            raise
        
        try:
            self.calibrator.setAtriPedFile(self.path_to_pedestal_file, 
                                           self.station_id)
        except:
            logging.critical("Setting the AtriPedFile failed")
            raise


    def get_calibrated_event(self, 
                             event_idx : int = None
                             ):

        logging.debug(f"Trying to fetch calibrated event {event_idx}")

        if event_idx is None:
            raise KeyError(f"Requested event index {event_idx} is invalid")
        if event_idx >= self.num_events:
            raise KeyError(f"Requested event index {event_idx} exceeds number of events in the run ({self.num_events})")
        if event_idx <0:
            raise KeyError(f"Requested event index {event_idx} is invalid (negative)")
        
        try:
            self.event_tree.GetEntry(event_idx)
            logging.debug(f"Called root get entry {event_idx}")
        except:
            logging.critical(f"Getting entry {event_idx} failed.")
            raise 

        calibrated_event = None
        try:
            calibrated_event = ROOT.UsefulAtriStationEvent(self.raw_event_ptr,
                                                           ROOT.AraCalType.kLatestCalib)
            logging.debug(f"Got calibrated event {event_idx}")
        except:
            logging.critical(f"Calibrating event index {event_idx} failed.")
            raise 
        
        return calibrated_event
