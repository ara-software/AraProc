import ctypes
import logging
import numpy as np
import os
import ROOT

from araproc.framework import waveform_utilities as wu
from araproc.analysis import dedisperse as dd

class AraDataset:

    """
    A class for representing an ARA dataset

    This is really a wrapper class that opens ARA root files, loads calibrators, etc.
    This wrapper currently only supports

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
    num_rf_channels: int
        number of RF channels for this dataset (usually just 16)
    rf_channel_indices : list of ints
        List of the channel indices (used for looping)

    """

    def __init__(self, 
                 path_to_data_file : str,
                 path_to_pedestal_file : str,
                 ):

        """
        Parameters
        ----------
        path_to_data_file : str
            The full path to the ARA event root data file
        path_to_pedestal_file : str
            The full path to the accompanying ARA pedestal file

        """
          
        self.path_to_root_file = None
        self.path_to_pedestal_file = None
        self.root_tfile = None
        self.event_tree = None
        self.run_number = None
        self.station_id = None
        self.num_events = None
        self.calibrator = None

        self.num_rf_channels = 16 # hard coded, but a variable
        self.rf_channel_indices = np.arange(0, self.num_rf_channels).tolist() # make indices
        
        # open the file, establish the tree, and its properties
        self.__sanitize_input_strings(path_to_data_file, path_to_pedestal_file)
        self.__open_tfile_load_ttree()
        self.num_events = self.event_tree.GetEntries()
        self.__establish_run_number() # set the run number
        self.__establish_station_id() # set the station ID
        self.__load_pedestal() # load the pedestals

        # establish the properties of the dedisperser
        self.__phase_spline = dd.load_arasim_phase_response_as_spline()

    def __sanitize_input_strings(self,
                           path_to_data_file,
                           path_to_pedestal_file):
        
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

    def __open_tfile_load_ttree(self):

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
            logging.debug("Successfully got eventTree")
        except:
            logging.critical("Loading the eventTree failed")
            self.root_tfile.Close() # close the file
            raise

        # load up the ARA event
        self.raw_event_ptr = ROOT.RawAtriStationEvent()
        try:
            self.event_tree.SetBranchAddress("event",ROOT.AddressOf(self.raw_event_ptr))
            logging.debug("Successfully assigned RawAtriStationEvent branch")
        except:
            logging.critical("Assigning the rawEventPtr in the eventTree failed")
            self.root_tfile.Close() # close the file
            raise

    def __establish_run_number(self):

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

    def __establish_station_id(self):

        """
        This function works to find the station id.
        It's a bit hacky. The reason we have to be a bit hacky is that
        pyroot doesn't seem to convert the UChar_t that's used in AraRoot
        to store the stationId. So we have to do a song and dance to get it.

        What this does is first uses the TTree Draw function to pull
        the information from the tree for event zero, but using "goff"
        to suppress the graphical output.
        GetV1 stores the array of values you requested for each entry in the ttree.
        Which we force root to internally treat as a float, and *that* seems
        to properly convert the UChar_t. We can then cast it to an int.

        """
        try:
            test = self.event_tree.Draw("abs(event.stationId)", "Entry$==0", "goff")
            self.station_id = int(np.frombuffer(self.event_tree.GetV1(), np.dtype('float'), test)[0])
            logging.debug(f"Got the station id {self.station_id}")
        except:
            logging.critical("Gettig the station id from the eventTree failed")
            raise

    def __load_pedestal(self):
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

    def get_useful_event(self, 
                             event_idx : int = None
                             ):
        
        """
        Fetch a specific calibrated event

        Parameters
        ----------
        event_idx : int
            The ROOT event index to be passed to GetEntry().
            Please note this is the ROOT TTree event index!
            Not the rawAtriEvPtr->eventNumber variable!

        Returns
        -------
        calibrated_event : UsefulAtriStationEvent
            A fully calibrated UsefulatriStationEvent
        """

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

        useful_event = None
        try:
            useful_event = ROOT.UsefulAtriStationEvent(self.raw_event_ptr,
                                                           ROOT.AraCalType.kLatestCalib)
            logging.debug(f"Got calibrated event {event_idx}")
        except:
            logging.critical(f"Calibrating event index {event_idx} failed.")
            raise 
        
        return useful_event
    
    def get_waveforms(self,
                      useful_event = None,
                      interp_tstep : float = 0.5, # ns
                      which_traces = "dedispersed"
                      ):
                     
        """
        Get waveforms for a calibrated event.

        Parameters
        ----------
        useful_event : UsefulAtriStationEvent
            The pointer to the UsefulAtriStationEvent
        
        which_traces : str
            The type of traces you want.
            Currently supports "calibrated", "interpolated", "dedispersed".
            If "interpolated" or "dedispersed" is selected, the returned
            waveform is itnerpolatd with a time bases of timestep "inter_tstep"
        
        interp_tstep : float
            The time sampling you want used for the interpolated traces

        Returns
        -------
        waveforms : 
            A dictionary, mapping RF channel ID to waveforms.
            Keys are channel id (an integer)
            Values are TGraphs
        """

        if which_traces not in ["calibrated", "interpolated", "dedispersed"]:
            raise KeyError(f"Requested waveform treatment ({which_traces}) is not supported")

        # first, get the standard calibrated waves
        cal_waves = {}
        for ch in self.rf_channel_indices:
            try:
                wave = useful_event.getGraphFromRFChan(ch)
                ROOT.SetOwnership(wave, True) # give python full ownership (see README.md for more info)
                cal_waves[ch] = wave
                logging.debug(f"Got channel {ch}")
            except:
                logging.critical(f"Getting the wave for ch {ch} failed")
                raise
        
        if which_traces == "calibrated":
            return cal_waves
    
        # for anything else, we need interpolation

        interp_waves = {}
        for chan_key, wave in cal_waves.items():
            try:
                interp_wave = ROOT.FFTtools.getInterpolatedGraph(wave,interp_tstep)
                ROOT.SetOwnership(interp_wave, True) # give python full ownership
                interp_waves[chan_key] = interp_wave
                logging.debug(f"Got and interpolated channel {chan_key}")
            except:
                logging.critical(f"Interpolating wave for ch {chan_key} failed")
                raise
        
        # if they wanted interpolated waves, just return those
        if which_traces == "interpolated":
            return interp_waves

        # and if they want a dedispersed wave, do that too
        dedispersed_waves = {}
        for chan_key, wave in interp_waves.items():
            try:
                times, volts = wu.tgraph_to_arrays(wave)
                times_dd, volts_dd = dd.dedisperse_wave(times, 
                                                        volts, 
                                                        self.__phase_spline
                                                        )
                
                dedispersed_waves[chan_key] = wu.arrays_to_tgraph(times_dd, volts_dd)
            except:
                logging.critical(f"Dedispersing wave for ch {chan_key} failed")
                raise
        
        if which_traces == "dedispersed":
            return dedispersed_waves
