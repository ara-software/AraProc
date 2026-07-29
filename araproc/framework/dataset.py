import array
import copy
import logging
import os
import numpy as np
import ROOT
import yaml
import warnings
from scipy import stats


from araproc.framework import waveform_utilities as wu
from araproc.framework import file_utilities as futil
from araproc.framework import constants as fconst
from araproc.analysis import dedisperse as dd
from araproc.analysis import cw_filter as cwf
from araproc.analysis.snr import get_snr
from araproc.framework.ray_tracer import find_ray_paths

import importlib.resources as pkg_resources
from . import config_files


def get_filters(station_id, analysis_config):

    """
    Load the yaml file containing the filter settings for this station.

    Parameters
    ----------
    station : int
        The station for which we want to load filter configurations.
    config : int
        The analysis config we want to load

    Returns
    -------
    cw_filters : dictionary
        A dictionary of filters.
        The keys are filter names (e.g. "filt1").
        The values are the FFTtools.sineSubtract filter objects, configured
        to do the filtering.
    """

    if station_id not in fconst.valid_station_ids:
        raise KeyError(f"Station {station_id} is not supported")

    file = pkg_resources.open_text(config_files, 
                                   "analysis_configs.yaml")
    file_content = yaml.safe_load(file)

    cw_filters = {}

    try:
        this_station_config = file_content[f"station{station_id}"][f"config{analysis_config}"]
    except:
        logging.error(f"Could not find station {station_id}, config {analysis_config} in the cw config file")
        raise
    
    if this_station_config["filters"] is not None:
        for filter_name, config_settings in this_station_config["filters"].items():
            the_filter = ROOT.FFTtools.SineSubtract(3,
                                                   config_settings["min_power_ratio"],
                                                   False)
            the_filter.setVerbose(False)
            the_filter.setFreqLimits(config_settings["min_freq"], 
                                    config_settings["max_freq"])
            #the_filter.setPeakFindingOption(0) # change to use global max as peak 
            the_filter.setPeakFindingOption(3) # change to only consider peaks above baseline computed by a Savitzky Golay filter 
            ROOT.SetOwnership(the_filter, True) # give python full ownership
            cw_filters[filter_name] = {}
            cw_filters[filter_name]["filter"] = the_filter
            cw_filters[filter_name]["min_freq"] = config_settings["min_freq"]
            cw_filters[filter_name]["max_freq"] = config_settings["max_freq"]
            cw_filters[filter_name]["min_power_ratio"] = config_settings["min_power_ratio"]

    file.close()

    return cw_filters

def find_likely_triggering_interaction(
    report_station,
    detector_station,
    station_id,
):
    """
    Identify the interaction and ray most consistently associated with the
    station trigger.

    Unlike selecting one highest-SNR channel, this function combines the
    SignalBin-versus-Trig_Pass timing evidence from all channels that
    participated in the applicable trigger.

    Parameters
    ----------
    report_station : ROOT AraSim Report::StationReport
        Usually ``report.stations[0]``.
    detector_station : ROOT AraSim Detector::ARA_station
        Usually ``detector.stations[0]``.
    station_id : int
        ARA station ID. Station 6 is treated as the phased-array station.

    Returns
    -------
    interaction_idx, ray_idx : int, int
        Most strongly supported interaction and ray indices.

    Raises
    ------
    RuntimeError
        If no triggered channel contains usable interaction/ray timing
        information.
    """

    triggered_antennas = []

    for string_idx in range(len(report_station.strings)):
        report_string = report_station.strings[string_idx]
        detector_string = detector_station.strings[string_idx]

        for antenna_idx in range(len(report_string.antennas)):
            report_antenna = report_string.antennas[antenna_idx]
            detector_antenna = detector_string.antennas[antenna_idx]

            if int(report_antenna.Trig_Pass) <= 0:
                continue

            antenna_type = int(detector_antenna.type)
            triggered_antennas.append(
                (string_idx, antenna_idx, antenna_type)
            )

    if not triggered_antennas:
        raise RuntimeError(
            "No antennas with a positive Trig_Pass were found."
        )

    if int(station_id) == 6:
        # AraSim's PA coincidence uses VPol antennas on string 0.
        antennas_to_analyze = [
            antenna
            for antenna in triggered_antennas
            if antenna[0] == 0 and antenna[2] == 0
        ]
    else:
        n_triggered_vpol = sum(
            antenna_type == 0
            for _, _, antenna_type in triggered_antennas
        )
        n_triggered_hpol = sum(
            antenna_type == 1
            for _, _, antenna_type in triggered_antennas
        )

        if n_triggered_vpol >= 3 and n_triggered_hpol >= 3:
            antennas_to_analyze = triggered_antennas
        elif n_triggered_vpol >= 3:
            antennas_to_analyze = [
                antenna
                for antenna in triggered_antennas
                if antenna[2] == 0
            ]
        elif n_triggered_hpol >= 3:
            antennas_to_analyze = [
                antenna
                for antenna in triggered_antennas
                if antenna[2] == 1
            ]
        else:
            # Retain useful timing evidence for unusual trigger settings whose
            # multiplicity is lower than the standard three-channel trigger.
            antennas_to_analyze = triggered_antennas

    if not antennas_to_analyze:
        raise RuntimeError(
            "No applicable triggered antennas were found for interaction "
            "selection."
        )

    # SignalBin and SignalExt are transient AraSim members and are not
    # available after reading the output ROOT file. Use the persistent
    # per-channel Likely_Sol result and combine the evidence across all
    # antennas that participated in the applicable station trigger.
    interaction_votes = {}
    interaction_ray_votes = {}

    for string_idx, antenna_idx, _ in antennas_to_analyze:
        antenna = report_station.strings[string_idx].antennas[antenna_idx]
        likely_solution = np.asarray(antenna.Likely_Sol, dtype=int)

        if likely_solution.size < 2:
            continue

        interaction_idx = int(likely_solution[0])
        ray_idx = int(likely_solution[1])

        if interaction_idx < 0 or ray_idx < 0:
            continue

        interaction_votes[interaction_idx] = (
            interaction_votes.get(interaction_idx, 0) + 1
        )

        interaction_ray_votes[(interaction_idx, ray_idx)] = (
            interaction_ray_votes.get((interaction_idx, ray_idx), 0) + 1
        )

    if not interaction_votes:
        raise RuntimeError(
            "Triggered antennas were found, but none had a valid saved "
            "Likely_Sol value."
        )

    # Select the interaction supported by the greatest number of applicable
    # triggering antennas. Use the lower index only as a deterministic tie
    # breaker.
    best_interaction = min(
        interaction_votes,
        key=lambda interaction_idx: (
            -interaction_votes[interaction_idx],
            interaction_idx,
        ),
    )

    rays_for_best_interaction = [
        (ray_idx, vote_count)
        for (interaction_idx, ray_idx), vote_count
        in interaction_ray_votes.items()
        if interaction_idx == best_interaction
    ]

    best_ray = min(
        rays_for_best_interaction,
        key=lambda item: (
            -item[1],
            item[0],
        ),
    )[0]

    return int(best_interaction), int(best_ray)

def find_avg_receipt_ang(
    vertex_position,
    sol_type,
    station_id,
    frequency=0.3,
    accuracy=0.2,
):
    """
    Calculate the average ray-traced receipt angle across a station.

    Each channel is ray traced independently because the antennas have
    different horizontal positions and depths.

    Parameters
    ----------
    vertex_position : array-like
        Vertex ``[x, y, z]`` in meters. The coordinates must be station-centric
        in x and y, with z measured relative to the local surface. Negative z
        therefore indicates a position below the surface.
    sol_type : int
        Ray solution selected by increasing time of flight:

        - 0: direct solution, having the shortest time of flight.
        - 1: reflected/refracted solution, having the second-shortest time of
          flight.
    station_id : int
        ARA station ID. Station 100 uses station 1 geometry, matching the
        existing AraSim convention.
    frequency : float, default 0.3 (to match AraSim)
        Frequency passed to the AraSim ray tracer. 
    accuracy : float, default 0.2 (to match AraSim)
        AraSim path-finding accuracy.

    Returns
    -------
    float
        Arithmetic mean of the available ``TraceRecord.receiptAngle`` values,
        in radians. 
    Raises
    ------
    ValueError
        If ``vertex_position`` is invalid or ``sol_type`` is not 0 or 1.
    RuntimeError
        If the requested ray solution is unavailable for every channel.
    """
    if sol_type not in (0, 1):
        raise ValueError(
            f"sol_type must be 0 (direct) or 1 (reflected); got {sol_type!r}."
        )

    vertex = np.asarray(vertex_position, dtype=float)

    if vertex.shape != (3,):
        raise ValueError(
            "vertex_position must contain exactly three coordinates "
            f"[x, y, z]; got shape {vertex.shape}."
        )

    if not np.all(np.isfinite(vertex)):
        raise ValueError(
            f"vertex_position must contain finite coordinates; got {vertex!r}."
        )

    # AraSim station 100 uses station 1 geometry.
    geometry_station_id = 1 if int(station_id) == 100 else int(station_id)

    geometry_tool = ROOT.AraGeomTool.Instance()
    station_info = geometry_tool.getStationInfo(geometry_station_id)

    receipt_angles = []

    for channel in range(fconst.num_rf_channels):
        antenna_info = station_info.getAntennaInfo(channel)
        antenna_position = [
            float(antenna_info.antLocation[0]),
            float(antenna_info.antLocation[1]),
            float(antenna_info.antLocation[2]),
        ]

        paths = find_ray_paths(
            source_xyz=vertex,
            receiver_xyz=antenna_position,
            frequency=frequency,
            accuracy=accuracy,
        )

        # Paths are sorted by time of flight:
        # index 0 is direct and index 1 is reflected/refracted.
        if len(paths) <= sol_type:
            continue

        receipt_angle = float(paths[sol_type].receiptAngle)

        if np.isfinite(receipt_angle):
            receipt_angles.append(receipt_angle)

    if not receipt_angles:
        solution_name = "direct" if sol_type == 0 else "reflected"
        return fconst.no_raytrace_solution

    return float(np.mean(receipt_angles))

class DataWrapper:

    """
    A class for representing an ARA dataset.

    This wraps around real ara data.

    ...

    Attributes
    ----------
    path_to_data_file : str
       the full path to the data root file (event data)
    path_to_pedestal_file : str
        the full path to the pedestal file to be used to calibrate this data file
    root_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened data event file
    event_tree : ROOT TTree
        the sound that the animal makes
    cw_id_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened cw id file
    cw_id_tree : ROOT TTree
        the sound the the mousetrap makes 
    cw_id_reader : ROOT TTreeReader
        the ROOT TTreeReader of the cw id tree
    always_on_min_cw_id_freq : float
        minimum frequency for CW ID (in GHz), all filters below this frequency are always activated
    always_on_max_cw_id_freq : float
        maximum frequency for CW ID (in GHz), all filters above this frequency are always activated
    run_number: int
        ARA run number for this dataset
        This will be inferred from the data itself
    station_id: 
        station id from 1->5, 100 (only A1-A5 supported for now)
    do_not_calibrate: bool
        a boolean that controls if the calibration (pedestals, voltage ad timing cal files) will be loaded and applied
        If True, then you can only access "raw" events, meaning no waveforms.
        If False, then AraProc will ask AraRoot to perform AraCalType::kLatestCalib (pedestals, and timing/volt cal). User will have access to waveforms.
        Most users will leave this False.
    path_to_cw_ids : str
        the full path to the file containing identified CW frequencies for this data file
    num_events: int
        number of data events in the data ROOT file
    raw_event_ptr : RawAtriStationEvent
        a RawAtriStationEvent ptr from AraRoot
    cw_id_reader_vals : dict of TTreeReaderValues
        dictionary of reader values for different branches
    qual_cuts : AraQualCuts
        An AraRoot AraQualCuts object
    calibrator: AraEventCalibrator
        the AraEventCalibrator that carries the appropriate pedestals for this run
    config : int
        The analysis config for this run, looked up from AraQualCuts
    """

    def __init__(self,
                 path_to_data_file : str = None,
                 path_to_pedestal_file : str = None,
                 station_id : int = None,
                 do_not_calibrate : bool = False,
                 path_to_cw_ids : str = None,
                 always_on_min_cw_id_freq : float = 0.120,
                 always_on_max_cw_id_freq : float = 0.350,
                 ):
        
        self.path_to_data_file = None
        self.path_to_pedestal_file = None
        self.do_not_calibrate = do_not_calibrate
        self.path_to_cw_ids = None
        self.always_on_min_cw_id_freq = always_on_min_cw_id_freq
        self.always_on_max_cw_id_freq = always_on_max_cw_id_freq
        self.root_tfile = None
        self.event_tree = None
        self.cw_id_tfile = None
        self.cw_id_tree = None
        self.cw_id_reader = None
        self.run_number = None
        self.station_id = None
        self.data_station_id = None
        self.num_events = None
        self.calibrator = None
        self.raw_event_ptr = None
        self.cw_id_readers = None
        self.config = None
        self.__num_rf_readout_blocks = None
        self.__num_soft_readout_blocks = None

        if station_id not in fconst.valid_station_ids:
            raise Exception(f"Station id {station_id} is not supported")
        self.station_id = station_id

        if path_to_data_file[:4] == "http":
            self.path_to_data_file = path_to_data_file
        elif futil.file_is_safe(path_to_data_file):
            self.path_to_data_file = path_to_data_file
        else:
            raise Exception(f"{path_to_data_file} has a problem!")
        
        self.__open_tfile_load_ttree()
        self.__establish_run_number() # set the run number
        self.__establish_data_station_id() # get station id from root file
        self.num_events = self.event_tree.GetEntries()

        # force AraRoot to load the proper sqlite table
        # we need to do this right away to make sure that for A3 > 2018
        # we get the right channel mapping
        self.event_tree.GetEntry(0)
        geo_tool = ROOT.AraGeomTool.Instance()
        geo_tool.getStationInfo(self.station_id, self.raw_event_ptr.unixTime)
        
        self.__assign_config()

        # now pedestals
        if not self.do_not_calibrate:
            self.path_to_pedestal_file = path_to_pedestal_file
            self.__load_pedestal(path_to_pedestal_file) # load the pedestals
        else:
            logging.warning(f"do_not_calibrate = {self.do_not_calibrate}. Data will not be calibrated...")
        
        # now load identified CW frequencies
        self.path_to_cw_ids = path_to_cw_ids
        self.__load_cw_ids(path_to_cw_ids)


    def __open_tfile_load_ttree(self):

        # open the TFile
        try:
            self.root_tfile = ROOT.TFile.Open(self.path_to_data_file, "READ")
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
        run_ptr = array.array('i', [0])  # Int_t buffer
        try:
            rc = self.event_tree.SetBranchAddress("run", run_ptr)
            if rc != 0:
                raise RuntimeError(f"SetBranchAddress('run', …) failed with rc={rc}")
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

        if int(run_ptr[0]) == 0:
            raise ValueError("Extracted run number was 0, something probably went wrong!")

        self.run_number = int(run_ptr[0])

    def __establish_data_station_id(self):

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

        This helper function isn't actually used, but it's nice to have it preserved.
        Since it's kind of clever (thanks Marco!)

        """
        try:
            test = self.event_tree.Draw("abs(event.stationId)", "Entry$==0", "goff")
            self.data_station_id = int(np.frombuffer(self.event_tree.GetV1(), np.dtype('float'), test)[0])
            logging.debug(f"Got the station id {self.station_id}")
        except:
            logging.critical("Getting the station id from the eventTree failed")
            raise

    def __assign_config(self):
        self.qual_cuts = ROOT.AraQualCuts.Instance()
        self.config = self.qual_cuts.getLivetimeConfiguration(self.run_number, self.station_id)

    def get_num_rf_readout_blocks(self):
        if self.__num_rf_readout_blocks is not None:
            # if this param is alread initiated, then use it
            return self.__num_rf_readout_blocks
        else:
            # otherwise, initialize it
            self.__establish_readout_params()
            return self.__num_rf_readout_blocks
    
    def get_num_soft_readout_blocks(self):
        if self.__num_soft_readout_blocks is not None:
            # if this param is alread initiated, then use it
            return self.__num_soft_readout_blocks
        else:
            # otherwise, initialize it
            self.__establish_readout_params()
            return self.__num_soft_readout_blocks

    def __establish_readout_params(self):
        
        print("    [LOG] Establishing the Length of RF and Software Blocks....")
        # Enable only needed branches
        self.event_tree.SetBranchStatus("*", 0)
        self.event_tree.SetBranchStatus("event.numReadoutBlocks", 1)
        self.event_tree.SetBranchStatus("event.triggerInfo[4]", 1)

        total_entries = self.event_tree.GetEntries()
        numReadoutBlocks_values = np.empty(total_entries, dtype=int)
        rf_events = np.empty(total_entries, dtype=bool)
        software_triggers = np.empty(total_entries, dtype=bool)

        # Read all data
        for i in range(total_entries):
            self.event_tree.GetEntry(i)
            
            # Get numReadoutBlocks directly
            numReadoutBlocks_values[i] = int(self.event_tree.event.numReadoutBlocks)

            # Get event types
            event_obj = self.event_tree.event
            rf_events[i] = event_obj.isRFTrigger()
            software_triggers[i] = event_obj.isSoftwareTrigger()
                    
        # turn all branches back on
        self.event_tree.SetBranchStatus("*", 1)
        
        # Get RF trigger events
        rf_readout_blocks = numReadoutBlocks_values[rf_events]
        if len(rf_readout_blocks) > 0:
            rf_mode = int(stats.mode(rf_readout_blocks).mode)//fconst.num_dda
        else:
            rf_mode = int(-10)

        # Get software trigger events  z
        sw_readout_blocks = numReadoutBlocks_values[software_triggers]
        if len(sw_readout_blocks) > 0:
            sw_mode = int(stats.mode(sw_readout_blocks).mode)//fconst.num_dda
        else:
            sw_mode = int(-10)
    
        # Return results (handle None cases)
        self.__num_rf_readout_blocks = rf_mode
        self.__num_soft_readout_blocks = sw_mode
        print(f"    [LOG] DONE: Length of RF and Software Blocks is {self.__num_rf_readout_blocks}, {self.__num_soft_readout_blocks}")

    def __load_pedestal(self, path_to_pedestal_file = None):

        # find the right pedestal
        # if the user provided one, then choose that
        # otherwise, try to find it in cvmfs
        # and if that doesn't work, raise an error
        if path_to_pedestal_file is not None:
            if futil.file_is_safe(path_to_pedestal_file):
                logging.info(f"Will try to load custom ped file: {path_to_pedestal_file}")
                self.path_to_pedestal_file = path_to_pedestal_file
            else:
                raise Exception(f"{path_to_pedestal_file} has a problem!")

        else:
            cvmfs_ped = futil.get_cvmfs_ped_file_name(self.station_id, self.run_number)
            if futil.file_is_safe(cvmfs_ped):
                logging.info(f"Will try to load cvmfs ped file: {cvmfs_ped}")
                self.path_to_pedestal_file = cvmfs_ped
            else:
                raise Exception(f"{cvmfs_ped} has a problem!")


        # with the correct pedestal found, we can actually load it
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

    def __load_cw_ids(self, path_to_cw_ids = None):

        # find the right cw id file
        # if the user provided one, then choose that
        # otherwise, try to find it in cvmfs
        # and if that doesn't work, raise an error
        if path_to_cw_ids is not None:
            if futil.file_is_safe(path_to_cw_ids):
                logging.info(f"Will try to load cw id file: {path_to_cw_ids}")
                self.path_to_cw_ids = path_to_cw_ids
            else:
                raise Exception(f"{path_to_cw_ids} has a problem!")
 
        else:
            print("Automatic loading of cw id files is not yet implemented. Please pass a file path.")
            print("Continuing without CW ID info...")
            return

        # now try to load the file
        try:
            self.cw_id_tfile = ROOT.TFile(self.path_to_cw_ids, "READ")
            logging.debug(f"Successfully opened {self.path_to_cw_ids}")
        except:
            logging.critical(f"Opening {self.path_to_cw_ids} failed")
            raise

        # set the TTree
        try:
            self.cw_id_tree = self.cw_id_tfile.Get("NewCWTree")
            logging.debug("Successfully got NewCWTree")
        except:
            logging.critical("Loading the NewCWTree failed")
            self.cw_id_tfile.Close() # close the file
            raise
       
        # Make the reader and store it
        self.cw_id_reader = ROOT.TTreeReader(self.cw_id_tree)
        # Dictionary to hold branch readers
        self.cw_id_reader_vals = {}
        
        # load up the event number info
        self.cw_id_reader_vals['event_num'] = ROOT.TTreeReaderValue["int"](self.cw_id_reader, "event_num")

        # load up the cw id info
        cw_id_info = ["badFreqs", "badSigmas"]
        scan_directions = ["fwd", "bwd"]
        polarizations = ["v", "h"]
        for info in cw_id_info:
            for direction in scan_directions:
                for pol in polarizations:

                    try:
                        key = f"{info}_{direction}_{pol}"
                        self.cw_id_reader_vals[key] = ROOT.TTreeReaderArray["double"](self.cw_id_reader, key)
                        
                        logging.debug(f"Successfully assigned cw id {key} branch")
                    except:
                        logging.critical(f"Assigning the {key} in the newCWTree failed")
                        self.cw_id_tfile.Close() # close the file
                        raise

        logging.debug("Successfully assigned all newCWTree branches")

    def get_raw_event(self, event_idx):
        """
        Fetch a specific *un*calibrated event

        Parameters
        ----------
        event_idx : int
            The ROOT event index to be passed to GetEntry().
            Please note this is the ROOT TTree event index!
            Not the rawAtriEvPtr->eventNumber variable!

        Returns
        -------
        raw_event : RawAtriStationEvent
            A *un*calibrated RawAtriStationEvent
        """
        if event_idx is None:
            raise KeyError(f"Requested event index {event_idx} is invalid")
        if event_idx >= self.num_events:
            raise KeyError(f"Requested event index {event_idx} exceeds number of events in the run ({self.num_events})")
        if event_idx <0:
            raise KeyError(f"Requested event index {event_idx} is invalid (negative)")

        try:
            nbytes = self.event_tree.GetEntry(event_idx)
            logging.debug(f"Called root get entry {event_idx}")
            if nbytes <= 0:
                raise RuntimeError(f"Bad raw entry loading {event_idx}! File may be corrupted...")
        except:
            logging.critical(f"Getting entry {event_idx} failed.")
            raise
        
        raw_event = None
        try:
            raw_event = copy.deepcopy(self.raw_event_ptr)
            logging.debug(f"Copied raw event {event_idx}")
        except:
            logging.critical(f"Copying Raw event index {event_idx} failed.")
            raise 
        ROOT.SetOwnership(raw_event, True)
        
        return raw_event

    def get_useful_event(self, event_idx):
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
        if event_idx is None:
            raise KeyError(f"Requested event index {event_idx} is invalid")
        if event_idx >= self.num_events:
            raise KeyError(f"Requested event index {event_idx} exceeds number of events in the run ({self.num_events})")
        if event_idx <0:
            raise KeyError(f"Requested event index {event_idx} is invalid (negative)")
        
        if self.do_not_calibrate:
            raise Exception("Dataset is not calibrated! You are not allowed to get a useful event!")

        try:
            nbytes = self.event_tree.GetEntry(event_idx)
            logging.debug(f"Called root get entry {event_idx}")
            if nbytes <= 0:
                raise RuntimeError(f"Bad useful entry loading {event_idx}! File may be corrupted...")
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
        ROOT.SetOwnership(useful_event, True)
 
        return useful_event

    def get_cw_ids(self, event_number):
        """
        Get CW ID info for the loaded event.

        Parameters
        ----------
        event_number : int
          Event number to load CW ID for.
        
        Returns
        -------
        cw_ids : tuple of numpy arrays (vpol, hpol)
           Tuple of numpy arrays containing identified bad frequencies in each
           polarization (aggregated over scan direction) 
        """

        cw_ids_v = []
        cw_ids_h = []       

        if self.cw_id_reader is None:
            return None

        # reset reader -- reader MUST move sequentially thru TTree
        self.cw_id_reader.Restart()
        
        # iterate thru reader until we find the right entry
        # we need to IMMEDIATELY readout the reader data, otherwise there's
        # no guarantee the reader will point to the right (or any) entry
        while self.cw_id_reader.Next():
            if self.cw_id_reader_vals["event_num"].Get()[0] == event_number:        
                scan_directions = ["fwd", "bwd"]
                for direction in scan_directions:
                    key = f"badFreqs_{direction}_v"
                    cw_ids_v.extend(self.cw_id_reader_vals[key])            

                    key = f"badFreqs_{direction}_h"
                    cw_ids_h.extend(self.cw_id_reader_vals[key])            

                cw_ids_v = np.unique(np.asarray(cw_ids_v))
                cw_ids_h = np.unique(np.asarray(cw_ids_h))
                cw_ids = (cw_ids_v, cw_ids_h)

                return cw_ids
  
        # if we get to this point, it means the requested event number was 
        # not in the cw id file
        raise Exception(f"Event number {event_number} was not found in the CW ID file!")

        return

    def get_event_index(self, event_number):
        """
        Quickly get the index for a specific event number.

        Parameters
        ----------
        event_number : int
            The event number requested.

        Returns
        -------
        event_idx : int
            The ROOT event index with the requested event number.
        """
        if event_number is None:
            raise KeyError(f"Requested event number {event_number} is invalid")
        if event_number < 0:
            raise KeyError(f"Requested event number {event_number} is invalid (negative)")

        self.event_tree.Draw("Entry$", f"eventNumber=={event_number}", "goff")
        if(self.event_tree.GetSelectedRows() == 1):
            event_idx = int(self.event_tree.GetV1()[0])
        elif(self.event_tree.GetSelectedRows() > 1):
            raise Exception(f"More than one entry in ROOT file has event number {event_number}!")
        else:
            raise KeyError(f"Requested event number {event_number} not found in ROOT file.")

        return event_idx

class SimWrapper:

    """
    A wrapper class for accessing AraSim files.

    ...

    Attributes
    ----------
    path_to_data_file : str
       the full path to the data root file (event data)
    root_tfile : ROOT TFile
        the pointer to the ROOT TFile that corresponds to the opened data event file
    event_tree : ROOT TTree
        the ARA event tree (looks like data)
    sim_tree : ROOT TTree
        the AraSim AraTree2 (exlcusive to simulation)
    run_number: int
        ARA run number for this dataset
        This will be inferred from the data itself
    station_id: 
        station id from 1->5, 100 (only A1-A5 supported for now)
    num_events: int
        number of data events in the data ROOT file
    """

    def __init__(self,
                 path_to_data_file : str = None,
                 station_id : int = None
                 ):
        
        self.path_to_data_file = None
        self.root_tfile = None
        self.event_tree = None
        self.sim_tree = None
        self.sim_settings_tree = None
        if 'merged' in path_to_data_file:
            basename = os.path.basename(path_to_data_file)
            sim_type = basename.split("_")[0]
            if sim_type == "veff": 
                # veff_lgE19.0_A2_config2_merged_0001.root 
                # veff_lgE20.0_A2_config1_neutrino_only_merged_0001.root
                if len(basename.split("_")) == 6: 
                    sim_type, lgE, station, config, merged, run_root = basename.split("_")
                elif len(basename.split("_")) == 8: 
                    sim_type, lgE, station, config, neutrino, only, merged, run_root = basename.split("_")
                else:
                    raise ValueError(f"Unable to parse run information from file name: {basename}")
                energy = float( lgE[3:] )
                run = int(run_root.split(".")[0])
            elif sim_type == 'nu':
                # nu_A4_config2_merged_0002.root
                # nu_PA_config1_merged_0003.root
                sim_type, station, config, merged, run_root = basename.split("_")
                energy = 0
                run = int(run_root.split(".")[0])
            config = config[-1]
            if station == "PA": 
                station_num = 6
            else:
                station_num = int(station[-1])
            self.run_number = station_num*1_0_000_0000 + int(float(config))*1_000_0000 + int( energy*10 )*1_0000 + run
        elif 'noise' in path_to_data_file:
            basename = os.path.basename(path_to_data_file)
            araout_noise, setup, station, config_txt_run_root = basename.split("_")
            config, txt, run, root = config_txt_run_root.split(".")
            config = config[-1]
            run = int(run[3:])
            energy = 0
            if station == "PA": 
                station_num = 6
            else:
                station_num = int(station[-1])
            self.run_number = station_num*1_0_000_0000 + int(float(config))*1_000_0000 + int( energy*10 )*1_0000 + run
        elif '.run' not in path_to_data_file:
            self.run_number = -10 # if there's no file we will crash anyway
        else:
            # arasim files come in the form /path/AraOut.[setupFile].run[runNo].root so first
            # split on '.run' and grab the last element to get [runNo].root then split on
            # '.root' and grab the first element to be left with [runNo]
            self.run_number = int(path_to_data_file.split('.run')[-1].split('.root')[0])
            
        self.station_id = None
        self.num_events = None
        self.config = None
        
        # AraSim specific stuff
        self.useful_event_ptr = None
        self.report_ptr = None
        self.event_ptr = None
        self.settings_ptr = None

        if station_id not in fconst.valid_station_ids:
            raise Exception(f"Station id {station_id} is not supported")
        self.station_id = station_id

        if futil.file_is_safe(path_to_data_file):
            self.path_to_data_file = path_to_data_file
        else:
            raise Exception(f"{path_to_data_file} has a problem!")
        
        self.__open_tfile_load_ttree()
        self.__assign_config()
        self.num_events = self.event_tree.GetEntries()

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
            self.sim_tree = self.root_tfile.Get("AraTree2")
            self.sim_settings_tree = self.root_tfile.Get("AraTree")
            logging.debug("Successfully got eventTree, AraTree, AraTree2")
        except:
            logging.critical("Loading the eventTree, AraTree, or AraTree2 failed")
            self.root_tfile.Close() # close the file
            raise

        # check to make sure eventTree and AraTree2 have matching number of entries
        # if this isn't true, understanding what's in the file is very, very hard
        if self.event_tree.GetEntries() != self.sim_tree.GetEntries():
            raise Exception("eventTree and AraTree2 have a different number"
                            f"of entries ({self.event_tree.GetEntries()} vs {self.sim_tree.GetEntries})."
                            "This is not supported by AraProc.")

        # load up the ARA event
        self.useful_event_ptr = ROOT.UsefulAtriStationEvent()
        self.report_ptr = ROOT.Report()
        self.event_ptr = ROOT.Event()
        self.detector_ptr = ROOT.Detector()
        self.icemodel_ptr = ROOT.IceModel()
        self.settings_ptr = ROOT.Settings()
        try:
            self.event_tree.SetBranchAddress("UsefulAtriStationEvent",ROOT.AddressOf(self.useful_event_ptr))
            self.sim_tree.SetBranchAddress("report", ROOT.AddressOf(self.report_ptr))
            self.sim_tree.SetBranchAddress("event", ROOT.AddressOf(self.event_ptr))
            self.sim_settings_tree.SetBranchAddress("detector", ROOT.AddressOf(self.detector_ptr))
            self.sim_settings_tree.SetBranchAddress("icemodel", ROOT.AddressOf(self.icemodel_ptr))
            self.sim_settings_tree.SetBranchAddress("settings", ROOT.AddressOf(self.settings_ptr))
            self.event_tree.GetEntry(0)
            self.sim_tree.GetEntry(0)
            self.sim_settings_tree.GetEntry(0)
            logging.debug("Successfully assigned UsefulAtriStationEvent, report, event, detector, icemodel, and settings branch")
        except:
            logging.critical("Assigning the useful_event_ptr, report_ptr, detector_ptr, icemodel_ptr, settings_ptr, or event_ptr failed")
            self.root_tfile.Close() # close the file
            raise

    def __assign_config(self):
        # The only way I seem to be able to reliable get at this variable is to pull it via "Scan".
        # For whatever reason, self.sim_settings_tree.GetEntry(0) makes it mad.
        # In order to suppress the print to the screen called by "Scan", 
        # we temporariliy redirect that ouput to /dev/null.
        # Sorry to be so convoluted...
        ROOT.gSystem.RedirectOutput("/dev/null");
        self.sim_settings_tree.Scan("DETECTOR_STATION_LIVETIME_CONFIG", "", "", 1) # only scan first entry to prevent full file from being loaded
        self.config = int(self.settings_ptr.DETECTOR_STATION_LIVETIME_CONFIG)
        ROOT.gROOT.ProcessLine("gSystem->RedirectOutput(0);")

    def __check_event_idx_sanity(self, event_idx):
        if event_idx is None:
            raise KeyError(f"Requested event index {event_idx} is invalid")
        if event_idx >= self.num_events:
            raise KeyError(f"Requested event index {event_idx} exceeds number of events in the run ({self.num_events})")
        if event_idx <0:
            raise KeyError(f"Requested event index {event_idx} is invalid (negative)")
        return True


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
        
        self.__check_event_idx_sanity(event_idx)

        try:
            nbytes = self.event_tree.GetEntry(event_idx)
            logging.debug(f"Called root get entry {event_idx}")
            if nbytes <= 0:
                raise RuntimeError(f"Bad entry loading {event_idx}! File may be corrupted...")
            for chan in range(16):
                gr = self.useful_event_ptr.getGraphFromRFChan(chan)
                ROOT.SetOwnership(gr, True)
                try:
                    if gr is None or gr.GetN() < 10:
                        raise RuntimeError(f"Bad waveform in {event_idx}, channel {chan}! File may be corrupted...")
                finally:
                    del gr
                 
        except:
            logging.critical(f"Getting entry {event_idx} failed.")
            raise 

        useful_event = copy.deepcopy(self.useful_event_ptr) # make a copy and pass that back
        ROOT.SetOwnership(useful_event, True)
        return useful_event
    
    def get_event_index(self, event_number):
        """
        Quickly get the index for a specific event number.

        Parameters
        ----------
        event_number : int
            The event number requested.

        Returns
        -------
        event_idx : int
            The ROOT event index with the requested event number.
        """
        if event_number is None:
            raise KeyError(f"Requested event number {event_number} is invalid")
        if event_number < 0:
            raise KeyError(f"Requested event number {event_number} is invalid (negative)")

        self.event_tree.Draw("Entry$", f"eventNumber=={event_number}", "goff")
        if(self.event_tree.GetSelectedRows() == 1):
            event_idx = int(self.event_tree.GetV1()[0])
        elif(self.event_tree.GetSelectedRows() > 1):
            raise Exception(f"More than one entry in ROOT file has event number {event_number}!")
        else:
            raise KeyError(f"Requested event number {event_number} not found in ROOT file.")

        return event_idx

    def get_sim_information(self, 
                            event_idx : int = None
                            ):
        """
        Fetch simulation / Monte Carlo truth about a simulated event

        Parameters
        ----------
        event_idx : int
            The ROOT event index to be passed to GetEntry().
            Please note this is the ROOT TTree event index!
            Not the rawAtriEvPtr->eventNumber variable!

        Returns
        -------
        sim_info : dictionary
            A python dictionary of useful information about this event.
            For example, the weight, etc.
            This is likely to change over time, and probably it's best to 
            just look at what objects are getting stored below.
        """

        self.__check_event_idx_sanity(event_idx)
    
        try:
            nbytes = self.sim_tree.GetEntry(event_idx)
            logging.debug(f"Called sim_tree get entry {event_idx}")
            if nbytes <= 0:
                raise RuntimeError(f"Bad entry loading {event_idx}! File may be corrupted...") 
        except:
            logging.critical(f"Getting entry {event_idx} in sim_tree failed.")
            raise 

        # find best guess for the interaction index
        likely_interaction, likely_ray = (
            find_likely_triggering_interaction(
                report_station=self.report_ptr.stations[0],
                detector_station=self.detector_ptr.stations[0],
                station_id=self.station_id,
            )
        )

        sim_info = {}
        sim_info["weight"] = self.event_ptr.Nu_Interaction[0].weight
        sim_info["enu"] = self.event_ptr.nu_prim_energy
        sim_info["pid"] = self.event_ptr.nu_prim_pid
        sim_info["evid"] = self.event_ptr.event_ID
        sim_info["likely_int_id"] = likely_interaction
        sim_info["vertex"] = self.get_AraSim_xyz_position(
            self.detector_ptr.stations[0],
            self.event_ptr.Nu_Interaction[likely_interaction].posnu)
        sim_info["sim_weight_vertex"] = self.get_sim_weight_position(
            self.event_ptr.Nu_Interaction[likely_interaction].posnu,
            self.detector_ptr.stations[0],
            self.icemodel_ptr)
        sim_info["direction"] = (
            self.event_ptr.Nu_Interaction[likely_interaction].nnu.Theta(),
            self.event_ptr.Nu_Interaction[likely_interaction].nnu.Phi()
        ) # Although this references the mostly likely triggering interaction, 
        # the direction of all particles in the same event should be the same

        sim_info['ray_solution'] = likely_ray 
        # 0 means the direction solution likely triggered the detector
        # 1 means the refracted/reflected solution likely triggered
        # -1 means a solution was not determined

        sim_info["is_noise"] = (int(self.settings_ptr.TRIG_ANALYSIS_MODE) == 2) # modes 0 & 1 are for signal, 2 is pure noise

        # get the receipt angle for the triggered antennas
        ray_sol = sim_info['ray_solution']
        sim_info["avg_rec_ang"] = self.find_avg_receipt_ang(sim_info["vertex"], ray_sol, self.station_id)

        return sim_info
    
    def get_AraSim_xyz_position(self, origin, position):
        """
        Return the XY displacement of a position from the origin, and the position depth with respect to
        the surface of the ice above the origin. 

        Parameters
        ----------
        origin : AraSim::Position
            Location serving as the XY origin. The Z-coordinate is taken relative to the ice above this point.
            Usually the station.
        position : AraSim::Position
            Location whose coordinates we are interested in. 
            Usually a cascade.

        Returns
        -------
        dx : float
            X displacement from `position` to `origin` in meters
        dy : float
            Y displacement from `position` to `origin` in meters
        depth : float
            Depth of `origin` with respect to the surface of the ice defined by `ROOT::IceModel`
        """

        # Get XY coordinates of provided origin
        x_origin = origin.GetX() 
        y_origin = origin.GetY()

        # Convert XY position coordinates
        x_position = position.GetX() 
        y_position = position.GetY() 

        # Get depth of the position -- careful here since arasim is pretending the Earth is flat...
        origin_surface = self.icemodel_ptr.Surface( origin.Lon(), origin.Lat()) - origin.R() # surface relative to origin
        origin_point = position.GetZ() - origin.GetZ() # position relative to origin
        z_position = origin_point - origin_surface # position relative to surface

        return x_position-x_origin, y_position-y_origin, z_position

    def get_sim_weight_position(self, vertex, station, icemodel):
        """
        Return the XY displacement of a vertex relative to the station-center for use in sim weight calculations.
        This displaces an event vertically if it is above the surface in AraSim to ensure it is given an invalid weight.

        Parameters
        ----------
        vertex : AraSim::Position
            Location whose coordinates we are interested in. 
            Usually a cascade.
        station : AraSim::ARA_station
            AraSim ARA_station object containing station information.            
        icemodel : AraSim::IceModel
            AraSim ice model object containing Earth geometry information.

        Returns
        -------
        x : float
            X displacement from vertex to station-center in meters 
        y : float
            Y displacement from vertex to station-center in meters
        z : float
            Z displacement from vertex to station-center in meters
        """
           
        # find station center in internal coordinates
        avgPos = ROOT.Position(0., 0., 0.)
        avgX = 0.
        avgY = 0.
        avgZ = 0.
        count = 0
        nStrings = int(station.strings.size())
        for i in range(nStrings):
            nAnts = int(station.strings[i].antennas.size())

            for j in range(nAnts):
                avgX += station.strings[i].antennas[j].GetX()
                avgY += station.strings[i].antennas[j].GetY()
                avgZ += station.strings[i].antennas[j].GetZ()
                count += 1
        
        avgX /= count
        avgY /= count
        avgZ /= count

        avgPos.SetXYZ(avgX, avgY, avgZ) # construct this way to ensure lon/lat are set properly

        # Get coordinates of station 
        station_x = avgPos.GetX()
        station_y = avgPos.GetY()
        station_z = avgPos.GetZ() 

        # get surface radius above station
        station_lon = avgPos.Lon()
        station_lat = avgPos.Lat()
        station_surface = icemodel.Surface(station_lon, station_lat)

        # Get coordinates of vertex 
        vertex_x = vertex.GetX()
        vertex_y = vertex.GetY()
        vertex_z = vertex.GetZ()

        # Convert vertex to station-centered coordinates
        x = vertex_x - station_x
        y = vertex_y - station_y
        z = vertex_z - station_z

        if station_surface < vertex.R(): # vertex above surface -- shift z-position to ensure chord length is 0 and sim weight is invalid
            z += 100e3 # shift up by 100 km

        return x, y, z 

class AnalysisDataset:

    """
    A class for representing an ARA dataset.

    This wraps around both data and simulation sets.

    ...

    Attributes
    ----------
    is_simulation : bool
        Is this a simulated dataset or not? 
    path_to_data_file : str
       the full path to the data root file (event data)
    path_to_pedestal_file : str
        the full path to the pedestal file to be used to calibrate this data file
        This argument is optional optional when building a dataset around real data.
        The DataWrapper will try to find the pedesetal files in cvmfs if this argument is None.
        If it's not None, the user-specified pedestals are used.
        This argument is imcompatible with is_simulation = True (since AraSim output is already calibrated.)
        An error will be thrown if those two things clash.
    do_not_calibrate : bool
        A boolean to control if you want real data to be calibrated or not.
        If True, then you can only access "raw" events, meaning no waveforms.
        If False, then AraProc will ask AraRoot to perform AraCalType::kLatestCalib (pedestals, and timing/volt cal). User will have access to waveforms.
        Most users will leave this False.
    path_to_cw_ids : str
        the full path to the file containing identified CW frequencies for this data file
    always_on_min_cw_id_freq : float
        minimum frequency for CW ID (in GHz), all filters below this frequency are always activated
    always_on_max_cw_id_freq : float
        maximum frequency for CW ID (in GHz), all filters above this frequency are always activated
    run_number: int
        ARA run number for this dataset
        This will be inferred from the data itself
    station_id:
        station id from 1->5, 100 (only A1-A5 supported for now)
    num_events: int
        number of data events in the data ROOT file
    num_rf_channels: int
        number of RF channels for this dataset (usually just 16)
    rf_channel_indices : list of ints
        List of the channel indices (used for looping)        
    interp_tstep : float
        The time step you want wavforms interpolated with for later.
        In nanoseconds.
    dataset_wrapper
        Either a DataWrapper or a SimWrapper based on is_simulation.
        This allows this meta class to access either a SimWrapper or the DataWrapper
        behind the scenes, and the user doesn't have to worry about which is 
        which, except for getting the "is_simulation" flag right.
    excluded_channels : array of ints
        A list of the RF channels you want excluded from the analysis
    exluded_channels_and_vpol: array of ints
        The excluded_channel list, plus channels 0->7 (inclusive).
        This list can be used to calculate hpol only variables, so you get e.g.
        bad channels and the exclusion of the vpol channels.
    exluded_channels_and_hpol: array of ints
        The excluded_channel list, plus channels 8->15 (inclusive).
        This list can be used to calculate pol only variables, so you get e.g.
        bad channels and the exclusion of the hpol channels.
    which_phase_file: str
        by default it loads simulated phase file
        choose which_phase_file = 'lab_measured' for the lab measured phase data 
    """

    def __init__(self, 
                 path_to_data_file : str,
                 station_id : int,
                 path_to_pedestal_file : str = None,
                 interp_tstep : float = 0.5 ,
                 is_simulation : bool = False,
                 do_not_calibrate : bool = False,
                 path_to_cw_ids : str = None,
                 always_on_min_cw_id_freq : float = 0.120,
                 always_on_max_cw_id_freq : float = 0.350,
                 which_phase_file = "ARA_TwoFilter",
                 ):
    
        self.is_simulation = is_simulation
        self.do_not_calibrate = do_not_calibrate
        self.path_to_cw_ids = path_to_cw_ids
        self.always_on_min_cw_id_freq = always_on_min_cw_id_freq
        self.always_on_max_cw_id_freq = always_on_max_cw_id_freq
        self.which_phase = which_phase_file
        if self.is_simulation and self.do_not_calibrate:
            raise Exception(f"Simulation (is_simulation = {self.is_simulation}) and uncalibrated data (do_not_calibrate = {self.do_not_calibrate}) are incompatible settings")
        if self.is_simulation and self.path_to_cw_ids:
            raise Exception(f"Simulation (is_simulation = {self.is_simulation}) and IDed CW filtering (path_to_cw_ids = {self.path_to_cw_ids}) are incompatible settings")

        self.path_to_data_file = None
        self.pedestal_file = None # not required if simulation
        self.run_number = None
        self.station_id = None
        self.data_station_id = None
        self.num_events = None
        self.num_rf_channels = None
        self.rf_channel_indices = None
        self.interp_tstep = None
        self.dataset_wrapper = None
        self.excluded_channels = None
        self.excluded_channels_and_hpol = None
        self.excluded_channels_and_vpol = None

        self.num_rf_channels = 16 # hard coded, but a variable
        self.rf_channel_indices = fconst.rf_channels_ids

        self.rf_channel_polarizations = [ ch//8 for ch in self.rf_channel_indices ]

        if (interp_tstep < 0) or not np.isfinite(interp_tstep):
            raise ValueError(f"Something is wrong with the requested interpolation time step: {interp_tstep}")
        if (interp_tstep < 0.3):
            raise ValueError("Requested interpolation time step is <0.3 ns. SineSubtract has been found to give\n"
                             "\t    incorrect/anomalous behavior for such small timesteps, and CW peaks are not properly\n"
                             "\t    removed. Comment this warning if you want, but PROCEED WITH CAUTION.")
        self.interp_tstep = interp_tstep
        
        if not self.is_simulation:
            self.dataset_wrapper = DataWrapper(path_to_data_file,
                                                 path_to_pedestal_file,
                                                 station_id=station_id,
                                                 do_not_calibrate = self.do_not_calibrate,
                                                 path_to_cw_ids = self.path_to_cw_ids,
                                                 always_on_min_cw_id_freq = self.always_on_min_cw_id_freq,
                                                 always_on_max_cw_id_freq = self.always_on_max_cw_id_freq,
                                             )
        else:
            self.dataset_wrapper = SimWrapper(path_to_data_file,
                                                station_id=station_id
                                                )
        
        if self.dataset_wrapper is None:
            raise Exception("Something went wrong with creating either the sim or data wrapper")
        
        self.run_number = self.dataset_wrapper.run_number
        self.station_id = self.dataset_wrapper.station_id
        if(not self.is_simulation):
          self.data_station_id = self.dataset_wrapper.data_station_id
        self.num_events = self.dataset_wrapper.num_events
        self.config = self.dataset_wrapper.config

        self.__get_excluded_channels() # after we set the station ID and config

        # establish the properties of the dedisperser
        self.__phase_spline = dd.load_arasim_phase_response_as_spline(self.which_phase)

        # and the properties of the CW filter, and the bandpass filters
        self.__cw_filters = get_filters(self.station_id, self.config)
        self.__setup_bandpass_filters()

    def __get_excluded_channels(self):
            
        file = pkg_resources.open_text(config_files, 
                                       "analysis_configs.yaml")
        file_content = yaml.safe_load(file)
        
        try:
            this_station_config = file_content[f"station{self.station_id}"][f"config{self.config}"]
        except:
            logging.error(f"Could not find station {self.station_id}, config {self.config} in the cw config file")
            raise

        self.excluded_channels = np.asarray(this_station_config["excluded_channels"])
        self.excluded_channels_and_vpol = np.concatenate((self.excluded_channels, np.asarray(fconst.vpol_channel_ids)))
        self.excluded_channels_and_hpol = np.concatenate((self.excluded_channels, np.asarray(fconst.hpol_channel_ids)))
        file.close()

    def __setup_bandpass_filters(self):
            
        # apparenty we need this line CUZ REASONS?
        # it must activate something in FFTW behind the scenes. This was very not obvious, 
        # and I don't understand how to fix it tbh...
        dummy = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.LOWPASS, 2, 100)
        ROOT.SetOwnership(dummy, True) # give python full ownership

        nyquist = 1./(2.*self.interp_tstep*1E-9) # interp speed in seconds
        freq_lopass = 875E6/nyquist # lowpass fitler at 875 MHz, in units of nyquist
        freq_hipass = 50E6/nyquist # highpass filter at 50 MHz, in units of nyquist

        lp = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.LOWPASS, 5, freq_lopass)
        hp = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.HIGHPASS, 5, freq_hipass)
        ROOT.SetOwnership(lp, True) # give python full ownership
        ROOT.SetOwnership(hp, True) # give python full ownership

        self.__lowpass_filter = lp
        self.__highpass_filter = hp

    def get_raw_event(self, 
                             event_idx : int = None
                             ):

        if self.is_simulation:
            raise Exception("You are working with simulation! Cannot get a raw event!")

        raw_event = self.dataset_wrapper.get_raw_event(event_idx)
        
        return raw_event  
    
    def get_useful_event(self, 
                             event_idx : int = None
                             ):
        
        if self.do_not_calibrate:
            raise Exception("Dataset is not calibrated! You are not allowed to get a useful event!")

        useful_event = self.dataset_wrapper.get_useful_event(event_idx)
        
        return useful_event
   
    def get_cw_ids(self, event_number):

        # CW ids currently not supported for simulation
        if self.is_simulation:
            return None 

        return self.dataset_wrapper.get_cw_ids(event_number)

    def get_event_cw(self, useful_event): 

        if self.is_simulation:
            return None

        return self.dataset_wrapper.get_event_cw(useful_event)
       
    def get_event_index(self, 
                            event_number : int = None):

        event_idx = self.dataset_wrapper.get_event_index(event_number)

        return event_idx

    def get_event_sim_info(self, 
                             event_idx : int = None
                             ):
        if not self.is_simulation:
            raise ValueError("You requested simulation information, but this dataset is marked as data. Abort!")

        sim_info = self.dataset_wrapper.get_sim_information(event_idx)
        return sim_info

    def get_num_rf_readout_blocks(self):
        if self.is_simulation:
            raise ValueError("You requested data information, but this dataset is marked as simulation. Abort!")
        return self.dataset_wrapper.get_num_rf_readout_blocks()


    def get_num_soft_readout_blocks(self):
        if self.is_simulation:
            raise ValueError("You requested data information, but this dataset is marked as simulation. Abort!")
        return self.dataset_wrapper.get_num_soft_readout_blocks()

    def get_wavepacket(self,
                      useful_event,
                      which_traces = "filtered",
                      crop_times = None,
                      ):
                     
        """
        Get waveforms for a calibrated event.

        Parameters
        ----------
        useful_event : UsefulAtriStationEvent
            The pointer to the UsefulAtriStationEvent
        
        which_traces : str
            The type of traces you want.
            Currently supports "calibrated", "interpolated", "dedispersed", "cw_filtered",
            "bandpassed", "filtered", and "cropped".
            
                "calibrated" : 
                    These are the waveforms direct from AraRoot,
                    e.g. the product of evt.getGraphFromRFChan.
                    They have no processing applied at all.
                "interpolated" :
                    These are calibrated waves with interpolation applied.
                    The interpoaltion timestep is given by "self.interp_tstep"
                    (the time step is set when the dataset object was created)
                "cw_filtered" :
                    These are NOT dedispersed waves and are only 
                    filtered of CW via the FFTtools SineSubtract filter.
                "dedispersed" : 
                    These are interpolated & CW filtered waves with dedispersion applied.
                    For the dedispersion, we assume the phase response
                    of the ARA system as found in AraSim, specifically
                    the "ARA_Electronics_TotalGain_TwoFilters.txt" file.
                    The same response is assumed for all channels.
                "bandpassed" :
                    These are dedispersed waves that are additionally
                    bandpass filtered to remove out of band noise, but are NOT cw filtered. 
                "filtered" : 
                    These are dedispersed waves that are additionally
                    filtered of CW via the FFTtools SineSubtract filter.
                    They are then bandpass filtered to remove out of band noise.
                    Most people should use the "filtered" traces,
                    and so they are the default argument.
                "cropped" :
                    These are filtered waveforms that are cropped to a time range
                    passed by the user through the crop_times argument. Outside this
                    range the waveform is set to 0. 
            As a cheat sheet, here are the steps each comprise:
                "calibrated" : cal
                "interpolated" : cal + int
                "cw_filtered" : cal + int + cwf
                "dedispersed" : cal + int + cwf + dd
                "bandpassed" : cal + int + dd + bp
                "filtered" : cal + int + cwf + dd + bp 
                "cropped" : cal + int + cwf + dd + bp + CROP

        crop_times : dict of float pairs (keys : int, channel number; values : tuple, time range)
            The time range (tmin, tmax) which traces should be cropped to for each channel. 
            To avoid cropping you can set tmin/tmax to a very negative/positive number. If
            a channel is missing from the dictionary, it will be skipped in the cropping process.          
        Returns
        -------
        wavepacket : dict
            A dict with three entries:
              "event" : int  
                Event number
              "waveforms" : dict
                Dict mapping RF channel ID to waveforms.
                Keys are channel id (an integer)
                Values are TGraphs
              "trace_type" : string
                Waveform type requested by which_trace
        """       

        if useful_event is None:
            raise KeyError("Passed useful event is None for some reason")

        if which_traces not in ["calibrated", "interpolated", "dedispersed", "cw_filtered", "bandpassed", "filtered", "cropped"]:
            raise KeyError(f"Requested waveform treatment ({which_traces}) is not supported")

        wavepacket = {}
        wavepacket["event"] = useful_event.eventNumber
        wavepacket["trace_type"] = which_traces

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
            wavepacket["waveforms"] = cal_waves
            return wavepacket
    
        # for anything else, we need interpolation

        interp_waves = {}
        for chan_key, wave in cal_waves.items():
            try:
                interp_wave = ROOT.FFTtools.getInterpolatedGraph(wave,self.interp_tstep)
                ROOT.SetOwnership(interp_wave, True) # give python full ownership
                interp_waves[chan_key] = interp_wave
                logging.debug(f"Got and interpolated channel {chan_key}")
            except:
                logging.critical(f"Interpolating wave for ch {chan_key} failed")
                raise
        
        # if they wanted interpolated waves, just return those
        if which_traces == "interpolated":
            wavepacket["waveforms"] = interp_waves
            return wavepacket

        # and if they want filtered waves
        if which_traces != "bandpassed":
            event_number = useful_event.eventNumber
            cw_ids = self.get_cw_ids(event_number)
            filtered_waves = cwf.apply_filters(self.__cw_filters, interp_waves, cw_ids, self.always_on_min_cw_id_freq, self.always_on_max_cw_id_freq)
        else:
            filtered_waves = interp_waves

        if which_traces == "cw_filtered":
            wavepacket["waveforms"] = filtered_waves
            return wavepacket

        # and if they want a dedispersed wave, do that too
        dedispersed_waves = {}
        for chan_key, wave in filtered_waves.items():
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
            wavepacket["waveforms"] = dedispersed_waves
            return wavepacket


        # and finally, apply some bandpass cleanup filters
        bandpassed_waves = {}
        for chan_key in list(dedispersed_waves.keys()):

            """
            This takes some explaining, why I'm not just calling
                digitalFilter.filterGraph(...)
            
            Basically, the DigitalFilter class function `filterGraph` calls
            the function `filter()`. The `filter` function calls the function
            `filterOut`, which returns a POINTER to the filtered array.
            Pyroot doesn't know how to clean this up properly.
            So here, I create a pointer to an array via `array.array`,
            and pass that pointer myself via the `filterOut` function.
            That way python has proper ownership of it, and can destroy it
            and avoid a memory leak.

            I discovered this leak because when I called `filterGraph`,
            I ended up with a huge memory leak.
            Beware these python <-> c++ interfaces, especially around pointers...
            """

            wave = dedispersed_waves[chan_key]
            n = wave.GetN()
            x = wave.GetX()
            y = wave.GetY()
            y_filt_lp = array.array("d", [0]*len(y))
            y_filt_hp = array.array("d", [0]*len(y))
            self.__lowpass_filter.filterOut(n, y, y_filt_lp)
            self.__highpass_filter.filterOut(n, y_filt_lp, y_filt_hp)

            filt_graph = ROOT.TGraph(n, x, y_filt_hp)
            ROOT.SetOwnership(filt_graph, True)

            del wave
            bandpassed_waves[chan_key] = filt_graph

        if which_traces in ["bandpassed", "filtered"]:
            wavepacket["waveforms"] = bandpassed_waves
            return wavepacket

        # crop waveforms
        if(crop_times is None):
            raise ValueError("crop_times argument must be passed for cropped waveforms!")

        cropped_waves = {}
        for chan_key in list(filtered_waves.keys()):
            
            x, y = wu.tgraph_to_arrays(filtered_waves[chan_key])

            # get the time range to crop to (zero pad outside this range)
            if chan_key in crop_times:
                tmin, tmax = crop_times[chan_key]
            
                mask = np.logical_or(x < tmin, x > tmax)
                y[mask] = 0.             

            cropped_wave = wu.arrays_to_tgraph(x, y)
            ROOT.SetOwnership(cropped_wave, True)

            cropped_waves[chan_key] = cropped_wave

        if which_traces == "cropped":
            wavepacket["waveforms"] = cropped_waves
            return wavepacket


