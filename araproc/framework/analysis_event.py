import ctypes
import logging
import numpy as np
import os
import ROOT

class AraAnalysisEvent:

    """
    A class for AraAnalysis

    This is really a wrapper class that takes a UsefulAtriStationEvent,
    and produces an "analysis ready" object.
    Including things like getting the waveforms interpolated, etc. 

    ...

    Attributes
    ----------
    waveform_sets : dictionary
       A dictionary of waveform sets

    """

    def __init__(self, 
                 station_id : int,
                 run_number : int,
                 useful_event : ROOT.UsefulAtriStationEvent(),
                 ):
          
        self.num_rf_channels = 16 # hard coded, but a variable
        self.rf_channel_indices = np.arange(0, self.num_rf_channels).tolist() # make indices

        self.station_id = station_id
        self.run_number = run_number
        self.event_number = useful_event.eventNumber # make sure the run knows who it is
        
        self.waveform_sets = {}
        self.waveform_sets["calibrated"] = self.get_calibrated_waveforms(useful_event)

    def get_calibrated_waveforms(self, useful_event):

        """
        Fetch the RF channel TGraphs for a calibrated event

        Parameters
        ----------
        event_idx : int
            The ROOT event index to be passed to GetEntry().
            Please note this is the ROOT TTree event index!
            Not the rawAtriEvPtr->eventNumber variable!

        Returns
        -------
        calibrated_waveforms : dictionary
            A dictionary containing calibrated waveforms.
            The keys are integers.
            The values are calibrated TGraphs.
        """

        calibrated_waveforms = {}

        for ch in self.rf_channel_indices:
            try:
                calibrated_waveforms[ch] = useful_event.getGraphFromRFChan(ch)
                logging.debug(f"Got channel {ch}")
            except:
                logging.critical(f"Getting the wave for ch {ch} failed")
                raise
        
        return calibrated_waveforms
