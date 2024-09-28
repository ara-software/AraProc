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

        # always load up the calibrated waveforms
        self.__get_calibrated_waveforms(useful_event)
        self.__get_interpolated_waveforms(useful_event)

    def __del__(self):

        """
        Tear down elegantly to prevent memory leaks.
        Most important thing to do here is make sure we wipe all the TGraphs
        that we created while analyzing the event.
        Python's garbage collector is *probably* smart enough to manage this cleanup on its own,
        but better safe than sorry, since who has control of the objects in memory
        is sometimes fuzzy in pyroot.
        """
        for waveset_i, waveset in self.waveform_sets.items():
            for wave_i, wave in waveset.items():
                del wave

    def __get_calibrated_waveforms(self, useful_event):

        """
        Fetch the RF channel TGraphs for a calibrated event.
        This method is intentionally private. The user shouldn't be calling this.
        They should access calibrated waveforms by doing `ana_event.waveform_sets["calibrated"]`.

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
        
        self.waveform_sets["calibrated"] = calibrated_waveforms

    def __get_interpolated_waveforms(self, 
                                   useful_event, 
                                   interp_tstep=0.5):

        """
        Get the calibrated waveforms, and interpolate them.
        This method is intentionally private. The user shouldn't be calling this.
        They should access interpolated waveforms by doing `ana_event.waveform_sets["interpolated"]`.


        Parameters
        ----------
        useful_event : ROOT.UsefulAtriStationEvent()
            The AraRoot UsefulAtriStationEvent.
        interp_tstep: float
            The interpolation timestep (as a float) you want used
        

        Returns
        -------
        None : None
            This private method has the job to set up the interpolated
            waveforms in self.waveform_sets["interpolated"].
        """

        self.waveform_sets["interpolated"] = {}

        if "calibrated" not in self.waveform_sets:
            raise KeyError(f"Calibrated waveforms are not in this event. Something has gone wrong.")
        cal_waves = self.waveform_sets["calibrated"]

        for chan_key, wave in cal_waves.items():
            try:
                print(type(wave))
                self.waveform_sets["interpolated"][chan_key] = ROOT.FFTtools.getInterpolatedGraph(wave,interp_tstep)
                logging.debug(f"Got and interpolated channel {chan_key}")
            except:
                logging.critical(f"Getting or interpolating wave for ch {chan_key} failed")
                raise
