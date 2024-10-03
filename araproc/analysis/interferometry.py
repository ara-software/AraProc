import array
import copy
import ctypes
import logging
import numpy as np
import os
import ROOT

from araproc.framework import waveform_utilities as wu
from araproc.analysis import dedisperse as dd
from araproc.analysis import cw_filter as cwf

class RayTraceCorrelatorWrapper:

    """
    A wrapper class around the AraRoot RayTraceCorrelator

    This things has the responsibility to manage interferometry
    via the AraRoot RayTraceCorrelator
    
    ...

    Attributes
    ----------
    station_id : int
       The station ID
    correlators : dictionary
        A dictionary.
        Keys = radii
        Values = the ray trace correlator instance

    """

    def __init__(self,
                 station_id : int,
                 num_channels : int =  16 # by default 16
                 ):
        if station_id not in [1, 2, 3, 4, 5]:
            raise KeyError(f"Station {station_id} is not supported")
        if not isinstance(num_channels, int):
            raise ValueError(f"{num_channels} needs to be an int")
        if (num_channels < 0) or not np.isfinite(num_channels):
            raise ValueError(f"{num_channels} is invalid (negative, not finite, etc.)")

        self.station_id = station_id
        self.num_channels = num_channels

        self.correlators = {}


    def __sanitize_inputs(self,
                           path_to_dir_file,
                           path_to_ref_file,
                           ):
        
        # check if they gave us strings
        if not isinstance(path_to_dir_file, str):
            raise TypeError("Path to dir file must be a string")
        if not isinstance(path_to_ref_file, str):
            raise TypeError("Path to ref file must be a string")

        # check if the file exists
        if not os.path.exists(path_to_dir_file):
            raise FileNotFoundError(f"Direct root File ({path_to_dir_file}) not found")
        if not os.path.exists(path_to_ref_file):
            raise FileNotFoundError(f"Ref root File ({path_to_ref_file}) not found")

        # check they gave us files and not directories
        if not os.path.isfile(path_to_dir_file):
            raise ValueError(f"{path_to_dir_file} looks like a directory, not a file")
        # check they gave us files and not directories
        if not os.path.isfile(path_to_ref_file):
            raise ValueError(f"{path_to_ref_file} looks like a directory, not a file")
        
        # now that we're sure these are clean, we can assign the variables
        self.path_to_data_file = path_to_dir_file
        self.path_to_pedestal_file = path_to_ref_file


    def add_rtc(self, 
                radius : float,
                path_to_dir_file : str,
                path_to_ref_file : str,
                ref_name = None, 
                angular_size : float  = 1., 
                ):

        """
        Parameters
        ----------
        ref_name : None
            The refernece name you want used as the key in the self.correlators dictionary.
            If you don't put anything in, the key will be set to the radius.
        radius : float
            The radius for this correlator in meters.
        path_to_dir_file : str
            the full path to the direct ray tracing file
        path_to_pedestal_file : str
            the full path to the reflected file
        angular_side : float
            Size of the angular grid in degrees.
            Default is 1 degrees.

        """
        
        # open the file, establish the tree, and its properties
        self.__sanitize_inputs(path_to_dir_file, path_to_ref_file)

        if ref_name is None:
            ref_name = radius

        if (radius < 0) or not np.isfinite(radius):
            raise ValueError(f"{radius} is invalid")

        if (angular_size < 0) or not np.isfinite(angular_size):
            raise ValueError(f"{angular_size} is invalid")

        the_correlator = ROOT.RayTraceCorrelator(self.station_id,
                                                 self.num_channels,
                                                 radius,
                                                 angular_size,
                                                 path_to_dir_file,
                                                 path_to_ref_file
                                                 )
        ROOT.SetOwnership(the_correlator, True) # take ownership
        the_correlator.LoadTables()
        self.correlators[ref_name] = the_correlator

    def get_pairs(pol : str = "vpol",
                  excluded_channels =  None
                  ):
        print("hello")

class StandardReco:


    def __init__(self,
                 station_id : int,
                 ):
        if station_id not in [1, 2, 3, 4, 5]:
            raise KeyError(f"Station {station_id} is not supported")
        self.station_id = station_id
        
        num_channels_library = {
            1 : 16,
            2 : 16,
            3 : 16,
            4 : 16,
            5 : 16
        }

        excluded_channels_library = {
            1 : np.array([], dtype=int),
            2 : np.array([15], dtype=int),
            3 : np.array([], dtype=int),
            4 : np.array([], dtype=int),
            5 : np.array([], dtype=int),
        }
        excluded_channels = excluded_channels_library[self.station_id]
        excluded_channels_vec = ROOT.std.vector("int")(excluded_channels)


        # each station has a slightly different distance for the cal pulser reco,
        # so look that up
        calpulser_r_library = {
            1 : "40.00",
            2 : "40.00",
            3 : "40.00",
            4 : "40.00",
            5 : "40.00"
        }

        self.num_channels = num_channels_library[station_id]
        self.excluded_channels = excluded_channels_library[station_id]
        self.station_id = station_id

        self.rtc_wrapper = RayTraceCorrelatorWrapper(self.station_id)

        # always add a "nearby" correlator for 41m away
        dir_path = os.path.join("/data/user/brianclark/rt_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_50_radius_{calpulser_r_library[station_id]}_angle_1.00_solution_0.root"
                                )
        ref_path = os.path.join("/data/user/brianclark/rt_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_50_radius_{calpulser_r_library[station_id]}_angle_1.00_solution_1.root"
                                )
        self.rtc_wrapper.add_rtc(ref_name = "nearby",
                radius=float(calpulser_r_library[station_id]),
                path_to_dir_file=dir_path,
                path_to_ref_file=ref_path
                )

        # always add a "distant" correlator for 300m away
        dir_path = os.path.join("/data/user/brianclark/rt_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_50_radius_300.00_angle_1.00_solution_0.root"
                                )
        ref_path = os.path.join("/data/user/brianclark/rt_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_50_radius_300.00_angle_1.00_solution_1.root"
                                )
        self.rtc_wrapper.add_rtc(ref_name = "distant",
                radius=300,
                path_to_dir_file=dir_path,
                path_to_ref_file=ref_path
                )
        
        # we need a geomtool
        geom_tool = ROOT.AraGeomTool.Instance()
        ROOT.SetOwnership(geom_tool, True)
        the_pairs_v = self.rtc_wrapper.correlators["distant"].SetupPairs(self.station_id,
                                                       geom_tool,
                                                       ROOT.AraAntPol.kVertical,
                                                       excluded_channels_vec
                                                       )
        ROOT.SetOwnership(the_pairs_v, True) # take posession
        self.pairs_v = the_pairs_v
        the_pairs_h = self.rtc_wrapper.correlators["distant"].SetupPairs(self.station_id,
                                                       geom_tool,
                                                       ROOT.AraAntPol.kHorizontal,
                                                       excluded_channels_vec
                                                       )
        ROOT.SetOwnership(the_pairs_h, True) # take posession
        self.pairs_h = the_pairs_h


    def go_basic_reco(self, waveform_bundle):

        # set up the waveform map the way AraRoot wants it
        # as a std::map<int, TGraph*>
        waveform_map = ROOT.std.map("int", "TGraph*")()
        ROOT.SetOwnership(waveform_map, True)
        for chan_i in waveform_bundle.keys():
            waveform_map[chan_i] = waveform_bundle[chan_i]
        
        # get the correlation functions
        corr_functions_v = self.rtc_wrapper.correlators["distant"].GetCorrFunctions(
                                                                    self.pairs_v,
                                                                    waveform_map
                                                                    )
        corr_functions_h = self.rtc_wrapper.correlators["distant"].GetCorrFunctions(
                                                                    self.pairs_h,
                                                                    waveform_map
                                                                    )

        dir_map = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
            0
        )
        ROOT.SetOwnership(dir_map, True)


        # It's really important to pass memory ownership of the correlation
        # function map back to python only AFTER we are completely done using them.
        # Otherwise this memory leaks, for reasons I don't totally understand.
        ROOT.SetOwnership(corr_functions_v, True)
        ROOT.SetOwnership(corr_functions_h, True)
        for c in corr_functions_v:
            ROOT.SetOwnership(c, True)
        for c in corr_functions_h:
            ROOT.SetOwnership(c, True)

        del corr_functions_v
        del corr_functions_h
        del waveform_map
        return dir_map
