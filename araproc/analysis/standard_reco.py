import ctypes
import numpy as np
import math
import os
import ROOT
from scipy.interpolate import Akima1DInterpolator
import itertools

from araproc.analysis import interferometry as interf
from araproc.framework import constants as const
from araproc.framework import map_utilities as mu
from araproc.framework import waveform_utilities as wfu
from araproc.analysis import snr

class StandardReco:

    """
    A "standard reconstruction" class.

    This class will handle the "standard" set of reconstructions we want to run on every event.
    For this reason, it needs to know some details about each station.
    So, each station has an entry in a dictionary which specifies 
        - the number of channels (this should always be 16),
        - the distance of the cal pulser for that station

    This function sets up the properties (num channels, channel exclusions, etc.).
    Based on this information, it also calcultes the appropriate waveform pairs.
    Then it instantiates the two ray trace correlators wwe are expecting to use.
    One at the pulser distance, with the key/name "nearby."
    And one at 300m, with the key/name "distant."
    They are stored in RTC wrapper which belongs to this class.
    See the list of attributes below.

    Attributes
    ----------
    num_channels : int
        The number of channels in the station to be reeconstructed (should always be 16)
    excluded_channels : array(int)
        The list of channels you want *exclued* from the interferometry
    station_id : int
        The id of the station you want to reco (only station 1-5 supported)
    rtc_wrapper : RayTraceCorrelatorWrapper
        A instance of the RayTraceCorrelatorWrapper described earlier in this module.
    the_pairs_v : std::map<int, std::vector<int>>
        The list of vpol pairs we want in the reconstruction
        This is returned directly by AraRoot.
        https://github.com/ara-software/AraRoot/blob/master/AraCorrelator/RayTraceCorrelator.h#L190C42-L190C52
    the_pairs_h : std::map<int, std::vector<int>>
        The list of hpol pairs we want in the reconstruction
        This is returned directly by AraRoot.
        https://github.com/ara-software/AraRoot/blob/master/AraCorrelator/RayTraceCorrelator.h#L190C42-L190C52
    __arrival_delays_v_distant : std::pair<bins, delays>
        This is a complicated struct, but it stores the delays between the vpol pairs.
        This is returned by AraRoot.
        See https://github.com/ara-software/AraRoot/blob/master/AraCorrelator/RayTraceCorrelator.h#L79
        This for the vpol distant hypothesis.
    __arrival_delays_h_distant : std::pair<bins, delays>
        This is a complicated struct, but it stores the delays between the hpol pairs.
        This is returned by AraRoot.
        See https://github.com/ara-software/AraRoot/blob/master/AraCorrelator/RayTraceCorrelator.h#L79
        This for the hpol distant hypothesis.
    __arrival_delays_v_nearby : std::pair<bins, delays>
        Ditto to above, just nearby distant hypothesis.
    __arrival_delays_h_nearby : std::pair<bins, delays>
        Ditto to above, just nearby distant hypothesis.
    """


    def __init__(self,
                 station_id : int,
                 excluded_channels = np.array([]) # by default no excluded channels
                 ):
        
        if station_id not in const.valid_station_ids:
            raise KeyError(f"Station {station_id} is not supported")
        self.station_id = station_id
    
        if not isinstance(excluded_channels, np.ndarray):
            raise KeyError(f"Excluded channel list needs to be a 1D numpy array")
        if excluded_channels.ndim != 1:
            raise AttributeError(f"Excluded channels has the wrong number of dimensions -- should be 1D only")

        num_channels_library = {
            100 : 16,
            2   : 16,
            3   : 16,
            4   : 16,
            5   : 16
        }

        excluded_channels_vec = ROOT.std.vector("int")()
        for i in excluded_channels:
            excluded_channels_vec.push_back(int(i))

        # each station has a slightly different distance for the cal pulser reco,
        # so look that up

        self.calpulser_r_library = {
            100 : "48.02",
            2 : "42.86",
            3 : "41.08",
            4 : "52.60",
            5 : "49.23"
        }

        # for distant events, we use the same radius
        self.distant_events_r_library = {
            1 : "300",
        }

        self.num_channels = num_channels_library[station_id]
        self.station_id = station_id

        self.rtc_wrapper = interf.RayTraceCorrelatorWrapper(self.station_id)

        # always add a "nearby" correlator for 41m away
        #dir_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
        #                        "PUBLIC/groups/arasoft/raytrace_timing_tables",
        dir_path = os.path.join("/data/ana/ARA/processing/support/raytrace_timing_tables/Updated_Tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_{self.calpulser_r_library[station_id]}_angle_1.00_solution_0.root"
                                )
        #ref_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
        #                        "PUBLIC/groups/arasoft/raytrace_timing_tables",
        ref_path = os.path.join("/data/ana/ARA/processing/support/raytrace_timing_tables/Updated_Tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_{self.calpulser_r_library[station_id]}_angle_1.00_solution_1.root"
                                )
        self.rtc_wrapper.add_rtc(ref_name = "nearby",
                radius=float(self.calpulser_r_library[station_id]),
                path_to_dir_file=dir_path,
                path_to_ref_file=ref_path
                )

        # always add a "distant" correlator for 300m away
        #dir_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
        #                        "PUBLIC/groups/arasoft/raytrace_timing_tables",
        dir_path = os.path.join("/data/ana/ARA/processing/support/raytrace_timing_tables/Updated_Tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_300.00_angle_1.00_solution_0.root"
                                )
        #ref_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
        #                        "PUBLIC/groups/arasoft/raytrace_timing_tables",
        ref_path = os.path.join("/data/ana/ARA/processing/support/raytrace_timing_tables/Updated_Tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_300.00_angle_1.00_solution_1.root"
                                )
        self.rtc_wrapper.add_rtc(ref_name = "distant",
                radius=float(self.distant_events_r_library[1]),
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

        # get the arrival delays
        self.__arrival_delays_v_distant = self.rtc_wrapper.correlators["distant"].GetArrivalDelays(self.pairs_v)
        self.__arrival_delays_h_distant  = self.rtc_wrapper.correlators["distant"].GetArrivalDelays(self.pairs_h)
        self.__arrival_delays_v_nearby = self.rtc_wrapper.correlators["nearby"].GetArrivalDelays(self.pairs_v)
        self.__arrival_delays_h_nearby  = self.rtc_wrapper.correlators["nearby"].GetArrivalDelays(self.pairs_h)


        self.__latest_event_num = -1

    def __calculate_cross_correlations(self, waveform_bundle, pairs, applyHilbert):

        """
        Wrapper for cross-correlation calculation.

        Parameters
        ----------
        waveform_bundle : dict
            A dictionary of the 16 waveforms.
            The key is the RF channel number.
            The value is a TGraph.
            There should be 16 entries, even if you don't intend to use all
            16 traces in your interferometry.
            The exclusions are handled further down under the excluded channels section.       
        pairs : std::map<int, std::vector<int> >
            Channel pairs to form cross-correlations from.
        applyHilbert : bool
            Boolean whether the cross-correlation function is Hilbert enveloped.
 
        Returns
        -------
        cross_correlations : std::vector<TGraph>

        """
        
        # set up the waveform map the way AraRoot wants it
        # as a std::map<int, TGraph*>
        wf_map = ROOT.std.map("int", "TGraph*")()
        ROOT.SetOwnership(wf_map, True)
        for chan_i in waveform_bundle.keys():
            wf_map[chan_i] = waveform_bundle[chan_i]
                
        cross_correlations = self.rtc_wrapper.correlators["distant"].GetCorrFunctions(pairs, wf_map, applyHilbert)
        
        del wf_map

        return cross_correlations

    def do_standard_reco(self, wavepacket, applyHilbert=True):
        
        """
        A function to do a standard set of reconstructions.
        
        This is the "standard" reconstruction processing expected in the analysis.
        This means this function does a suite of things.
        The first thing it does is get the waveform map in the way AraRoot
        wants it (a std::map<int, TGraph*>).
        Then it gets the correlation functions between all the v antennas,
        and all the correlation functions between all the h antennas.
        It then calculates a suite of (currently) four maps, with the following names:
            "pulser_v" : vpol x-corr map at the radius of the cal pulser
            "pulser_h" : hpol x-corr map at the radius of the cal pulser
            "distant_v" : vpol x-corr map at the radius of 300 m
            "distant_h" : hpol x-corr map at the radius of 300 m 
        
        For each of these maps, it locates the peak corr and the phi/theta of the peak.
        For each of these maps, it returns a dictionary of results, with a key (a string),
        corresponding to the correlation, theta, phi, and the map itself.

        These four results each become their own entry in a master "reco_results"
        dictionary. See more information below on what the output looks like,
        and how to access elements of it.


        Parameters
        ----------
        wavepacket : dict
            A dict with three entries:
              "event" : int  
                Event number
              "waveforms" : dict
                A dictionary of the 16 waveforms.
                The key is the RF channel number.
                The value is a TGraph.
                There should be 16 entries, even if you don't intend to use all
                16 traces in your interferometry.
                The exclusions are handled further down under the excluded channels section.
              "trace_type" : string
                Waveform type requested by which_trace
        applyHilbert : bool
            Boolean whether the cross-correlation function is Hilbert enveloped.
 
        Returns
        -------
        reco_results : dict
            This is a dictionary of results.
            Each result is also a dictionary of results.
            (Dicts of dicts!)
            
            reco_results = {
                "pulser_v" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                              "which_distance": which_distance, "solution": solution, "polarization": polarization}
                "pulser_h" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                              "which_distance": which_distance, "solution": solution, "polarization": polarization}
                "distant_v_dir" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                                   "which_distance": which_distance, "solution": solution, "polarization": polarization}
                "distant_v_ref" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                                   "which_distance": which_distance, "solution": solution, "polarization": polarization}
                "distant_h_dir" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                                   "which_distance": which_distance, "solution": solution, "polarization": polarization}
                "distant_h_ref" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D, 
                                   "which_distance": which_distance, "solution": solution, "polarization": polarization}
            }

            So for example, if you want to find the peak theta value for the map
            calculated in vpol at the pulser distance, you should do:
                reco_results["pulser_v"]["theta"]
            
            The map you could get at by doing e.g.
                the_map = reco_results["pulser_v"]["map"]
            
        """

        reco_results = {}

        event_number = wavepacket["event"]
        waveform_bundle = wavepacket["waveforms"]

        # update correlation functions if needed
        if(self.__latest_event_num != event_number):
          self.__latest_event_num = event_number

        # get the correlation functions
        self.__corr_functions_v_fullphase = self.__calculate_cross_correlations(waveform_bundle, self.pairs_v, False)
        self.__corr_functions_h_fullphase = self.__calculate_cross_correlations(waveform_bundle, self.pairs_h, False)
        self.__corr_functions_v = self.__calculate_cross_correlations(waveform_bundle, self.pairs_v, True)
        self.__corr_functions_h = self.__calculate_cross_correlations(waveform_bundle, self.pairs_h, True) 

        # choose right correlation function for requested maps
        if(applyHilbert):
          this_corr_functions_v = self.__corr_functions_v
          this_corr_functions_h = self.__corr_functions_h
        else:
          this_corr_functions_v = self.__corr_functions_v_fullphase
          this_corr_functions_h = self.__corr_functions_h_fullphase
 
        ############################
        ####### VPol Pulser ########
        ############################

        # check the cal pulser in V
        pulser_map_v = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
            self.pairs_v, this_corr_functions_v, self.__arrival_delays_v_nearby, 0,)
        
        corr_pulser_v, phi_pulser_v, theta_pulser_v = mu.get_corr_map_peak(pulser_map_v)
        reco_results["pulser_v"] = {
            "corr" : corr_pulser_v,  "theta" : theta_pulser_v, "phi" : phi_pulser_v,
            "map" : pulser_map_v, "radius" : self.rtc_wrapper.correlators["nearby"].GetRadius(),
            'which_distance': 'nearby', 'solution': 0, 'polarization': 0,
        }

        ############################
        ####### VPol Maps ##########
        ############################

        # make a 300 m map in V (Direct rays)
        distant_map_v_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v, this_corr_functions_v, self.__arrival_delays_v_distant, 0,)

        # make a 300 m map in V (Refracted/Reflected rays)
        distant_map_v_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v, this_corr_functions_v, self.__arrival_delays_v_distant, 1,)

        # Get the correlation, phi, and theta for both maps
        corr_distant_v_dir, phi_distant_v_dir, theta_distant_v_dir = mu.get_corr_map_peak(distant_map_v_dir)
        corr_distant_v_ref, phi_distant_v_ref, theta_distant_v_ref = mu.get_corr_map_peak(distant_map_v_ref)

        # Store the direct rays results
        reco_results["distant_v_dir"] = {
            "corr": corr_distant_v_dir,  "theta": theta_distant_v_dir, "phi": phi_distant_v_dir,
            "map": distant_map_v_dir, "radius": self.rtc_wrapper.correlators["distant"].GetRadius(),
            'which_distance': 'distant', 'solution': 0, 'polarization': 0,
        }

        # Store the refracted/reflected rays results
        reco_results["distant_v_ref"] = {
            "corr": corr_distant_v_ref,  "theta": theta_distant_v_ref, "phi": phi_distant_v_ref,
            "map": distant_map_v_ref, "radius": self.rtc_wrapper.correlators["distant"].GetRadius(),
            'which_distance': 'distant', 'solution': 1, 'polarization': 0,
        }

        ############################
        ####### HPol Maps ##########
        ############################

        # make a 300 m map in H (Direct rays)
        distant_map_h_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h, this_corr_functions_h, self.__arrival_delays_h_distant, 0, )

        # make a 300 m map in H (Refracted/Reflected rays)
        distant_map_h_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h, this_corr_functions_h, self.__arrival_delays_h_distant, 1, )

        # Get the correlation, phi, and theta for both maps
        corr_distant_h_dir, phi_distant_h_dir, theta_distant_h_dir = mu.get_corr_map_peak(distant_map_h_dir)
        corr_distant_h_ref, phi_distant_h_ref, theta_distant_h_ref = mu.get_corr_map_peak(distant_map_h_ref)

        # Store the direct rays results
        reco_results["distant_h_dir"] = {
            "corr": corr_distant_h_dir,  "theta": theta_distant_h_dir, "phi": phi_distant_h_dir,
            "map": distant_map_v_dir, "radius": self.rtc_wrapper.correlators["distant"].GetRadius(),
            'which_distance': 'distant', 'solution': 0, 'polarization': 1,
        }

        # Store the refracted/reflected rays results in a separate dictionary
        reco_results["distant_h_ref"] = {
            "corr": corr_distant_h_ref,  "theta": theta_distant_h_ref, "phi": phi_distant_h_ref,
            "map": distant_map_v_ref, "radius": self.rtc_wrapper.correlators["distant"].GetRadius(),
            'which_distance': 'distant', 'solution': 1, 'polarization': 1,
        }

        return reco_results

    def lookup_arrival_time(self, 
                            channel : int, 
                            theta : float,
                            phi : float,
                            which_distance : str = "distant",
                            solution : int = 0
                            ):
        
        """
        A function to lookup the arrival times at a given antenna

        Parameters
        ----------
        channel : int
            The RF channel (0->16) about which you want information
        theta : float
            The theta value, in standard AraRoot RayTraceCorrelator interferometric coordinates, 
            meaning theta [-90 -> 90, 0 = horizontal],  of the sky point you want looked up.
        phi : float
            The phi value, in standard AraRoot RayTraceCorrelator interferometric coordinates, 
            meaning phi [-180 -> 180],  of the sky point you want looked up.
        which_distance : str
            Because this is our standard reco set, you can choose between "nearby" (at the cal pulser distance),
            or "distant" (at 300 m).
        solution : int
            Which solution do you want. 0 = direct, 1 = reflected/refracted.
        

        Returns
        -------
        arrival_time : float
            The arrival time at the antenna in nanoseconds, 
            assuming t=0 is when the ray leaves the vertex point of (theta, phi, R).
            This should be a positive number. A large negative number, e.g. -1000,
            indicates that no solution was found for that point on the sky.
        """

        if (channel < 0) or (not np.isfinite(channel)) or (channel > self.num_channels):
            raise ValueError(f"Channel number {channel} is not a valid request. Abort! ")
    
        if (abs(theta)>90) or (not np.isfinite(theta)):
            raise ValueError(f"Theta value requested ({theta}) is not physical. Make sure your number is in degrees, from -90 to 90")
        
        if (abs(phi)>180) or (not np.isfinite(phi)):
            raise ValueError(f"Phi value requested ({phi}) is not physical. Make sure your number is in degrees, from -180 to 180")
        
        if which_distance not in self.rtc_wrapper.correlators.keys():
            raise ValueError(f"Distance value requested ({which_distance}) is not loaded.")
        
        if solution not in [0, 1]:
            raise ValueError(f"Solution requested requested ({solution}) is not loaded.")

        arrival_time = self.rtc_wrapper.correlators[which_distance].LookupArrivalTime(
            channel, solution, theta, phi)

        return arrival_time

    def __get_correlation_function(self, ch1, ch2, wavepacket, applyHilbert):
        """
        Returns the correlation function for the channel pair ch1-ch2.

        Parameters
        ----------
        ch1 : int
          Number of first channel in correlation pair.
        ch2 : int
          Number of second channel in correlation pair.
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
        applyHilbert : bool
            Boolean whether the cross-correlation function is Hilbert enveloped.

        Returns
        -------
        corr_func : TGraph
          Cross-correlation function for requested channel pair.
        """

        event_number = wavepacket["event"]
        waveform_bundle = wavepacket["waveforms"]

        # update correlation functions if needed
        if(self.__latest_event_num != event_number):
          self.__latest_event_num = event_number

          # get the correlation functions
          self.__corr_functions_v_fullphase = self.__calculate_cross_correlations(waveform_bundle, self.pairs_v, False)
          self.__corr_functions_h_fullphase = self.__calculate_cross_correlations(waveform_bundle, self.pairs_h, False)
          self.__corr_functions_v = self.__calculate_cross_correlations(waveform_bundle, self.pairs_v, True)
          self.__corr_functions_h = self.__calculate_cross_correlations(waveform_bundle, self.pairs_h, True) 

        if((ch1 in const.vpol_channel_ids and ch2 in const.hpol_channel_ids) 
            or (ch1 in const.hpol_channel_ids and ch2 in const.vpol_channel_ids)):
          raise Exception("Correlation functions only available for like-polarization channels. Abort.")

        # Vpols
        if(ch1 in const.vpol_channel_ids):
          idx = self.get_pair_index(ch1, ch2, self.pairs_v)
          if(applyHilbert):
            return self.__corr_functions_v[idx]
          else:
            return self.__corr_functions_v_fullphase[idx]

        idx = self.get_pair_index(ch1, ch2, self.pairs_h)
        if(applyHilbert):
          return self.__corr_functions_h[idx]
        return self.__corr_functions_h_fullphase[idx]

    def get_pair_index(self, ch1, ch2, pairs):
        """
        Gets index of channel pair.

        Parameters
        ----------
        ch1 : int
          Number of first channel in pair.
        ch2 : int
          Number of second channel in pair.
        pairs : std::map<int, std::vector<int> >
          Channel pairs to get index from.

        Returns
        -------
        idx : int
          Index corresponding to channel pair.
        """

        for it in pairs:
          idx = it[0]
          c1 = it[1][0]
          c2 = it[1][1]

          # check both ways in case pairs aren't sorted
          if((ch1 == c1 and ch2 == c2) or (ch1 == c2 and ch2 == c1)):
            return idx

        raise Exception("Requested channel pair not found. Abort.")
        
        return -1 # useless but return anyway

    def get_surface_corr_max(self, corr_map, z_thresh=-10):
        """
        Calculates the maximum surface correlation value within a specified theta range.
        The theta range is determined based on the average antenna z-coordinate and a threshold.

        Parameters
        ----------
        corr_map : dict
            Dictionary containing:
            - 'map': ROOT.TH2D histogram with correlation values mapped by theta and phi.
            - 'corr': Overall peak correlation value.
            - 'radius': Correlation radius for the station (in meters).
        z_thresh : float, optional
            Depth threshold for defining the theta range lower bound. Default is -10 m.

        Returns
        -------
        max_surf_corr : float
            Maximum surface correlation within the theta range.
        max_theta : float
            Theta angle corresponding to max surface correlation.
        max_phi : float
            Phi angle corresponding to max surface correlation.
        """
        # Extract the radius and validate
        radius = corr_map.get("radius", None)
        if radius is None:
            raise ValueError("Radius not found in corr_map.")
        radius = float(radius)

        # Calculate average antenna z coordinates
        _, _, avg_z = mu.calculate_avg_antenna_xyz(self.station_id, self.num_channels)

        # Check if surface is visible at the given radius
        if radius < (abs(avg_z) + z_thresh):
            return -np.inf, 0, 0 # surface not visible, this map has no surface corr max

        # Define theta range -- interested in theta values above z_thresh, so values >= theta_thresh (using elevation angle!)
        theta_thresh = math.degrees(math.asin((abs(avg_z) + z_thresh) / radius))

        # Access correlation map and filter values within the theta range
        hist = corr_map.get("map", None)
        if hist is None:
            raise ValueError("The 'map' key was not found in corr_map.")
        
        # set the range of the histogram to only be where it's relevant
        binwidth = hist.GetYaxis().GetBinWidth()
        theta_thresh_bin = np.ceil((theta_thresh - hist.GetYaxis().GetXmin())/binwidth + 0.5) 
        hist.GetYaxis().SetRange(theta_thresh_bin, hist.GetYaxis().fNbins)
        
        # locate the peak in that range, and return the location
        _peakZ = ctypes.c_int()
        _peakTheta = ctypes.c_int()
        _peakPhi = ctypes.c_int()
        hist.GetMaximumBin(_peakPhi, _peakTheta, _peakZ)
        max_phi = hist.GetXaxis().GetBinCenter(_peakPhi.value)
        max_theta = hist.GetYaxis().GetBinCenter(_peakTheta.value)
        max_surf_corr = hist.GetMaximum()

        # reset the range
        hist.GetYaxis().SetRange(0,0)

        return max_surf_corr, max_theta, max_phi

    def get_surface_corr_max_multiple(self, reco_results, z_thresh=-10, skip_maps={}):
        """
        Applies get_surface_corr_max to multiple maps and identifies the map with
        the highest surface correlation value.

        Parameters
        ----------
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        z_thresh : float, optional
            Depth threshold to define theta range lower bound. Default is -10 m.
        skip_maps : set
            Will skip every reconstruction result key listed in this set. 

        Returns
        -------
        results : list of dicts
            List with results for each map:
            - 'max_corr': Maximum surface correlation.
            - 'theta': Theta angle of max surface correlation.
            - 'phi': Phi angle of max surface correlation.
        max_result : dict
            Dictionary with the max surface correlation result across all maps.

        Usage example
        -------

        results, max_surf_result = standard_reco.get_surface_corr_max_multiple(
            reco_results["distant_v_dir"], reco_results["distant_v_ref"], z_thresh=-10
        )

        """
        results = []
        
        skip_maps = set(skip_maps) # make a set to speed up compilation

        for idx, (name, corr_map) in enumerate(reco_results.items()):
            if name in skip_maps: 
                continue
            max_corr, theta, phi = self.get_surface_corr_max(corr_map, z_thresh=z_thresh)
            results.append({
                'max_corr': max_corr, 'theta': theta, 'phi': phi, 'map index': idx, 
                'name': name, 'which_distance': corr_map['which_distance'],
                'solution': corr_map['solution'], 'polarization': corr_map['polarization']})

        valid_results = [res for res in results if res is not None]
        if valid_results:
            max_result = max(valid_results, key=lambda x: x['max_corr'])
        else:
            raise RuntimeError("No valid results to determine max surface correlation.")
        
        return {result['name']: result for result in results if result is not None}, max_result

    def find_map_with_max_corr(self, reco_results, skip_maps={}):
        """
        Finds the map with the highest overall 'corr' value.

        Parameters
        ----------
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        skip_maps : set
            Will skip every reconstruction result key listed in this set. 

        Returns
        -------
        max_corr_result : dict
            Dictionary with max 'corr' value and corresponding 'theta', 'phi'.

        Usage example:
        -------

        max_corr_info = standard_reco.find_map_with_max_corr(
            reco_results["distant_v_dir"], reco_results["distant_v_ref"]
        )

        """
        max_corr_result = {'max_corr': -float('inf')}
        
        skip_maps = set(skip_maps) # make a set to speed up compilation

        for idx, (name, corr_map) in enumerate(reco_results.items()):
            if name in skip_maps: 
                continue

            corr = corr_map['corr']
            if corr > max_corr_result['max_corr']:
                max_corr_result.update({
                    'max_corr': corr, 'theta': corr_map['theta'], 'phi': corr_map['phi'], 
                    'map_index': idx, 'name': name, 'which_distance': corr_map['which_distance'],
                    'solution': corr_map['solution'], 'polarization': corr_map['polarization']})
                    
        if max_corr_result['max_corr'] == -float('inf'):
            raise RuntimeError("No valid correlation values found.")
        
        return max_corr_result


    def calculate_surface_corr_ratio(self, reco_results, z_thresh=-10, skip_maps={}):
        """
        Calculates surface correlation ratio by comparing max surface correlation with max overall correlation.

        Parameters
        ----------
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        z_thresh : float, optional
            Depth threshold to define theta range lower bound. Default is -10 m.
        skip_maps : set
            Will skip every reconstruction result key listed in this set. 

        Returns
        -------
        surf_corr_ratio : float
            Ratio of max surface correlation to max overall correlation.
        max_surf_corr_result : dict
            Result for map with max surface correlation.
        max_corr_result : dict
            Result for map with max overall correlation.

        Usage example
        -------
        surf_corr_ratio, max_surf_corr_result, max_corr_result = standard_reco.calculate_surface_corr_ratio(
            reco_results["distant_v_dir"], reco_results["distant_v_ref"], z_thresh=-10
        )
        
        """
        _, max_surf_corr_result = self.get_surface_corr_max_multiple(
            reco_results, z_thresh=z_thresh, skip_maps=skip_maps)
        max_corr_result = self.find_map_with_max_corr(reco_results, skip_maps=skip_maps)
        
        if max_corr_result['max_corr'] != 0:
          surf_corr_ratio = max_surf_corr_result['max_corr'] / max_corr_result['max_corr']  
        else:
          surf_corr_ratio = np.inf
        
        return surf_corr_ratio, max_surf_corr_result, max_corr_result

    def min_frac_corr_depth(self, corr_map, max_corr, fraction=0.6):
        """
        Finds the shallowest depth where correlation meets a specified fraction of the maximum correlation,
        using 'theta' to calculate depth and adjusting by the average z-coordinate of antennas.

        Parameters
        ----------
        corr_map : dict
            Contains:
            - 'map': ROOT.TH2D histogram with correlation values by theta and phi.
            - 'radius': Correlation radius in meters.
        max_corr : float
            The maximum correlation to use in setting the correlation threshold.
        fraction : float
            Fraction of the max correlation to set the threshold.
            Default is 0.6 (or 60%).

        Returns
        -------
        min_depth : float
            Shallowest depth (measured from ice surface) where correlation meets the specified fraction of max.
        """
        # Get correlation threshold
        threshold_corr = fraction * max_corr

        # Access correlation map and radius
        hist = corr_map.get("map")
        radius = corr_map.get("radius")
        if hist is None or radius is None:
            raise ValueError("Map or radius not found in corr_map.")
        radius = float(radius)

        # Calculate average antenna z-coordinate
        _, _, avg_z = mu.calculate_avg_antenna_xyz(self.station_id, self.num_channels)

        # Find the shallowest bin above threshold.
        # We want the last bin ABOVE because theta goes from -90 to 90,
        # and we want the SHALLOWEST (so closets to 90).
        last_bin = hist.FindLastBinAbove(threshold_corr, 2) # look along the Y (theta) axis, so axis=2
        if last_bin != -1:
            # -1 to ROOT means no bin was found above the threshold
            theta = hist.GetYaxis().GetBinCenter(last_bin)
            min_depth = radius * math.sin(math.radians(theta)) + avg_z
        else:
            min_depth = -float('inf')  

        return min_depth

    def min_frac_corr_depth_multiple(self, reco_results, fraction=0.6, skip_maps={}):
        """
        Finds the shallowest depth across multiple correlation maps where correlation meets
        a specified fraction of the maximum correlation. Uses the min_frac_corr_depth function for each map,
        then returns the minimum depth found and the index of the corresponding map.

        Parameters
        ----------
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        fraction : float
            Fraction of the max correlation to set the threshold. Default is 0.6 (60%).
        skip_maps : set
            Will skip every reconstruction result key listed in this set. 

        Returns
        -------
        min_depth : float
            Shallowest depth (measured from ice surface) across all maps where correlation meets the specified fraction of max.
        min_depth_index : int
            Index of the map that contains the shallowest depth.
        """
        min_depth = -float('inf')  # Initialize to find the shallowest depth (least negative)
        min_depth_index = -1  # Track the index of the map with the minimum depth

        # Make a set to speed up compilation
        skip_maps = set(skip_maps)

        # Get max correlation across all maps
        max_corr = self.find_map_with_max_corr(reco_results, skip_maps=skip_maps)['max_corr']

        # Iterate over each map and find the shallowest depth using min_frac_corr_depth
        for idx, (name, corr_map) in enumerate(reco_results.items()):
            if name in skip_maps: 
                continue
            # Use the existing function to find the minimum depth for this map
            depth = self.min_frac_corr_depth(corr_map, max_corr, fraction=fraction)


            # Update min_depth if the current depth is shallower (closer to the surface)
            if depth > min_depth:
                min_depth = depth
                min_depth_index = idx

        if min_depth == -float('inf'):
            raise RuntimeError("No maps meet the fractional correlation threshold within the z_thresh.")

        return min_depth, min_depth_index

    def get_corr_snr(self, ch1, ch2, wave_packet):

        """
        Calculates channel-pair correlation SNR
    
        Parameters
        ----------
        ch1: int
            First channel number for correlation.
        ch2: int
            Second channel number for correlation.
        wavepacket: dict
            A dict with three entries:
              "event": int  
                Event number
              "waveforms": dict
                Dict mapping RF channel ID to waveforms.
                Keys are channel id (an integer)
                Values are TGraphs
              "trace_type": string
                Waveform type requested by which_trace.
 
        Returns
        -------
        corr_snr: float
           Channel pair correlation SNR.
        """

        corr_func = self.__get_correlation_function(ch1, ch2, wave_packet, False)
        corr_snr = snr.get_snr(corr_func)

        return corr_snr

    def get_avg_corr_snr(self, wave_packet, excluded_channels = []):

        """
        Computes channel-pair average correlation SNR for VPols and HPols.

        Parameters
        ----------
        wavepacket: dict
            A dict with three entries:
              "event": int  
                Event number
              "waveforms": dict
                Dict mapping RF channel ID to waveforms.
                Keys are channel id (an integer)
                Values are TGraphs
              "trace_type": string
                Waveform type requested by which_trace.
        excluded_channels: list
            List of dictionary keys to exclude from average.

        Returns
        -------
        avg_vpol_corr_snr : float
            The VPol average channel-pair correlation SNR.
        avg_hpol_corr_snr : float
            The HPol average channel-pair correlation SNR.
        """

        wave_bundle = wave_packet["waveforms"]
        chans = sorted(list(wave_bundle.keys()))

        good_vpol_chs, good_hpol_chs = [], []

        for chan in chans:
            if chan in excluded_channels:
                continue

            if chan in const.vpol_channel_ids:
                good_vpol_chs.append(chan)
            if chan in const.hpol_channel_ids:
                good_hpol_chs.append(chan)

        avg_vpol_corr_snr, avg_hpol_corr_snr = [], []
        
        for pair_V in list(itertools.combinations(good_vpol_chs, 2)):
            corr_snr_V = self.get_corr_snr(pair_V[0], pair_V[1], wave_packet)
            avg_vpol_corr_snr.append(corr_snr_V)

        for pair_H in list(itertools.combinations(good_hpol_chs, 2)):
            corr_snr_H = self.get_corr_snr(pair_H[0], pair_H[1], wave_packet)
            avg_hpol_corr_snr.append(corr_snr_H)

        avg_vpol_corr_snr = np.mean(avg_vpol_corr_snr)
        avg_hpol_corr_snr = np.mean(avg_hpol_corr_snr)

        return avg_vpol_corr_snr, avg_hpol_corr_snr  

    def get_arrival_delays_reco(
        self, reco_results, excluded_channels, reference_ch, 
        which_distance, solution
    ):
        """
        Extract the arrival time of the signals from the reconstruction then 
        calculate the time delays between each channel and the reference channel

        Parameters
        ----------
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        excluded_channels : set
            A set of channels IDs to exclude from arrival time calcualtion
        reference_ch : int
            Reference channel for use calculating arrival delays
        which_distance : str
            The distance of the reconstruction for analysis, following the 
            convention of other functions in this class. 
        solution : int
            0 for direct ray solution. 1 for reflected/refracted ray solution.

        Returns
        -------
        arrival_delays : dict
            Dictionary where keys are channel numbers and values are arrival 
            delays calculated based on expected arrival times from the 
            reconstructed event vertex.
        """

        # Make the excluded_channels object a set to speed up computation
        excluded_channels = set(excluded_channels)

        # Get all the arrival times calculated based on the reconstructed
        #   event vertex and the time tables
        arrival_times = {
            ch: 0.0 for ch in const.rf_channels_ids if ch not in excluded_channels}
        for ch_ID in arrival_times.keys():
            arrival_times[ch_ID] = self.lookup_arrival_time(
                ch_ID, reco_results['theta'], reco_results['phi'], 
                which_distance=which_distance, solution=solution
            )
        
        # Calculate time difference between each channel's arrival time
        #   and the reference channel's arrival time
        reference_arrival_time = arrival_times[reference_ch]
        arrival_delays = {ch: 0.0 for ch in arrival_times.keys()}
        for ch_ID in arrival_delays.keys():
            arrival_delays[ch_ID] = arrival_times[ch_ID] - reference_arrival_time

        return arrival_delays

    def get_arrival_delays_AraRoot_xcorr(
        self, wavepacket, excluded_channels, reference_ch, reco_delays, 
        is_software, zoom_window=40
    ):
        """
        Determine the arrival delays from each channel by finding the time of
        max correlation between each channel and the reference channel in a 
        `zoom_window` nanosecond window around the expected arrival delay of a 
        signal given the results from reconstructing the event.

        Parameters
        ----------
        wavepacket : dict
            AraProc wavepacket (keys include waveforms and the event number)
        excluded_channels : set
            A set of channels IDs to exclude from arrival time calcualtion
        reference_ch : int
            Reference channel for use calculating arrival delays
        reco_delays : dict
            Dictionary where keys are channel numbers and values are arrival 
            delays calculated based on expected arrival times from the 
            reconstructed event vertex.
        is_software : bool
            Boolean if the event is a software trigger or not. 
        zoom_window : float
            Total size of the window (in nanoseconds) we will use to find the 
            maximum cross correlation value around the expected arrival delay 
            from the reconstructed event vertex.
        
        Returns
        -------
        delays : dict
            Dictionary where keys are channel numbers and values are arrival
            delays calculated from cross correlating each channel with 
            the reference channel.
        """

        # Make the excluded_channels object a set to speed up computation
        excluded_channels = set(excluded_channels)

        # Calculate the time delay between each channel and the reference channel
        delays = {}
        for ch_ID in wavepacket['waveforms'].keys():

            # If this channel is to be excluded, skip it.
            if ch_ID in excluded_channels: 
                continue

            if ch_ID == reference_ch: 
                # Delay between the reference channel and itself will be 0
                delay = 0
            else: 
                
                # Load  the cross correlation for this channel and the reference 
                xcorr_times, xcorr_volts = wfu.tgraph_to_arrays(
                    self.__get_correlation_function(
                        ch_ID, reference_ch, wavepacket, applyHilbert=False
                    )
                )

                # Identify the `zoom_window` nanosecond window around the 
                #   reconstructed delay
                zoomed_indices = np.where(
                    (np.abs( xcorr_times - reco_delays[ch_ID] )) < zoom_window // 2
                )[0]

                # Calculate the time of maximum correlation from this
                #   window of expected signal delay.
                if len(zoomed_indices) == 0: 

                    # Software triggers are so short, a channel may not have
                    #   signal during the time window where signal is expected.
                    #   As a result, `zoom_indices` will have no entries. 
                    # Just find the time of peak correlation for the full 
                    #   waveform in this scenario.
                    delay = xcorr_times[ np.argmax(xcorr_volts) ]

                    # If the event can't find the time window to zoom in on
                    #   and its not a software trigger, warn user
                    if not is_software: 
                        print(
                            f"Trouble calculating csw for non-software event with "
                            f"channels {ch_ID} and {reference_ch}")

                else: 
                    delay = xcorr_times[ 
                        np.argmax(xcorr_volts[zoomed_indices]) # index of max xcorr in zoomed array
                        + zoomed_indices[0] # Adjusted by first zoomed_index
                    ]

            # AraRoot always compares channels with smaller IDs to channels with 
            #   larger IDs but we always want to compare to the reference channel.
            # Correct for this if the current ch_ID is larger than the reference ch_ID
            if ch_ID > reference_ch: 
                delay *= -1
            
            # Save the calculated arrival delay
            delays[ch_ID] = delay

        return delays

    def get_csw(
        self, wavepacket, is_software, is_calpulser, solution, polarization, reco_results,
        excluded_channels, which_distance='distant'
    ):
        """
        Build the Coherently Summed Waveform (CSW) for the given `useful_event`
          by shifting each waveform accoriding to the time delay calculated 
          after finding the max cross correlation between two channels 
          in a window around the expected arrival delay. 

        Parameters
        ----------
        data : AnalysisDataset
        wavepacket : dict
            AraProc wavepacket (keys include waveforms and the event number)
        is_software : bool
            Boolean for whether this event is a software trigger or not.
        solution : int
            0 for direct ray solution. 1 for reflected/refracted ray solution.
        polarization : int
            0 for VPol, 1 for HPol
        reco_results : dict
            Reco results already with the specific reconstruction reqeusted (so
            the keys of this object include 'theta' and 'phi')
        excluded_channels : set
            Any channels to exclude that aren't already excluded by 
            livetime configurations.
        which_distance : str
            The distance of the reconstruction for analysis, following the 
            convention of other functions in this class. 

        Returns
        -------
        csw : ROOT.TGraph
            The coherently summed waveform. Only returned if `return_csw==True`
        warning : int
            If non-zero, contains information about something that went wrong
              in the CSW calculation.
        """

        # Initialize warning object as 0
        warning = 0

        # Check that reco_results have been passed in properly
        if 'theta' not in reco_results.keys():
            raise ValueError(
                "reco_result object passed into this function does not appear "
                "to be for a specific reconstruction. Should have a 'theta' key "
                "but only has the following:", reco_results.keys())

        # Make the excluded_channels object a set to speed up computation
        excluded_channels = set(excluded_channels)

        # Determine which channels should be used in the csw
        channels_to_csw = []
        for ch_ID in wavepacket['waveforms'].keys(): 
            if ch_ID not in excluded_channels:
                channels_to_csw.append(ch_ID)

        # Determine the "reference channel" to base the CSW around as the
        #   channel with the maximum voltage
        reference_ch = -123456
        reference_ch_max_voltage = -1
        for ch_ID in channels_to_csw:
            this_max_voltage = np.max(np.asarray(wavepacket['waveforms'][ch_ID].GetY()))
            if this_max_voltage > reference_ch_max_voltage:
                reference_ch_max_voltage = this_max_voltage
                reference_ch = ch_ID 

        # Get arrival delays relative to the reference channel based on
        #   expected arrival times from reconstruction results
        arrival_delays_reco = self.get_arrival_delays_reco(
            reco_results, excluded_channels, reference_ch, which_distance, solution)
        
        # Get arrival delays by zooming in on the cross correlation between
        #   each channel and the reference channel around the expected
        #   arrival delay from reconstruction results
        arrival_delays = self.get_arrival_delays_AraRoot_xcorr(
            wavepacket, excluded_channels, reference_ch, 
            arrival_delays_reco, is_software)  
        
        # Calculate expected signal time as when the reference channel saw 
        #   its maximum signal
        expected_signal_time = np.asarray(wavepacket['waveforms'][reference_ch].GetX())[
            np.argmax( np.asarray(wavepacket['waveforms'][reference_ch].GetY()) )
        ]

        # Initialize the final CSW waveform time and voltage arrays using the
        #   reference channel's time array resized to size of the channel with 
        #   the shortest waveform's waveform
        shortest_wf_ch = 123456
        shortest_wf_length = np.inf
        for ch_ID in channels_to_csw:
            if wavepacket['waveforms'][ch_ID].GetN() < shortest_wf_length: 
                shortest_wf_length = wavepacket['waveforms'][ch_ID].GetN() 
                shortest_wf_ch = ch_ID
        csw_values = np.zeros((1, shortest_wf_length))
        csw_times = np.asarray(
            wavepacket['waveforms'][reference_ch].GetX())[:shortest_wf_length]
        csw_dt = csw_times[1] - csw_times[0]

        # Roll the waveform from each channel so the starting time of each
        #   waveform lines up. Then add the waveform to the CSW.
        for ch_ID in channels_to_csw: 

            # Load this channel's voltage and time arrays. Shift time by arrival delay
            values = np.asarray(wavepacket['waveforms'][ch_ID].GetY())
            times = np.asarray(wavepacket['waveforms'][ch_ID].GetX()) - arrival_delays[ch_ID]

            # Determine if the time binning of this channel matches the time binning
            #   of the CSW. This assumes all waveforms have the same time
            #   bin size. Calculating arrival delays via the cross correlation
            #   should fix this already, so just raise a warning if
            #   the bins aren't lined up.
            # Ex: If the CSW has times of [0.9, 1.4, 1.9], waveform 1 has a
            #   time array of [3.4, 3.9, 4.4 ], waveform 2 has a time array 
            #   of [1.6, 2.1, 2.6], and waveform 3 has a time array of 
            #   [0.5, 1.0, 1.5] then waveform 1 doesn't need rebinning since 
            #   each time step starts and ends on the same fraction of a time 
            #   bin as the CSW does. However, waveform 2 needs to be shifted back
            #   by 0.2 nanoseconds to get a fitting time aray of [1.4, 1.9, 2.4]
            #   and waveform 3 needes to be shifted back by 0.1
            rebinning_shift = (
                (csw_times[0] - times[0])
                % csw_dt
                # Take the remainder of the start time difference with csw_dt.
                #   If this is ~0, the waveforms have the same binning 
                #   otherwise, this reveals how much to shift the waveform by. 
                # For example, waveform 1 yeilds: (3.4 - 0.9) % 0.5 = 0
                #   waveform 2 yeilds (1.6 - 1.4) % 0.5 = 0.2 
                #   and waveform 3 yeilds (0.9 - 0.5) % 0.5 = 0.4 
            ) 
            if csw_dt - 0.0001 > abs(rebinning_shift) > 0.0001: 
                warning += 10_00_00

            # Trim this waveform's length to match the CSW length
            if len(times) > len(csw_times):
                trim_ammount = len(times) - len(csw_times)
                if (
                    ( times[0] - csw_times[0] < 0 ) # this wf has an earlier start time than the CSW
                    and ( times[-1] - csw_times[-1] <= csw_dt/2) # this wf has a earlier or equal end time than the CSW
                ): # We need to trim from the beginning of the waveform
                    times  = times [trim_ammount:]
                    values = values[trim_ammount:]
                elif (
                    ( times[0] - csw_times[0] > -csw_dt/2) # this wf has a later or equal start time than the CSW
                    and (times[-1] - csw_times[-1] > 0) # this wf has a later end time than the CSW
                ): # we need to trim from the end of the waveform
                    times  = times [:-trim_ammount]
                    values = values[:-trim_ammount]
                elif (
                    ( times[0] - csw_times[0] < 0 ) # this wf starts earlier than the CSW 
                    and ( times[-1] - csw_times[-1] > 0 ) # this wf ends later than the CSW
                ): # we need to trim from both ends of the waveform
                    leading_trimmable = np.argwhere( 
                        np.round(times,5) < np.round(csw_times[0], 5) )
                    trailing_trimmable = np.argwhere( 
                        np.round(times, 5) > np.round(csw_times[-1], 5) )
                    times  = times [ len(leading_trimmable) : -len(trailing_trimmable) ] 
                    values = values[ len(leading_trimmable) : -len(trailing_trimmable) ] 
                else: 
                    # This waveform either has the same start and end time and 
                    #   has no need to be trimmed, or starts and ends within
                    #   the CSW time array, so it should be the shortest
                    #   waveform unless the binning changed between waveforms.
                    raise ValueError(
                        f"Unexpected trimming request for channel {ch_ID}. "
                        f"The CSW has a length of {len(csw_times)} and a dt of {csw_dt} "
                        f"which may be incompatible with channel {ch_ID}s waveform "
                        f"with a length of {len(times)} and a dt of {times[1]-times[0]}. "
                        f"The CSW runs from {csw_times[0]} to {csw_times[-1]} and "
                        f"this waveform runs from {times[0]} to {times[-1]}." 
                    )

            # Roll the waveform so that the start and end times of the waveform
            #   line up exactly with the CSW
            roll_shift_bins = (csw_times[0] - times[0]) / csw_dt
            roll_shift_time = roll_shift_bins*(times[1] - times[0])
            if abs(roll_shift_bins) % 1.0 > 0.0001: 
                # roll_shift is not close to an integer. Add to the warning
                warning += 10
            roll_shift_bins = int(roll_shift_bins)
            if abs(roll_shift_bins)>len(times) and not is_software:
                # More waveform to roll than there is time in the waveform,
                #   so add to the warning tracker. 
                # Software triggers are so short, this sometimes occurs for them.
                #   Don't warn in this scenario.
                warning += 10_00_00_00
            if roll_shift_bins > 0 and abs(roll_shift_bins)<len(times): 
                # Rolling from front to back, check that signal region isn't in the front
                if times[0] <= expected_signal_time <= times[roll_shift_bins]: 
                    warning += 10_00
            elif roll_shift_bins < 0 and abs(roll_shift_bins)<len(times): 
                # Rolling from back to front, check that signal region isn't in the back
                if  times[roll_shift_bins]  <= expected_signal_time <= times[-1]: 
                    warning += 10_00
            rolled_wf = np.roll( values, -roll_shift_bins )
            rolled_times = np.linspace(
                times[0] + roll_shift_time,
                times[-1] + roll_shift_time,
                len(times)
            )

            # Add this channel's waveform to the CSW
            csw_values = np.sum( np.dstack( (csw_values[0], rolled_wf) ), axis=2) 

        # Un-nest the csw. csw.shape was (1,len(csw_times)) but is now len(csw_times)
        csw_values = np.squeeze(csw_values)
        
        return wfu.arrays_to_tgraph(csw_times, csw_values), warning



