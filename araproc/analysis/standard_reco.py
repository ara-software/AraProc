import ctypes
import numpy as np
import math
import os
import ROOT

from araproc.analysis import interferometry as interf
from araproc.framework import constants as const
from araproc.framework import map_utilities as mu


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
        dir_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
                                "PUBLIC/groups/arasoft/raytrace_timing_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_{self.calpulser_r_library[station_id]}_angle_1.00_solution_0.root"
                                )
        ref_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
                                "PUBLIC/groups/arasoft/raytrace_timing_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_{self.calpulser_r_library[station_id]}_angle_1.00_solution_1.root"
                                )
        self.rtc_wrapper.add_rtc(ref_name = "nearby",
                radius=float(self.calpulser_r_library[station_id]),
                path_to_dir_file=dir_path,
                path_to_ref_file=ref_path
                )

        # always add a "distant" correlator for 300m away
        dir_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
                                "PUBLIC/groups/arasoft/raytrace_timing_tables",
                                f"arrivaltimes_station_{self.station_id}_icemodel_40_radius_300.00_angle_1.00_solution_0.root"
                                )
        ref_path = os.path.join("/cvmfs/icecube.osgstorage.org/icecube",
                                "PUBLIC/groups/arasoft/raytrace_timing_tables",
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

        self.__latest_event_num = -1

    def __calculate_cross_correlations(self, waveform_bundle, pairs):

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
 
        Returns
        -------
        cross_correlations : std::vector<TGraph>

        """
        
        # set up the waveform map the way AraRoot wants it
        # as a std::map<int, TGraph*>
        waveform_map = ROOT.std.map("int", "TGraph*")()
        ROOT.SetOwnership(waveform_map, True)
        for chan_i in waveform_bundle.keys():
            waveform_map[chan_i] = waveform_bundle[chan_i]
                
        cross_correlations = self.rtc_wrapper.correlators["distant"].GetCorrFunctions(pairs, wf_map)
        
        del waveform_map

        return cross_correlations

    def do_standard_reco(self, waveform_bundle, event_number):
        
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
        waveform_bundle : dict
            A dictionary of the 16 waveforms.
            The key is the RF channel number.
            The value is a TGraph.
            There should be 16 entries, even if you don't intend to use all
            16 traces in your interferometry.
            The exclusions are handled further down under the excluded channels section.
        event_number : int
            The number of the event corresponding to the traces in wave_bundle. Required
            to access event-specific information like the correlation functions.       
 
        Returns
        -------
        reco_results : dict
            This is a dictionary of results.
            Each result is also a dictionary of results.
            (Dicts of dicts!)
            
            reco_results = {
                "pulser_v" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D}
                "pulser_h" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D}
                "distant_v" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D}
                "distant_h" : {"corr" : peak_corr, "theta" : peak_theta, "phi" : peak_phi, "map" : ROOT.TH2D}
            }

            So for example, if you want to find the peak theta value for the map
            calculated in vpol at the pulser distance, you should do:
                reco_results["pulser_v"]["theta"]
            
            The map you could get at by doing e.g.
                the_map = reco_results["pulser_v"]["map"]
            
        """

        reco_results = {}

        # update correlation functions if needed
        if(self.__latest_event_num != event_number):
          self.__latest_event_num = event_number

          # get the correlation functions
          self.__corr_functions_v = __calculate_cross_correlations(waveform_bundle, self.pairs_v)
          self.__corr_functions_h = __calculate_cross_correlations(waveform_bundle, self.pairs_h) 
        
        # check the cal pulser in V
        pulser_map_v = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
            self.pairs_v,
            self.__corr_functions_v,
            0
        )
        corr_pulser_v, phi_pulser_v, theta_pulser_v = mu.get_corr_map_peak(pulser_map_v)
        reco_results["pulser_v"] = {"corr" : corr_pulser_v, 
                                    "theta" : theta_pulser_v,
                                    "phi" : phi_pulser_v,
                                    "map" : pulser_map_v,
                                    "radius" : self.rtc_wrapper.correlators["nearby"].GetRadius()
                                    }

        # # check the cal pulser in H
        # pulser_map_h = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
        #     self.pairs_h,
        #     self.__corr_functions_h,
        #     0
        # )
        # corr_pulser_h, phi_pulser_h, theta_pulser_h = mu.get_corr_map_peak(pulser_map_h)
        # reco_results["pulser_h"] = {"corr" : corr_pulser_h, 
        #                             "theta" : theta_pulser_h,
        #                             "phi" : phi_pulser_h,
        #                             "map" : pulser_map_h,
        #                             "radius" : self.rtc_wrapper.correlators["nearby"].GetRadius()
        #                             }

        # make a 300 m map in V (Direct rays)
        distant_map_v_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            self.__corr_functions_v,
            0
        )

        # make a 300 m map in V (Refracted/Reflected rays)
        distant_map_v_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            self.__corr_functions_v,
            1
        )

        # Get the correlation, phi, and theta for both maps
        corr_distant_v_dir, phi_distant_v_dir, theta_distant_v_dir = mu.get_corr_map_peak(distant_map_v_dir)
        corr_distant_v_ref, phi_distant_v_ref, theta_distant_v_ref = mu.get_corr_map_peak(distant_map_v_ref)

        # Store the direct rays results
        reco_results["distant_v_dir"] = {
            "corr": corr_distant_v_dir, 
            "theta": theta_distant_v_dir,
            "phi": phi_distant_v_dir,
            "map": distant_map_v_dir,
            "radius": self.rtc_wrapper.correlators["distant"].GetRadius()
        }

        # Store the refracted/reflected rays results
        reco_results["distant_v_ref"] = {
            "corr": corr_distant_v_ref, 
            "theta": theta_distant_v_ref,
            "phi": phi_distant_v_ref,
            "map": distant_map_v_ref,
            "radius": self.rtc_wrapper.correlators["distant"].GetRadius()
        }

        # make a 300 m map in H (Direct rays)
        distant_map_h_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h,
            self.__corr_functions_h,
            0
        )

        # make a 300 m map in H (Refracted/Reflected rays)
        distant_map_h_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h,
            self.__corr_functions_h,
            1
        )

        # Get the correlation, phi, and theta for both maps
        corr_distant_h_dir, phi_distant_h_dir, theta_distant_h_dir = mu.get_corr_map_peak(distant_map_h_dir)
        corr_distant_h_ref, phi_distant_h_ref, theta_distant_h_ref = mu.get_corr_map_peak(distant_map_h_ref)

        # Store the direct rays results
        reco_results["distant_h_dir"] = {
            "corr": corr_distant_h_dir, 
            "theta": theta_distant_h_dir,
            "phi": phi_distant_h_dir,
            "map": distant_map_v_dir,
            "radius": self.rtc_wrapper.correlators["distant"].GetRadius()
        }

        # Store the refracted/reflected rays results in a separate dictionary
        reco_results["distant_h_ref"] = {
            "corr": corr_distant_h_ref, 
            "theta": theta_distant_h_ref,
            "phi": phi_distant_h_ref,
            "map": distant_map_v_ref,
            "radius": self.rtc_wrapper.correlators["distant"].GetRadius()
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

        theta_bin = ctypes.c_int()
        phi_bin = ctypes.c_int()
        self.rtc_wrapper.correlators[which_distance].ConvertAngleToBins(theta, phi,
                                                                        theta_bin, phi_bin)
        arrival_time = self.rtc_wrapper.correlators[which_distance].LookupArrivalTimes(channel,
                                                                                       solution,
                                                                                       theta_bin.value,
                                                                                       phi_bin.value
                                                                                       )
        return arrival_time

    def __get_correlation_function(self, ch1, ch2, event_number):
        """
        Returns the correlation function for the channel pair ch1-ch2.

        Parameters
        ----------
        ch1 : int
          Number of first channel in correlation pair.
        ch2 : int
          Number of second channel in correlation pair.
        event_number : int
          Number of event whose correlation function is being requested.

        Returns
        -------
        corr_func : TGraph
          Cross-correlation function for requested channel pair.
        """

        # update correlation functions if needed
        if(self.__latest_event_num != event_number):
          self.__latest_event_num = event_number

          # get the correlation functions
          self.__corr_functions_v = __calculate_cross_correlations(waveform_bundle, self.pairs_v)
          self.__corr_functions_h = __calculate_cross_correlations(waveform_bundle, self.pairs_h) 

        if(ch1//8 != ch2//8):
          raise Exception("Correlation functions only available for like-polarization channels. Abort.")

        # Vpols
        if(ch1//8 == 0):
          idx = get_pair_index(ch1, ch2, self.pairs_v)
          return self.__corr_functions_v[idx]

        idx = get_pair_index(ch1, ch2, self.pairs_h)
        return self.__corr_functions_h[idx]

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
            raise RuntimeError("Surface not visible: radius is smaller than average antenna z-coordinate.")

        # Define theta range
        theta_surface = math.degrees(math.asin(abs(avg_z) / radius))
        theta_thresh = math.degrees(math.asin((abs(avg_z) + z_thresh) / radius))
        theta_min, theta_max = min(theta_surface, theta_thresh), max(theta_surface, theta_thresh)

        # Access correlation map and filter values within the theta range
        hist = corr_map.get("map", None)
        if hist is None:
            raise ValueError("The 'map' key was not found in corr_map.")
        
        surf_corr_values = []
        corr_bins = []

        # Loop through the bins and store values within theta range
        for x_bin in range(1, hist.GetNbinsX() + 1):
            for y_bin in range(1, hist.GetNbinsY() + 1):
                theta = hist.GetYaxis().GetBinCenter(y_bin)
                phi = hist.GetXaxis().GetBinCenter(x_bin)
                corr = hist.GetBinContent(x_bin, y_bin)
                if theta_min <= theta <= theta_max:
                    surf_corr_values.append(corr)
                    corr_bins.append((x_bin, y_bin))

        # Identify maximum surface correlation
        if not surf_corr_values:
            raise RuntimeError("No correlation values found in the specified theta range.")
        
        # Check for NaNs
        if np.isnan(surf_corr_values).any():
            raise ValueError("surf_corr_values contains NaN values.")

        # Calculate maximum and index
        max_surf_corr = np.max(surf_corr_values)
        max_idx = np.argmax(surf_corr_values)
        max_bin = corr_bins[max_idx]
        
        max_theta = hist.GetYaxis().GetBinCenter(max_bin[1])
        max_phi = hist.GetXaxis().GetBinCenter(max_bin[0])

        return max_surf_corr, max_theta, max_phi

    def get_surface_corr_max_multiple(self, *maps, z_thresh=-10):
        """
        Applies get_surface_corr_max to multiple maps and identifies the map with
        the highest surface correlation value.

        Parameters
        ----------
        *maps : dict
            Variable number of dictionaries, each containing:
            - 'map': ROOT.TH2D histogram with correlation values.
            - 'corr': Overall peak correlation value.
            - 'radius': Correlation radius for the station.
        z_thresh : float, optional
            Depth threshold to define theta range lower bound. Default is -10 m.

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
        for idx, corr_map in enumerate(maps):
            try:
                max_corr, theta, phi = self.get_surface_corr_max(corr_map, z_thresh=z_thresh)
                results.append({'max_corr': max_corr, 'theta': theta, 'phi': phi, 'map index': idx})
            except (ValueError, RuntimeError) as e:
                print(f"Error processing map at index {idx}: {e}")
                results.append(None)

        valid_results = [res for res in results if res is not None]
        if valid_results:
            max_result = max(valid_results, key=lambda x: x['max_corr'])
        else:
            raise RuntimeError("No valid results to determine max surface correlation.")
        
        return results, max_result

    def find_map_with_max_corr(self, *maps):
        """
        Finds the map with the highest overall 'corr' value.

        Parameters
        ----------
        *maps : dict
            Variable number of dictionaries, each containing:
            - 'corr': Peak correlation value.
            - 'theta': Theta angle for peak correlation.
            - 'phi': Phi angle for peak correlation.

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
        max_corr_result = {'max_corr': -float('inf'), 'theta': None, 'phi': None, 'map_index': None}
        
        for idx, corr_map in enumerate(maps):
            try:
                corr = corr_map['corr']
                if corr > max_corr_result['max_corr']:
                    max_corr_result.update({'max_corr': corr, 'theta': corr_map['theta'], 
                                            'phi': corr_map['phi'], 'map_index': idx})
            except KeyError as e:
                print(f"Key error in map at index {idx}: {e}")
            except Exception as e:
                print(f"Error processing map at index {idx}: {e}")
        
        if max_corr_result['max_corr'] == -float('inf'):
            raise RuntimeError("No valid correlation values found.")
        
        return max_corr_result


    def calculate_surface_corr_ratio(self, *maps, z_thresh=-10):
        """
        Calculates surface correlation ratio by comparing max surface correlation with max overall correlation.

        Parameters
        ----------
        *maps : dict
            Variable number of dictionaries, each containing:
            - 'map': ROOT.TH2D histogram with correlation values.
            - 'corr': Overall peak correlation value.
            - 'radius': Correlation radius for the station (assumed to be in meters).
        z_thresh : float, optional
            Depth threshold to define theta range lower bound. Default is -10 m.

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
        _, max_surf_corr_result = self.get_surface_corr_max_multiple(*maps, z_thresh=z_thresh)
        max_corr_result = self.find_map_with_max_corr(*maps)
        
        surf_corr_ratio = max_surf_corr_result['max_corr'] / max_corr_result['max_corr'] if max_corr_result['max_corr'] != 0 else float('inf')
        
        return surf_corr_ratio, max_surf_corr_result, max_corr_result

    def min_frac_corr_depth(self, corr_map, fraction=0.6, z_thresh=0):
        """
        Finds the shallowest depth where correlation meets a specified fraction of the maximum correlation,
        using 'theta' to calculate depth and adjusting by the average z-coordinate of antennas.
        Ensures the result remains below the specified depth threshold (z_thresh).

        Parameters
        ----------
        corr_map : dict
            Contains:
            - 'map': ROOT.TH2D histogram with correlation values by theta and phi.
            - 'radius': Correlation radius in meters.
        fraction : float
            Fraction of the max correlation to set the threshold.
            Default is 0.6 (or 60%).
        z_thresh : float
            Maximum allowable depth (relative to the surface). Default is -10 meters.

        Returns
        -------
        min_depth : float
            Shallowest depth (measured from ice surface) where correlation meets the specified fraction of max and below z_thresh.
        """
        # Get maximum correlation and threshold
        max_corr = float(corr_map.get("corr"))
        threshold_corr = fraction * max_corr

        # Access correlation map and radius
        hist = corr_map.get("map")
        radius = corr_map.get("radius")
        if hist is None or radius is None:
            raise ValueError("Map or radius not found in corr_map.")
        radius = float(radius)

        # Calculate average antenna z-coordinate
        _, _, avg_z = mu.calculate_avg_antenna_xyz(self.station_id, self.num_channels)

        min_depth = -float('inf')  # Initialize to find the shallowest depth (least negative)

        # Check correlation values against threshold and calculate depth
        for x_bin in range(1, hist.GetNbinsX() + 1):
            for y_bin in range(1, hist.GetNbinsY() + 1):
                corr_value = hist.GetBinContent(x_bin, y_bin)
                if corr_value >= threshold_corr:
                    # Calculate depth from theta and adjust by avg_z
                    theta = hist.GetYaxis().GetBinCenter(y_bin)
                    depth = radius * math.sin(math.radians(theta)) + avg_z
                    

                    # Ensure min_depth does not exceed z_thresh
                    if depth > min_depth and depth <= z_thresh:
                        min_depth = depth

        if min_depth == -float('inf'):
            raise RuntimeError("No depth meets the fractional correlation threshold within the z_thresh.")
        
        return min_depth

    def min_frac_corr_depth_multiple(self, *maps, fraction=0.6, z_thresh=0):
        """
        Finds the shallowest depth across multiple correlation maps where correlation meets
        a specified fraction of the maximum correlation. Uses the min_frac_corr_depth function for each map,
        then returns the minimum depth found and the index of the corresponding map.

        Parameters
        ----------
        *maps : dict
            A variable number of correlation maps, each containing:
            - 'map': ROOT.TH2D histogram with correlation values by theta and phi.
            - 'radius': Correlation radius in meters.
        fraction : float
            Fraction of the max correlation to set the threshold. Default is 0.6 (60%).
        z_thresh : float
            Maximum allowable depth (relative to the surface). Default is 0 meters.

        Returns
        -------
        min_depth : float
            Shallowest depth (measured from ice surface) across all maps where correlation meets the specified fraction of max.
        min_depth_index : int
            Index of the map that contains the shallowest depth.
        """
        min_depth = -float('inf')  # Initialize to find the shallowest depth (least negative)
        min_depth_index = -1  # Track the index of the map with the minimum depth

        # Iterate over each map and find the shallowest depth using min_frac_corr_depth
        for idx, corr_map in enumerate(maps):
            try:
                # Use the existing function to find the minimum depth for this map
                depth = self.min_frac_corr_depth(corr_map, fraction=fraction, z_thresh=z_thresh)

                # Update min_depth if the current depth is shallower (closer to the surface)
                if depth > min_depth:
                    min_depth = depth
                    min_depth_index = idx

            except RuntimeError as e:
                print(f"Map {idx} did not meet the threshold: {e}")

        if min_depth == -float('inf'):
            raise RuntimeError("No maps meet the fractional correlation threshold within the z_thresh.")

        return min_depth, min_depth_index
