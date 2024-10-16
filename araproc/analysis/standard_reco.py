import ctypes
import numpy as np
import math
import os
import ROOT

from araproc.analysis import interferometry as interf
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
        
        if station_id not in [1, 2, 3, 4, 5]:
            raise KeyError(f"Station {station_id} is not supported")
        self.station_id = station_id
    
        if not isinstance(excluded_channels, np.ndarray):
            raise KeyError(f"Excluded channel list needs to be a 1D numpy array")
        if excluded_channels.ndim != 1:
            raise AttributeError(f"Excluded channels has the wrong number of dimensions -- should be 1D only")

        num_channels_library = {
            1 : 16,
            2 : 16,
            3 : 16,
            4 : 16,
            5 : 16
        }

        excluded_channels_vec = ROOT.std.vector("int")(excluded_channels)

        # each station has a slightly different distance for the cal pulser reco,
        # so look that up
        self.calpulser_r_library = {
            1 : "48.02",
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

    def do_standard_reco(self, waveform_bundle):
        
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
        
        # check the cal pulser in V
        pulser_map_v = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
            0
        )
        corr_pulser_v, phi_pulser_v, theta_pulser_v = mu.get_corr_map_peak(pulser_map_v)
        reco_results["pulser_v"] = {"corr" : corr_pulser_v, 
                                    "theta" : theta_pulser_v,
                                    "phi" : phi_pulser_v,
                                    "map" : pulser_map_v,
                                    "radius" : self.calpulser_r_library[self.station_id]
                                    }

        # # check the cal pulser in H
        # pulser_map_h = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
        #     self.pairs_h,
        #     corr_functions_h,
        #     0
        # )
        # corr_pulser_h, phi_pulser_h, theta_pulser_h = mu.get_corr_map_peak(pulser_map_h)
        # reco_results["pulser_h"] = {"corr" : corr_pulser_h, 
        #                             "theta" : theta_pulser_h,
        #                             "phi" : phi_pulser_h,
        #                             "map" : pulser_map_h,
        #                             "radius" : self.calpulser_r_library[self.station_id]
        #                             }

        # make a 300 m map in V (Direct rays)
        distant_map_v_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
            0
        )

        # make a 300 m map in V (Refracted/Reflected rays)
        distant_map_v_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
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
            "radius": self.distant_events_r_library[1]
        }

        # Store the refracted/reflected rays results
        reco_results["distant_v_ref"] = {
            "corr": corr_distant_v_ref, 
            "theta": theta_distant_v_ref,
            "phi": phi_distant_v_ref,
            "map": distant_map_v_ref,
            "radius": self.distant_events_r_library[1]
        }

        # make a 300 m map in H (Direct rays)
        distant_map_h_dir = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h,
            corr_functions_h,
            0
        )

        # make a 300 m map in H (Refracted/Reflected rays)
        distant_map_h_ref = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_h,
            corr_functions_h,
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
            "radius": self.distant_events_r_library[1]
        }

        # Store the refracted/reflected rays results in a separate dictionary
        reco_results["distant_h_ref"] = {
            "corr": corr_distant_h_ref, 
            "theta": theta_distant_h_ref,
            "phi": phi_distant_h_ref,
            "map": distant_map_v_ref,
            "radius": self.distant_events_r_library[1]
        }

        del waveform_map
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

    def get_surface_correlation_ratio(self, corr_map, z_thresh=-10):
        """
        Calculates the ratio of the maximum surface correlation value within a specific
        theta range to the overall peak correlation value. The theta range is determined
        based on the average antenna z-coordinate and a threshold.

        Parameters
        ----------
        corr_map : dict
            A dictionary containing:
            - 'map': ROOT.TH2D histogram with correlation values mapped by theta and phi.
            - 'corr': The overall peak correlation value.
            - 'radius': The correlation radius for the station. Assumed to be in meters.
        z_thresh : float, optional
            The depth threshold value to define the lower bound of the theta range. z_thresh=-10 means "10 m below surface"
            Ensure that the unit of `z_thresh` is consistent with the unit of `radius` in `corr_map`.

        Returns
        -------
        surf_corr_ratio : float
            The ratio of the maximum surface correlation value within the theta range 
            to the overall peak correlation.
        max_surf_theta : float
            The theta angle corresponding to the maximum surface correlation value.
        max_surf_phi : float
            The phi angle corresponding to the maximum surface correlation value.
        """
        
        # Extract the radius from corr_map
        corr_radius = float(corr_map.get("radius", None))
        if corr_radius is None:
            raise ValueError("Radius not found in corr_map.")

        # Initialize AraGeomTool instance to retrieve antenna information
        geomTool = ROOT.AraGeomTool.Instance()

        # Collect z-coordinates of all antennas to calculate the average z-coordinate
        z = []
        for ant in range(self.num_channels):
            # Retrieve each antenna's location for the given station
            ant_location = geomTool.getStationInfo(self.station_id).getAntennaInfo(ant).antLocation
            z.append(ant_location[2])  # Append only the z-coordinate

        # Calculate the average z-coordinate of the antennas
        avg_z = sum(z) / len(z)

        # Check if the radius is large enough to visualize the surface
        if corr_radius < (abs(avg_z) + z_thresh):
            raise RuntimeError("Surface not visible: radius is smaller than average z coordinate of the antennas.")

        # Calculate the angles at the surface and z threshold relative to the horizontal
        theta_at_surface = math.degrees(math.asin(abs(avg_z) / corr_radius))
        theta_at_z_thresh = math.degrees(math.asin((abs(avg_z) + z_thresh) / corr_radius))

        # Access the ROOT.TH2D histogram from 'map' in corr_map
        map = corr_map.get("map", None)
        if map is None:
            raise ValueError("The 'map' key was not found in corr_map.")

        # Initialize lists to store correlation values within the specified theta range
        surf_corr = []
        tot_corr_idx = []

        # Define the theta range based on theta_at_surface and theta_at_z_thresh
        theta_min = min(theta_at_surface, theta_at_z_thresh)
        theta_max = max(theta_at_surface, theta_at_z_thresh)

        # Iterate over each bin in the map to filter by theta range
        nbins_x = map.GetNbinsX()
        nbins_y = map.GetNbinsY()
        for binx in range(1, nbins_x + 1):
            for biny in range(1, nbins_y + 1):
                # Get theta and phi for the current bin and retrieve its correlation value
                theta = map.GetYaxis().GetBinCenter(biny)
                phi = map.GetXaxis().GetBinCenter(binx)
                corr_value = map.GetBinContent(binx, biny)

                # Store correlation values and bin indices if within the theta range
                if theta_min <= theta <= theta_max:
                    surf_corr.append(corr_value)
                    tot_corr_idx.append((binx, biny))

        # Ensure that there are correlation values within the theta range
        if not surf_corr:
            raise RuntimeError("No correlation values found in the specified theta range!")

        # Calculate the maximum correlation within the specified theta range
        max_surf_corr = np.nanmax(surf_corr)
        max_surf_corr_index = np.nanargmax(surf_corr)

        # Identify bin location of the maximum correlation value
        max_surf_bin = tot_corr_idx[max_surf_corr_index]
        max_surf_bin_x, max_surf_bin_y = max_surf_bin
        max_surf_theta = map.GetYaxis().GetBinCenter(max_surf_bin_y)
        max_surf_phi = map.GetXaxis().GetBinCenter(max_surf_bin_x)

        # Access the overall peak correlation value from 'corr' in corr_map
        corr_value = corr_map.get("corr", None)
        if corr_value is None:
            raise ValueError("The 'corr' key was not found in corr_map.")

        # Calculate the ratio of the max surface correlation to the peak correlation
        surf_corr_ratio = max_surf_corr / corr_value

        return surf_corr_ratio, max_surf_theta, max_surf_phi