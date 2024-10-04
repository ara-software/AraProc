import numpy as np
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
        - a list of channels which should be exclued from the reconstruction
        - the distance of the cal pulser for that station
    These values are likely config dependent, so we expect this function to change.
    At some point, it may make sense to offload some of these into config files.
    But let's grow as we go.

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
    exclued_channels : array(int)
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

        self.rtc_wrapper = interf.RayTraceCorrelatorWrapper(self.station_id)

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
        ROOT.SetOwnership(pulser_map_v, True)

        corr_pulser_v, phi_pulser_v, theta_pulser_v = mu.get_corr_map_peak(pulser_map_v)
        reco_results["pulser_v"] = {"corr" : corr_pulser_v, 
                                    "theta" : theta_pulser_v,
                                    "phi" : phi_pulser_v,
                                    "map" : pulser_map_v
                                    }

        # check the cal pulser in H
        pulser_map_h = self.rtc_wrapper.correlators["nearby"].GetInterferometricMap(
            self.pairs_h,
            corr_functions_h,
            0
        )
        ROOT.SetOwnership(pulser_map_h, True)

        corr_pulser_h, phi_pulser_h, theta_pulser_h = mu.get_corr_map_peak(pulser_map_h)
        reco_results["pulser_h"] = {"corr" : corr_pulser_h, 
                                    "theta" : theta_pulser_h,
                                    "phi" : phi_pulser_h,
                                    "map" : pulser_map_h
                                    }

        # make a 300 m map in V
        distant_map_v = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
            0
        )
        ROOT.SetOwnership(distant_map_v, True)
        corr_distant_v, phi_distant_v, theta_distant_v = mu.get_corr_map_peak(distant_map_v)
        reco_results["distant_v"] = {"corr" : corr_distant_v, 
                                    "theta" : theta_distant_v,
                                    "phi" : phi_distant_v,
                                    "map" : distant_map_v
                                    }

        # make a 300 m map in H
        distant_map_h = self.rtc_wrapper.correlators["distant"].GetInterferometricMap(
            self.pairs_v,
            corr_functions_v,
            0
        )
        ROOT.SetOwnership(distant_map_h, True)
        corr_distant_h, phi_distant_h, theta_distant_h = mu.get_corr_map_peak(distant_map_h)
        reco_results["distant_h"] = {"corr" : corr_distant_h, 
                                    "theta" : theta_distant_h,
                                    "phi" : phi_distant_h,
                                    "map" : distant_map_h
                                    }

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
        return reco_results
