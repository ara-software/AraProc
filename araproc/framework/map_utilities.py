import ctypes
import os
import numpy as np
import ROOT
from functools import partial
from araproc.framework import constants

def get_corr_map_peak(the_map  = None):
        
    """
    Little utility helper function

    Parameters
    ----------
    the_map : ROOT.TH2D
        A 2D map you want the peak found
        
        
    Returns
    peak_corr : float
        The peak in the correlation  map
    peak_phi : int
        The theta value of the peak of the correlation map
    peak_phi : int
        The phi value of the peak of the correlation map
    """
        
    if the_map is None:
       raise Exception("No TH2D map was passed")
    
    _peakZ = ctypes.c_int()
    _peakTheta = ctypes.c_int()
    _peakPhi = ctypes.c_int()
    the_map.GetMaximumBin(_peakPhi, _peakTheta, _peakZ)

    peak_corr = the_map.GetMaximum()
    peak_phi = the_map.GetXaxis().GetBinCenter(_peakPhi.value)
    peak_theta = the_map.GetYaxis().GetBinCenter(_peakTheta.value)
    return peak_corr, peak_phi, peak_theta


def initialize_ARA_geometry(station_id):
    """
    Initialize the geometry tool and station information, providing access to all functions.

    Parameters:
    - station_id (int): ID of the station.

    Returns:
    - Object with dynamically bound functions for use.
    """
    # Load required ARA libraries
    ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR') + "/lib/libAraEvent.so")

    # Constants
    FEET_TO_METER = 0.3048  # Conversion factor from feet to meter
    #translator vector from the Surveyers Coordinate System to the ARA Global Coordinate System
    SUR_TO_GLOB = np.array([22400, 53900, -36.02])  # These are in feets

    # Initialize the Geometry Tool
    geomTool = ROOT.AraGeomTool.Instance()

    # Initialize the station informations
    st_info = geomTool.getStationInfo(station_id)

    # Define all functions
    def get_station_center(excluded_channels=None):
        """
        Calculates the average x, y, and z coordinates of antennas for a given station.

        Parameters
        ----------
        excluded_channels : list
           list of channels (antennas) to exclude.

        Returns
        av_ant_position:tuple
           station center (x,y,z) for given station 
        """ 
        if excluded_channels is None:
            excluded_channels = []

        antenna_coordinates = [[] for _ in range(3)]  # x, y, z
        total_antennas = constants.num_RF_channels

        for antenna_index in range(total_antennas):
            if antenna_index in excluded_channels:
                continue

            antenna_info = st_info.getAntennaInfo(antenna_index)
            for axis in range(3):
                antenna_coordinates[axis].append(antenna_info.antLocation[axis])
        av_ant_position = tuple(np.average(coords) for coords in antenna_coordinates)
        return av_ant_position

    def get_local_CP(cal_pulser_index):
        """
        Parameters
        ----------
        cal_pulser_index: int 
            The index of calpulser (0-3) 
        
        Returns
        local_calpulser: list 
            The [x,y,z] position of calpulser
        """
        cal_pulser_info = st_info.getCalAntennaInfo(cal_pulser_index)
        local_calpulser = [cal_pulser_info.antLocation[i] for i in range(3)]
        return local_calpulser

    def get_survey_to_global_coords(easting, northing, elevation=None):
        """
        Converts Surveyor's coordiates to ARA Global (Array) coordinates.

        Parameters
        ----------
        easting:
          The easting coordiates from surveyour's note in feets
        northing:
          The northing coordiates from surveyour's note in feets

        Returns
        -------        
        The ARA Global Coordinatesi (array) in feets
        """               

        if elevation is not None:
            return np.asarray([easting, northing, elevation]) - SUR_TO_GLOB
        return np.asarray([easting, northing]) - SUR_TO_GLOB[:2]

    def get_global_to_station_centric(array_coords):
        """
        Convert global coordinates to station-centric coordinates.

        Parameters:
        - array_coords (list or numpy array): 
          Global coordinates [x, y, z].

        Returns:
        - output_coords: numpy array
          Station-centric coordinates [x', y', z'].
        """
        # Use numpy array and pass as three arguments
        input_coords = np.array(array_coords, dtype=np.float64)
        output_coords = np.zeros(3, dtype=np.float64)
        geomTool.convertArrayToStationCoords(station_id, input_coords, output_coords)
        return output_coords

    def get_southpole_landmarks(landmark_type):
        """
        Currently We have included three known landmarks at southopole
        Parameters
        ----------
        landmark_type : str
           ICL : stands for IceCube Lab
           WT : Wind Turbine
           SPT : South Pole Telescope   
        Returns
        -------
        landmarks : array
            The landmark location in ARA Gloabal (Array) Coordinate system
        """ 
        landmarks = []
        if landmark_type == "ICL":
            for i in range(1): ## This can go range(4)
                corner = geomTool.getICLCorner(i)
                landmarks.append([corner[0], corner[1], corner[2]])
        elif landmark_type == "WT": ## This can go range(1,4)
            for i in range(3, 4):
                corner = geomTool.getWindTurbine(i)
                landmarks.append([corner[0], corner[1], corner[2]])
        elif landmark_type == "SPT":
            corner = geomTool.getSouthPoleTelescope()
            landmarks.append([corner[0], corner[1], corner[2]])
        else:
            raise ValueError(f"Unknown landmark type: {landmark_type}")
        return landmarks

    def get_distant_pulsers(pulser_name,spice_depth = None):
        """
        returns the station centric coordinate of given distant pulser
        Parameters
        ----------
        pulser_name: str
          IC1S : IceCube 1S pulser
          IC22S : IceCube 22S pulser
          Spice : Spicecore pulser 
        Returns
        -------
        st_centric_coords: array
          X,Y,Z location of pulser in station centric
        """
        if pulser_name == "Spice" and spice_depth is None:
           raise ValueError(f"Unknown Spice pulser depth")
        pulser_coords = {
            "IC1S": [45659.6457, 50490.4199],
            "IC22S": [44884.416, 51444.8819],
            "Spice": [42600, 48800],
        }
        pulser_depths = {"IC1S": -1400, "IC22S": -1450.47,"Spice": spice_depth}
        if pulser_name not in pulser_coords:
            raise ValueError(f"Unknown pulser name: {pulser_name}")

        easting, northing = pulser_coords[pulser_name]
        global_coords = get_survey_to_global_coords(easting, northing) * FEET_TO_METER
        global_coords = [float(global_coords[0]), float(global_coords[1]), 0]
        st_centric_coords = get_global_to_station_centric(global_coords)
        return st_centric_coords[0], st_centric_coords[1], pulser_depths[pulser_name]

    def get_cartesian_to_spherical(coord1, coord2): 
        """
        Parameters
        ----------
        coord1: array 
           coordinate of one location
        coord2: array
           coordinate of second location
        Returns
        -------
        r (radius), 90-theta (elevation) and phi (azimuth)              
        """  
        delta = np.array(coord2) - np.array(coord1)
        x, y, z = delta
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.degrees(np.arccos(z / r))
        phi = np.degrees(np.arctan2(y, x))
        return r, 90 - theta, phi

    return type(
        "GeometryFunctions",
        (),
        {
            "get_station_center": partial(get_station_center),
            "get_local_CP": partial(get_local_CP),
            "get_survey_to_global_coords": partial(get_survey_to_global_coords),
            "get_global_to_station_centric": partial(get_global_to_station_centric),
            "get_southpole_landmarks": partial(get_southpole_landmarks),
            "get_distant_pulsers": partial(get_distant_pulsers),
            "get_cartesian_to_spherical": partial(get_cartesian_to_spherical),
        },
    )()


def get_known_landmarks(station_id,list_of_landmarks = None,cal_pulser_index=None,spice_depth = None):
    """
    Parameters
    ----------    
    station_id : int
    cal_pulser_index : int (0-3)   
      which calpulser you want to see in your skymap
    list_of_landmarks: list
      which landmarks you want to see in your skymap  

    Returns
    -------
    collection : dict
       a dictionary that contains radius , elevation and azimuth of different sources
    """

    if list_of_landmarks is None:
       list_of_landmarks = ["IC22S","ICL"]
    if cal_pulser_index is None:
       cal_pulser_index = 2
 
    collect = {}
    geo = initialize_ARA_geometry(station_id)
    cp = geo.get_local_CP(cal_pulser_index)
    station_center = geo.get_station_center()
    x,y,z = station_center
    rCal, tCal, pCal = geo.get_cartesian_to_spherical(station_center, cp)
    collect["CP"] = [rCal, tCal, pCal]
    del rCal,tCal,pCal 

    for pulser in ["IC1S","IC22S","Spice"]:
        if pulser in list_of_landmarks:
           this_pulser = geo.get_distant_pulsers(pulser,spice_depth)
           r, t, p = geo.get_cartesian_to_spherical(station_center, this_pulser)
           collect[pulser] = [r, t, p]
           del r,t,p
    for known_loc in ["ICL", "WT","SPT"]:
        if known_loc in list_of_landmarks:
           this_loc = geo.get_southpole_landmarks(known_loc)
           this_loc = np.array(this_loc)[0]
           st_centric = geo.get_global_to_station_centric(this_loc[:3])
           r, t, p = geo.get_cartesian_to_spherical(station_center, st_centric)
           collect[known_loc] = [r, t, p]
           del r,t,p
    return collect

