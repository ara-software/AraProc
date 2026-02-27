import ctypes
import numpy as np
import ROOT
from araproc.framework import constants


def get_corr_map_peak(the_map  = None):
    """
    Little utility helper function

    Parameters
    ----------
    the_map : ROOT.TH2D
        A 2D map you want the peak found
        
    Returns
    -------
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



def load_raytrace_model(
    arasimsrc="/cvmfs/icecube.opensciencegrid.org/users/abishop/AraSim-main"
):
    """
    Load AraSim RayTrace into ROOT interpreter.
    
    """

    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace_IceModels.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/Vector.h"')

    ROOT.gInterpreter.ProcessLine(
        'auto _rt_atten = boost::shared_ptr<basicAttenuationModel>'
        '(new basicAttenuationModel);'
    )

    ROOT.gInterpreter.ProcessLine(
        f'auto _rt_refr = boost::shared_ptr<exponentialRefractiveIndex>'
        f'(new exponentialRefractiveIndex({1.3260},{1.780},{0.0202}));'   # no, nf and l from PA paper
    )
    
    ROOT.gInterpreter.ProcessLine(
        'RayTrace::TraceFinder _rt_tf(_rt_refr, _rt_atten);'
    )

    ROOT.gInterpreter.ProcessLine('Vector _rt_src; Vector _rt_rec;')



load_raytrace_model()


class AraGeom:
    """
    A class to initialize ARA geometry and provide associated utility functions.
    """


    def __init__(self, station_id):
        """
        Initialize the geometry tool and station information.
        """
        # Constants
        self.FEET_TO_METER = 0.3048  # Conversion factor from feet to meter
        self.SUR_TO_GLOB = np.array([22400, 53900, -36.02])  # Survey to global offset in feet

        self.station_id = station_id
        # Initialize the Geometry Tool
        self.geomTool = ROOT.AraGeomTool.Instance()

        # Station information
        self.st_info = self.geomTool.getStationInfo(self.station_id)


    def get_station_center(self, excluded_channels=None):
        """
        Calculates the average x, y, and z coordinates of antennas for a given station.

        Parameters
        ----------
        excluded_channels : list
           list of channels (antennas) to exclude.

        Returns
        -------
        av_ant_position:tuple
           station center (x, y, z) for given station
        """

        if excluded_channels is None:
            excluded_channels = []

        antenna_coordinates = [[] for _ in range(3)]  # x, y, z
        total_antennas = constants.num_rf_channels

        for antenna_index in range(total_antennas):
            if antenna_index in excluded_channels:
                continue

            antenna_info = self.st_info.getAntennaInfo(antenna_index)
            for axis in range(3):
                antenna_coordinates[axis].append(antenna_info.antLocation[axis])
        av_ant_position = tuple(np.average(coords) for coords in antenna_coordinates)
        return av_ant_position


    def get_local_CP(self, list_of_cal_pulser_indices):
        """
        Parameters
        ----------
        list_of_cal_pulser_indices: list
            The indices of calpulser [0,1,2,3] Note: 0 & 2 are HPols and 1 & 3 are VPols (0&1 are in one string and 2&3 are in another string)

        Returns
        -------
        calpulsers: list
            A dictionary of [x, y, z] position of calpulsers
        """

        calpulsers = {}
        for cal_index in list_of_cal_pulser_indices:
            cal_pulser_info = self.st_info.getCalAntennaInfo(int(cal_index))
            calpulsers[f'{cal_index}'] = [cal_pulser_info.antLocation[i] for i in range(3)]
        return calpulsers


    def get_survey_to_global_coords(self, easting, northing, elevation=None):
        """
        Converts Surveyor's coordiates to ARA Global (Array) coordinates.

        Parameters
        ----------
        easting:
          The easting coordiate from surveyour's note in feets
        northing:
          The northing coordiate from surveyour's note in feets
        elevation:
          The elevation coordiate from surveyour's note in feets

        Returns
        -------
        The ARA Global Coordinates (array) in feets
        """

        if elevation is not None:
            return np.asarray([easting, northing, elevation]) - self.SUR_TO_GLOB
        return np.asarray([easting, northing]) - self.SUR_TO_GLOB[:2]
    

    def get_array_from_lat_long(self,latitude,longitude):
        easting = self.geomTool.getEastingFromLatLong(latitude,longitude)
        northing = self.geomTool.getNorthingFromLatLong(latitude,longitude)
        return np.array([northing, easting,0])

    
    def get_global_to_station_centric(self, array_coords):
        """
        Convert global coordinates to station-centric coordinates.

        Parameters
        ----------
        - array_coords (list or numpy array):
          Global coordinates [x, y, z].

        Returns
        -------
        - output_coords: numpy array
          Station-centric coordinates [x', y', z'].
        """

        input_coords = np.array(array_coords, dtype=np.float64)
        output_coords = np.zeros(3, dtype=np.float64)
        self.geomTool.convertArrayToStationCoords(self.station_id, input_coords, output_coords)
        return output_coords


    def get_southpole_landmarks(self, landmark_type):
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
            for i in range(1):  ## This can go range(4)
                corner = self.geomTool.getICLCorner(i)
                landmarks.append([corner[0], corner[1], corner[2]])
        elif landmark_type == "WT":  ## This can go range(1,4)
            for i in range(3, 4):
                corner = self.geomTool.getWindTurbine(i)
                landmarks.append([corner[0], corner[1], corner[2]])
        elif landmark_type == "SPRESSO":
            landmarks.append(self.get_array_from_lat_long(-89.93120, 144.51249))
        elif landmark_type == "SPT":
            corner = self.geomTool.getSouthPoleTelescope()
            landmarks.append([corner[0], corner[1], corner[2]])
        else:
            raise ValueError(f"Unknown landmark type: {landmark_type}")
        return landmarks


    def get_distant_pulsers(self, pulser_name, spice_depth=None):
        """
        returns the station centric coordinate of given distant pulser
        Parameters
        ----------
        pulser_name: str
          IC1S : IceCube 1S pulser
          IC22S : IceCube 22S pulser
          SPIce : South Pole IceCore pulser
        northing:
          The northing coordiates from surveyour's note in feets

        Returns
        -------
        st_centric_coords: array
          X, Y, Z location of pulser in station centric
        """

        if pulser_name == "SPIce" and spice_depth is None:
            raise ValueError(f"Unknown SPIce pulser depth")
        pulser_coords = {
            "IC1S": [45659.6457, 50490.4199],
            "IC22S": [44884.416, 51444.8819],
            "SPIce": [42600, 48800],
        }
        pulser_depths = {"IC1S": -1400, "IC22S": -1450.47, "SPIce": spice_depth}
        if pulser_name not in pulser_coords:
            raise ValueError(f"Unknown pulser name: {pulser_name}")

        easting, northing = pulser_coords[pulser_name]
        global_coords = self.get_survey_to_global_coords(easting, northing) * self.FEET_TO_METER
        global_coords = [float(global_coords[0]), float(global_coords[1]), 0]
        st_centric_coords = self.get_global_to_station_centric(global_coords)
        return st_centric_coords[0], st_centric_coords[1], pulser_depths[pulser_name]


    def get_relative_cartesian_to_spherical(self, reference_point, target_point):
        """
        Compute the spherical coordinates of a vector defined by two Cartesian points.

        This function calculates the vector connecting `reference_point` to `target_point`,
        transforms the coordinate system to center on `reference_point`, and converts the
        resulting vector into spherical coordinates.

        Parameters
        ----------
        reference_point : array-like
            The reference point in Cartesian coordinates [x, y, z].
        target_point : array-like
            The target point in Cartesian coordinates [x, y, z].

        Returns
        -------
        r : float
            The magnitude of the vector (distance between the points).
        elevation : float
            The elevation angle (90 - theta) in degrees.
        azimuth : float
            The azimuth angle (phi) in degrees, measured from the x-axis in the xy-plane.
        """

        delta = np.array(target_point) - np.array(reference_point)
        x, y, z = delta
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.degrees(np.arccos(z / r))
        phi = np.degrees(np.arctan2(y, x))
        return r, 90 - theta, phi



    def get_critical_angle(self):
        """
        Returns
        -------
        critical angle in terms of elevation angle

        """
        station_depth = self.get_station_center()[2]
        n_air = constants.get_index_of_refraction(1e-3)
        n_ice = constants.get_index_of_refraction(station_depth)
        critical_angle = np.arcsin(n_air/n_ice)
        critical_angle *= (180/np.pi)

        return 90 - critical_angle  ## In terms of elevation angle



    #get ray-traced angle at the center of the station

    def get_raytraced_arrival_angles(self, source_xyz, R_sphere, accuracy=0.001, solution=None):

        src = np.array(source_xyz, dtype=float)
        station_center = np.array(self.get_station_center())

        
        z_rec = float(station_center[2])

        dx = src[0] - station_center[0]
        dy = src[1] - station_center[1]

        r_src = float(np.sqrt(dx**2 + dy**2))
        az_phi = float(np.arctan2(dy, dx))

        ROOT._rt_src.SetXYZ(r_src, 0.0, float(src[2]))
        ROOT._rt_rec.SetXYZ(0.0, 0.0, z_rec)

        sol_cnt = ctypes.c_int()
        sol_err = ctypes.c_int()

        allowed = ROOT.RayTrace.NoReflection if solution == 0 else ROOT.RayTrace.AllReflections

        paths = ROOT._rt_tf.findPaths(
            ROOT._rt_src,
            ROOT._rt_rec,
            0.5,
            ROOT.TMath.Pi()/2,
            sol_cnt,
            sol_err,
            allowed,
            accuracy
        )

        if not paths or len(paths) == 0:
            return []

        solutions = []

        for i, path in enumerate(paths):

            
            zenith_deg = np.degrees(path.receiptAngle)
            elevation  = zenith_deg - 90.0 
            azimuth    = np.degrees(az_phi)

            solutions.append({
                'elevation': elevation,
                'azimuth': azimuth,
                'tof': float(path.pathTime) * 1e9,
                'path_len': float(path.pathLen),
                'reflectionAngle': float(path.reflectionAngle),
            })

        return solutions



    def get_known_landmarks(self, list_of_landmarks=None, excluded_channels=None, list_of_cal_pulser_indices=None, spice_depth=None, solution=None):

        """
        Parameters
        ----------
        list_of_cal_pulser_indices : list ## example [0,1,2,3]
          which calpulsers you want to see in your skymap
        list_of_landmarks: list ## example ['ICL','IC22S','SPT','IC1S','SPIce','WT']
          which landmarks you want to see in your skymap
        spice_depth : int/float
          the depth of spice pulser ## example -1451.3 

        Returns
        -------
        collect : dict
           a dictionary that contains radius, elevation, and azimuth of different sources
        """

        if list_of_landmarks is None:
            list_of_landmarks = ["IC22S", "ICL", "SPRESSO"]
        elif list_of_landmarks == ['all']:
             list_of_landmarks = ['ICL','IC22S','SPT','IC1S','WT','SPRESSO']   
        if list_of_cal_pulser_indices is None:
            list_of_cal_pulser_indices = [1,3]
        elif list_of_cal_pulser_indices == ['all']:
             list_of_cal_pulser_indices = [0,1,2,3]
        if spice_depth is not None:
           list_of_landmarks.append('SPIce')

        station_center = self.get_station_center()
        collect = {}

    
        def _resolve(name, source_xyz, R_sphere, solution=None):
            """
            Compute straight-line AND ray-traced angles depending solution choice 0 or 1

            return [r, elevation, azimuth] using ray-traced if available] 
            If solution is not availble, fall back to straight line solution

            """
            r_sl, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(station_center, source_xyz)
            rt_sols = self.get_raytraced_arrival_angles(
            source_xyz, R_sphere, solution=solution)

            if not rt_sols:
                print(f"  [{name}]  RT  :  NO SOLUTION — fallback to SL")
                return [r_sl, elev_sl, az_sl]

            # Label each path D or R using AraSim's noReflection sentinel (100.0 rad)
            NO_REFL = 100.0
            for s in rt_sols:
                s['label'] = 'D' if abs(float(s.get('reflectionAngle', NO_REFL)) - NO_REFL) < 0.01 else 'R'

            # Prefer R if solution=1, fallback to D; if solution=0 or None, take shortest patth
            if solution == 1:
                preferred = sorted([s for s in rt_sols if s['label'] == 'R'], key=lambda x: x['path_len'])
                fallback  = sorted([s for s in rt_sols if s['label'] == 'D'], key=lambda x: x['path_len'])
            elif solution == 0:
                preferred = sorted(rt_sols, key=lambda x: x['path_len'])
                fallback  = []
            else:
                preferred = sorted(rt_sols, key=lambda x: x['path_len'])
                fallback  = []

            chosen = preferred[0] if preferred else (fallback[0] if fallback else None)

            if chosen is None:
                print(f"  [{name}]  no suitable solution — fallback to SL")
                return [r_sl, elev_sl, az_sl]

            return [r_sl, chosen['elevation'], chosen['azimuth']]



        #Adding landmark for calpulsers

        R_nearby   = float(constants.calpulser_r_library[self.station_id]) 

        calpulser = self.get_local_CP(list_of_cal_pulser_indices)
        
        for cp in calpulser:
            collect[f"CP{cp}"] = _resolve(f"CP{cp}", calpulser[cp], R_nearby, solution=0) #CPs hardcoded to direct solution,
            
        del calpulser


        #Adding landmark for Distant Landmarrks

        R_distant = float(constants.distant_events_r_library[self.station_id])
          
        for pulser in ["IC1S", "IC22S", "SPIce"]:
            if pulser in list_of_landmarks:
                this_pulser = self.get_distant_pulsers(pulser, spice_depth)
                collect[pulser] = _resolve(pulser, this_pulser, R_distant, solution)
                

        for known_loc in ["ICL", "WT", "SPRESSO", "SPT"]:
            if known_loc in list_of_landmarks:
                this_loc = self.get_southpole_landmarks(known_loc)
                this_loc = np.array(this_loc)[0]
                st_centric = self.get_global_to_station_centric(this_loc[:3])
                collect[known_loc] = _resolve(known_loc, st_centric, R_distant, solution)
                
        collect['critical_angle'] = self.get_critical_angle()
        
        return collect

