import ctypes
import numpy as np
import ROOT
from araproc.framework import constants
from itertools import combinations

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
        #print(" exc chan = ",excluded_channels)
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


    def get_raytraced_theta_phi(self, xyz_rel, excluded_channels, radius_map, reco,
                                     which_distance="distant", solution=0):
        """
        Find the ray-traced elevation (90-theta) and azimuth (phi) where a landmark
        intersects the reconstruction sphere of radius `radius_map`.

        Parameters
        ----------
        xyz_rel : array-like
        [x, y, z] of landmark relative to station center (SC) in meters.
        radius_map : float
        Radius of reconstruction sphere (m)
        reco : StandardReco 
        which_distance : str "nearby" or "distant"
        solution : int
        0 = direct ray, 1 = reflected/refracted ray

        Returns
        -------
        theta_RT : float
        Ray-traced elevation (degrees)
        phi_RT : float
        Azimuth (degrees, from x-axis)
        """
        xyz_np = np.array(xyz_rel, dtype=float)
        r_lm = np.linalg.norm(xyz_np) 

        # Straight-line projection to sphere of radius = radius_map
        xyz_proj = xyz_np * (radius_map / r_lm) 
        x, y, z = xyz_proj

        phi_RT = np.degrees(np.arctan2(y, x))  # azimuth
        theta_SL = np.degrees(np.arcsin(np.clip(z / radius_map, -1.0, 1.0)))  # straight-line elevation, should be between -1 to +1 
        #print("   excluded_channels are ",excluded_channels) # should exclude channels when considering the dt to theta,phi scanning

        # Get all valid VPol and HPol channels after removing excluded channels
        vpol_channels = [ch for ch in constants.vpol_channel_ids if ch not in excluded_channels]
        hpol_channels = [ch for ch in constants.hpol_channel_ids if ch not in excluded_channels]

        # Make all possible pairs within each polarization
        valid_vpol_pairs = list(combinations(vpol_channels, 2))
        valid_hpol_pairs = list(combinations(hpol_channels, 2))
        
        #print("Valid VPol pairs:", valid_vpol_pairs,"  no of VPOL pair = ",len(valid_vpol_pairs))
        #print("Valid HPol pairs:", valid_hpol_pairs,"   no of HPOL pair = ",len(valid_hpol_pairs))
        valid_pairs = valid_vpol_pairs + valid_hpol_pairs

        if not valid_pairs:
            print("WARNING!!! No valid inter-string VPol pairs; falling back to straight line solution for Landmark.")
            return theta_SL, phi_RT

        theta_list = []
        theta_grid = np.arange(-90, 91, 1.0)  # lookup grid 
        #consider all pairs of same polarization for the expected time difference
        for ch_a, ch_b in valid_pairs:
            pos_a = np.array(self.st_info.getAntennaInfo(ch_a).antLocation)
            pos_b = np.array(self.st_info.getAntennaInfo(ch_b).antLocation)
            #n(z) at average depth using the antenna pairs' depths
            n_eff = (constants.get_index_of_refraction(pos_a[2]) +
                    constants.get_index_of_refraction(pos_b[2])) / 2.0
            v_eff = constants.speed_of_light / n_eff # velocity v(z)

            dt_target = (np.linalg.norm(xyz_proj - pos_a) - np.linalg.norm(xyz_proj - pos_b)) / v_eff # expected time difference

            # You can lookup Brian's ray trace table over theta grid for reference
            # phi_RT  ~constant, n(z) varies with depth z

            dt_table = np.array([
            reco.lookup_arrival_time(ch_a, t, phi_RT, which_distance, solution) -
            reco.lookup_arrival_time(ch_b, t, phi_RT, which_distance, solution)
            for t in theta_grid
            ]) # ray-traced arrival time for a (given r=radius_map, theta, phi) per  channel

            residual = dt_table - dt_target
            #print('residual   = ',residual)
            sign_changes = np.where(np.diff(np.sign(residual)))[0] 

            if len(sign_changes) == 0:
                theta_list.append(theta_SL) # fall back if no sign change = straight line solution as done before
                continue

            # pick crossing closest to straight-line theta
            idx = sign_changes[np.argmin(np.abs(theta_grid[sign_changes] - theta_SL))]
            t0, t1 = theta_grid[idx], theta_grid[idx + 1]
            r0, r1 = residual[idx], residual[idx + 1]
            theta_pair = t0 - r0 * (t1 - t0) / (r1 - r0)  # assuming theta_pair linearly between t0 and t1
            theta_list.append(theta_pair)

        # Final ray-traced theta = median of all inter-string pairs
        theta_RT = float(np.median(theta_list))
        return theta_RT, phi_RT



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


    def get_known_landmarks(self, list_of_landmarks=None, list_of_excluded_chans=None, list_of_cal_pulser_indices=None,
                            spice_depth=None, reco=None, solution=0): 
        """
        Parameters
        ----------
        list_of_cal_pulser_indices : list ## example [0,1,2,3]
          which calpulsers you want to see in your skymap
        list_of_landmarks: list ## example ['ICL','IC22S','SPT','IC1S','SPIce','WT']
          which landmarks you want to see in your skymap
        spice_depth : int/float
          the depth of spice pulser ## example -1451.3 
        Note: Added ray-traced solution for distant pulsers [IC22S, IC1S, SPICE] and southpole landmarks [ICL, WT, SPRESSO and SPT]. Use ray-tracing to find
              where the ray coming from these sources interset the reconstruction sphere (with Radius R = 300, refer to ray-tracing tables generated by Brian for direct
              and ref. ray-traced arrival time per chan (VPOL and HPol) for theta (-90 to +90) deg, and phi (-180, 180). At the intersection point, theta, phi are 
              estimated and plotted on the skymap. 
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
        collect = {}
        
        calpulser = self.get_local_CP(list_of_cal_pulser_indices)
        station_center = self.get_station_center()
        for cp in calpulser:
            rCal, tCal, pCal = self.get_relative_cartesian_to_spherical(station_center, calpulser[cp])
            collect[f"CP{cp}"] = [rCal, tCal, pCal]
            del rCal, tCal, pCal
        del calpulser
        # Adding ray-traced landmarks for southpole landmarks and distant pulser
        # first: distant pulser 
        
        R_distant_map = float(reco.rtc_wrapper.correlators["distant"].GetRadius()) if reco is not None else None #The radius of distant map
        for pulser in ["IC1S", "IC22S", "SPIce"]:
            if pulser in list_of_landmarks:
                this_pulser = self.get_distant_pulsers(pulser, spice_depth)
                r, t, p = self.get_relative_cartesian_to_spherical(station_center, this_pulser)
                if reco is not None:
                    xyz_rel = np.array(this_pulser) - np.array(station_center)
                    t, p = self.get_raytraced_theta_phi(
                            xyz_rel, list_of_excluded_chans, R_distant_map, reco, 
                            which_distance="distant",
                            solution=solution)
                collect[pulser] = [r, t, p] # ray traced solution for distant pulser landmark
                del r,t,p 
        #next : known south pole landmarks (ray-traced solution)

        for known_loc in ["ICL", "WT", "SPRESSO", "SPT"]:
            if known_loc in list_of_landmarks:
                this_loc = self.get_southpole_landmarks(known_loc)
                this_loc = np.array(this_loc)[0]
                st_centric = self.get_global_to_station_centric(this_loc[:3])
                r, t, p = self.get_relative_cartesian_to_spherical(station_center, st_centric)
                if reco is not None:
                    xyz_rel = np.array(st_centric) - np.array(station_center)
                    t, p = self.get_raytraced_theta_phi(
                            xyz_rel, list_of_excluded_chans, R_distant_map, reco,
                            which_distance="distant",
                            solution=solution)
                collect[known_loc] = [r, t, p]
                del r,t,p
        collect['critical_angle'] = self.get_critical_angle()
        return collect

