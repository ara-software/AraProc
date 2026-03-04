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
    arasimsrc="/cvmfs/icecube.opensciencegrid.org/users/abishop/AraSim-main"):
    
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace_IceModels.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/Vector.h"')
    ROOT.gInterpreter.ProcessLine('#include <vector>')

    
    ROOT.gInterpreter.ProcessLine(
        'auto _rt_atten = boost::shared_ptr<basicAttenuationModel>'
        '(new basicAttenuationModel);'
    )
    ROOT.gInterpreter.ProcessLine(
        f'auto _rt_refr = boost::shared_ptr<exponentialRefractiveIndex>'
        f'(new exponentialRefractiveIndex({1.3260},{1.780},{0.0202}));' # same as PA analysis paper
    )
    ROOT.gInterpreter.ProcessLine('RayTrace::TraceFinder _rt_tf(_rt_refr, _rt_atten);')
    ROOT.gInterpreter.ProcessLine('Vector _rt_src; Vector _rt_rec;')

    ROOT.gInterpreter.Declare(r'''
        #include <vector>

        static std::vector<double> _rt_x;
        static std::vector<double> _rt_z;
        static std::vector<int>    _rt_stepType;

        struct RTPointSaver {
            void operator()(const RayTrace::fullRayPosition& p, RayTrace::RKStepType t) const {
                _rt_x.push_back(p.x);
                _rt_z.push_back(p.z);
                _rt_stepType.push_back((int)t);
            }
        };

        RayTrace::TraceRecord rt_doTrace_savePoints(
            double src_depth,
            double launch_theta,
            const RayTrace::rayTargetRecord& target,
            unsigned short allowedReflections,
            double frequency,
            double polarization,
            int &sol_error
        ){
            _rt_x.clear();
            _rt_z.clear();
            _rt_stepType.clear();

            RTPointSaver cb;

            return _rt_tf.doTrace<RayTrace::fullRayPosition, RTPointSaver>(
                src_depth, launch_theta, target, allowedReflections,
                frequency, polarization, sol_error, cb
            );
        }

        const std::vector<double>& rt_get_x(){ return _rt_x; }
        const std::vector<double>& rt_get_z(){ return _rt_z; }
        const std::vector<int>&    rt_get_stepType(){ return _rt_stepType; }
    ''')

    

    if not hasattr(ROOT, "rt_doTrace_savePoints"):
        raise RuntimeError(
            "Failed to declare rt_doTrace_savePoints in ROOT. ")

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



    def get_raytraced_sphere_crossing_angles(
        self,
        source_xyz,
        R_sphere,
        accuracy=0.001,
        solution=None,          
        frequency=0.5,
        polarization=0.0):
    
        """
        Ray-traced arrival direction using the sphere-crossing method, from the returned path from ray traced solution, find that particular (x,y,z) that intersects with radius map 

        direction is StationCenter(receiver) to crossing point on a sphere of radius R_sphere.

        
        - solution==0: try NoReflection first (if map type is everything except distant_*_ref, solution = 0 = direcr Ray Tracing). If no path, try AllReflections.
                        choose D if present else R.

        - solution==1: use AllReflections. choose R if present else choose D to have RT solution always. 

        -iff you do not get D or R from ray tracing, print a message and fallback to SL (never had a failed case like this except for landmark = SPT, all other landmark returns
             Ray-Traced solutions for all stations) 

        Returns:
          list of solutions (ordered by preference, shortest path_len first),
          each dict has keys:
            label ('D'/'R'), elevation, azimuth, tof, path_len, reflectionAngle, allowedUsed
        """

        src = np.array(source_xyz, dtype=float)
        rec = np.array(self.get_station_center(), dtype=float)

        z_src = float(src[2])
        z_rec = float(rec[2])

        dx = float(src[0] - rec[0])
        dy = float(src[1] - rec[1])

        r_src  = float(np.sqrt(dx*dx + dy*dy))     
        az_phi = float(np.arctan2(dy, dx))         


        ROOT._rt_src.SetXYZ(r_src, 0.0, z_src)
        ROOT._rt_rec.SetXYZ(0.0, 0.0, z_rec)

        def _find_paths_with_allowed(allowed_reflections):
            sol_cnt = ctypes.c_int()
            sol_err = ctypes.c_int()
            paths = ROOT._rt_tf.findPaths(
                ROOT._rt_src,
                ROOT._rt_rec,
                float(frequency),
                ROOT.TMath.Pi()/2,
                sol_cnt,
                sol_err,
                int(allowed_reflections),
                float(accuracy)
            )
            return paths, sol_cnt.value, sol_err.value



        #   Allowed reflections selection 
        if solution == 0: 
            #Prefer direct-only first
            paths, nsol, serr = _find_paths_with_allowed(ROOT.RayTrace.NoReflection)
            allowed_used = int(ROOT.RayTrace.NoReflection)

            
            if (paths is None) or (nsol == 0) or (len(paths) == 0):
            # Fallback to allow Allreflections for the Direct solution 
                paths, nsol, serr = _find_paths_with_allowed(ROOT.RayTrace.AllReflections)
                allowed_used = int(ROOT.RayTrace.AllReflections)
        else: 
            paths, nsol, serr = _find_paths_with_allowed(ROOT.RayTrace.AllReflections)
            allowed_used = int(ROOT.RayTrace.AllReflections)

        if (nsol == 0) or (not paths) or (len(paths) == 0):
            return []


        NO_REFL_SENTINEL = 100.0  # AraSim convention when reflections are allowed
        labeled = []
        for p in paths:
            if allowed_used == int(ROOT.RayTrace.NoReflection):
                label = "D"
            else:
                refl = float(p.reflectionAngle)
                label = "D" if abs(refl - NO_REFL_SENTINEL) < 0.01 else "R"
            labeled.append((p, label))

        
        def _best_of(lbl):
            c = [(p, lab) for (p, lab) in labeled if lab == lbl]
            if not c:
                return None
            c.sort(key=lambda t: float(t[0].pathLen))
            return c[0]

        cand_D = _best_of("D")
        cand_R = _best_of("R")

        to_trace = []
        if cand_D is not None:
            to_trace.append(cand_D)
        if cand_R is not None:
            to_trace.append(cand_R)

        R = float(R_sphere)
        out = []

        # Exact intersection of segment with circle of radius = ray traced map created by Brian 

        def _segment_circle_intersect(p0, p1, R):
            x0, z0 = p0
            x1, z1 = p1
            dxs = x1 - x0
            dzs = z1 - z0

            a = dxs*dxs + dzs*dzs
            b = 2.0*(x0*dxs + z0*dzs)
            c = (x0*x0 + z0*z0) - R*R

            disc = b*b - 4.0*a*c
            if disc < 0.0 or a == 0.0:
                return None

            sq = np.sqrt(disc)
            t1 = (-b - sq) / (2.0*a)
            t2 = (-b + sq) / (2.0*a)

            for t in (t1, t2):
                if 0.0 <= t <= 1.0:
                    return (x0 + t*dxs, z0 + t*dzs, t)
            return None

        for (path, label) in to_trace:
            target = ROOT.RayTrace.rayTargetRecord(z_rec, r_src)
            sol_error = ctypes.c_int(0)

            _ = ROOT.rt_doTrace_savePoints(
                z_src,
                float(path.launchAngle),
                target,
                int(allowed_used),
                float(frequency),
                float(polarization),
                sol_error
            )

            xs = np.array(list(ROOT.rt_get_x()), dtype=float)
            zs = np.array(list(ROOT.rt_get_z()), dtype=float)

            if sol_error.value != 0 or len(xs) < 2: # doTrace failed:
                continue

            
            xp = r_src - xs
            dz = zs - z_rec
            d_to_rec = np.sqrt(xp*xp + dz*dz)

            # Find crossing closest to receiver: outside to inside when scanning from receiver side
            i_cross = None
            for i in range(len(d_to_rec) - 1, 0, -1):
                if (d_to_rec[i] <= R) and (d_to_rec[i-1] > R):
                    i_cross = i
                    break

            if i_cross is None:
                
                continue

            p_out = (float(xp[i_cross-1]), float(dz[i_cross-1]))  # outside
            p_in  = (float(xp[i_cross]),   float(dz[i_cross]))    # inside

            inter = _segment_circle_intersect(p_out, p_in, R) # crossing point 
            if inter is None:
                
                continue

            xp_c, dz_c, _t = inter
            dcheck = np.sqrt(xp_c*xp_c + dz_c*dz_c)

            vx = xp_c * np.cos(az_phi)
            vy = xp_c * np.sin(az_phi)
            vz = dz_c

            horiz = np.sqrt(vx*vx + vy*vy)
            elev = np.degrees(np.arctan2(vz, horiz))
            az   = np.degrees(np.arctan2(vy, vx))

            

            out.append({
                "label": label,
                "elevation": float(elev),
                "azimuth": float(az),
                "tof": float(path.pathTime) * 1e9,
                "path_len": float(path.pathLen),
                "reflectionAngle": float(path.reflectionAngle),
                "allowedUsed": int(allowed_used),
            })

        if not out:
            return []

        # Apply preference for the ray traced solution based on solution flag (if not distant_*_ref, solution for all other map_type is = Direct = solution 0)
        if solution == 1: 
            preferred = [s for s in out if s["label"] == "R"]
            fallback  = [s for s in out if s["label"] == "D"]
        else:
            preferred = [s for s in out if s["label"] == "D"]
            fallback  = [s for s in out if s["label"] == "R"]

        preferred.sort(key=lambda s: s["path_len"])
        fallback.sort(key=lambda s: s["path_len"])

        return preferred if preferred else fallback









    def get_known_landmarks(self, list_of_landmarks=None, list_of_cal_pulser_indices=None, spice_depth=None):

        """
        IMPORTANT: 
          - We ALWAYS expect an RT solution. If not, print a message and fall back to SL
            so the script continues, but you will immediately see the failure per landmkark. 

        solution:
          
          - solution==0: prefer D (NoReflection first, fallback AllReflections is soluution 0 could not find D)
          - solution==1: prefer R (AllReflections, fallback D, iff R solution could not find R from ray tracing)
        
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
        collect = {}

        station_center = self.get_station_center(excluded_channels=excluded_channels)
        
        def _resolve(name, source_xyz, R_sphere, solution=None):

            r_sl, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(station_center, source_xyz)

            rt_sols = self.get_raytraced_sphere_crossing_angles(
                source_xyz,
                R_sphere,
                accuracy=0.001,
                solution=solution,
                )

            if not rt_sols: # only if Ray Traced results can not find D or R solution, it chooses Straight line solution for the skymap as before
                # may want to keep the prrint statements below if you want to know any landmark returned Straight Line solution
                #print(f"[CRITICAL] NO RAY-TRACED SOLUTION for landmark [{name}] "
                #      f"(solution={solution}, R_sphere={R_sphere}). THIS SHOULD NOT HAPPEN for any station_id EXCEPT for landmark = SPT ")
                #print("Falling back to Straight-Line ONLY so the script can continue.")
            
                
                return [r_sl, elev_sl, az_sl]

            chosen = rt_sols[0]

            return [r_sl, chosen["elevation"], chosen["azimuth"]]

        R_nearby   = float(constants.calpulser_r_library[self.station_id])
        calpulser = self.get_local_CP(list_of_cal_pulser_indices)

        
        for cp in calpulser:
            collect[f"CP{cp}"] = _resolve(f"CP{cp}", calpulser[cp], R_nearby, solution=solution)
            
        del calpulser
        
        # Distant landmarks

        for pulser in ["IC22S", "IC1S", "SPIce"]:
            if pulser in list_of_landmarks:
                this_pulser = self.get_distant_pulsers(pulser, spice_depth)
                collect[pulser] = _resolve(pulser, this_pulser, R_distant, solution=solution)

        # South Pole landmarks

        for known_loc in ["ICL", "WT", "SPRESSO", "SPT"]:
            if known_loc in list_of_landmarks:
                this_loc = self.get_southpole_landmarks(known_loc)
                this_loc = np.array(this_loc)[0]
                st_centric = self.get_global_to_station_centric(this_loc[:3])
                collect[known_loc] = _resolve(known_loc, st_centric, R_distant, solution=solution)


        collect['critical_angle'] = self.get_critical_angle()
        return collect

