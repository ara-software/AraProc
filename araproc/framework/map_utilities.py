import ctypes
import numpy as np
import ROOT
from araproc.framework import constants


def get_corr_map_peak(the_map=None):
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

    _peakZ     = ctypes.c_int()
    _peakTheta = ctypes.c_int()
    _peakPhi   = ctypes.c_int()
    the_map.GetMaximumBin(_peakPhi, _peakTheta, _peakZ)

    peak_corr  = the_map.GetMaximum()
    peak_phi   = the_map.GetXaxis().GetBinCenter(_peakPhi.value)
    peak_theta = the_map.GetYaxis().GetBinCenter(_peakTheta.value)
    return peak_corr, peak_phi, peak_theta


def load_raytrace_model(
    arasimsrc="/cvmfs/icecube.opensciencegrid.org/users/abishop/AraSim-main"):
    """
    Load and initialize the AraSim ray-tracing model into ROOT's interpreter.

    Sets up the exponential refractive index model n(z) = A - B*exp(-C*|z|)
    with A=1.780, B=0.454 (= A - 1.326), C=0.0202 /m, and registers the
    rt_doTrace_savePoints function that saves all ray-path (x, z) points
    for sphere-crossing calculations. A, B, C parameters are same as in PA analysis paper

    Parameters
    ----------
    arasimsrc : str
        Path to the AraSim source directory containing RayTrace.h,
        RayTrace_IceModels.h, and Vector.h.
    """
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
        f'(new exponentialRefractiveIndex({1.3260},{1.780},{0.0202}));'
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
        raise RuntimeError("Failed to declare rt_doTrace_savePoints in ROOT.")


load_raytrace_model()


class AraGeom:

    """
    A class to initialize ARA geometry and provide associated utility functions.
    """

    def __init__(self, station_id):
        """
        Initialize the geometry tool and station information.

        Parameters
        ----------
        station_id : int
            ARA station number (e.g. 1, 2, 3, 4, 5).
        """
        self.FEET_TO_METER = 0.3048
        self.SUR_TO_GLOB   = np.array([22400, 53900, -36.02])

        self.station_id = station_id
        self.geomTool   = ROOT.AraGeomTool.Instance()
        self.st_info    = self.geomTool.getStationInfo(self.station_id)

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

        antenna_coordinates = [[] for _ in range(3)]
        for antenna_index in range(constants.num_rf_channels):
            if antenna_index in excluded_channels:
                continue
            antenna_info = self.st_info.getAntennaInfo(antenna_index)
            for axis in range(3):
                antenna_coordinates[axis].append(antenna_info.antLocation[axis])

        return tuple(np.average(coords) for coords in antenna_coordinates)

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

    def get_array_from_lat_long(self, latitude, longitude):
        """
        Convert latitude/longitude to ARA array coordinates.

        Parameters
        ----------
        latitude : float
            Latitude in decimal degrees.
        longitude : float
            Longitude in decimal degrees.

        Returns
        -------
        numpy.ndarray
            Array coordinates [northing, easting, 0], in meters.
        """

        easting  = self.geomTool.getEastingFromLatLong(latitude, longitude)
        northing = self.geomTool.getNorthingFromLatLong(latitude, longitude)
        return np.array([northing, easting, 0])

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
          Station-centric coordinates [x', y', z'] in meters.
        """

        input_coords  = np.array(array_coords, dtype=np.float64)
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
        elif landmark_type == "WT": ## This can go range(1,4)
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
            "IC1S":  [45659.6457, 50490.4199],
            "IC22S": [44884.416, 51444.8819],
            "SPIce": [42600, 48800],
        }
        pulser_depths = {"IC1S": -1400, "IC22S": -1450.47, "SPIce": spice_depth}

        if pulser_name not in pulser_coords:
            raise ValueError(f"Unknown pulser name: {pulser_name}")

        easting, northing = pulser_coords[pulser_name]
        global_coords = self.get_survey_to_global_coords(easting, northing) * self.FEET_TO_METER
        global_coords = [float(global_coords[0]), float(global_coords[1]), 0]
        st_centric    = self.get_global_to_station_centric(global_coords)
        return st_centric[0], st_centric[1], pulser_depths[pulser_name]

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


        delta    = np.array(target_point) - np.array(reference_point)
        x, y, z  = delta
        r        = np.sqrt(x**2 + y**2 + z**2)
        theta    = np.degrees(np.arccos(z / r))
        phi      = np.degrees(np.arctan2(y, x))
        return r, 90 - theta, phi

    def get_critical_angle(self):
        """
        Returns
        -------
        critical angle in terms of elevation angle

        """
        
        station_depth  = self.get_station_center()[2]
        n_air          = constants.get_index_of_refraction(1e-3)
        n_ice          = constants.get_index_of_refraction(station_depth)
        critical_angle = np.arcsin(n_air / n_ice)
        return 90 - np.degrees(critical_angle)

    def get_raytraced_sphere_crossing_angles(
        self,
        source_xyz,
        R_sphere,
        accuracy=0.001,
        solution=0,
        frequency=0.5,
        polarization=0.0,
        station_center=None,
    ):
        """
        Find the apparent direction of a source on a reconstruction sphere via ray tracing.

        The source is ray traced to the station center using AraSim's TraceFinder.
        The selected ray branch is then intersected with a sphere of radius R_sphere(depending on --map_type)
        centered on the station. This sphere crossing defines the apparent direction
        (elevation, azimuth) that should be overlaid on the sky map.

        If the traced path does not itself cross R_sphere (e.g. source is closer than
        R_sphere, or the path curves away), the outermost traced point and its local
        direction are used to extend the ray outward until it intersects the sphere
        (backward extension). If the result is still unphysical (elevation outside
        [-90, 90]), the straight-line direction is returned as a fallback.

        Ray solutions are identified by propagation-time ordering:
          - Shortest TOF  → label 'D' (direct)
          - Second TOF    → label 'R' (surface-reflected / refracted)

        Parameters
        ----------
        source_xyz : array-like
            Source position [x, y, z] in station-centric Cartesian coordinates, in meters.
        R_sphere : float
            Radius of the reconstruction sphere, in meters.
        accuracy : float, optional
            Ray-tracing convergence accuracy passed to AraSim findPaths. Default 0.001.
        solution : int, optional
            Requested solution type:
              - 0 : direct path (shortest propagation time)
              - 1 : reflected/refracted path (second-shortest propagation time)
            If solution=1 and findPaths returns only 1 solution, returns [] immediately
            so the caller can fall back to the straight-line direction.
        frequency : float, optional
            Signal frequency passed to AraSim ray tracer, in GHz. Default 0.5 GHz.
        polarization : float, optional
            Signal polarization passed to AraSim ray tracer, in radians. Default 0.0.
        station_center : array-like or None, optional
            Station center [x, y, z] in station-centric Cartesian coordinates, in meters.
            

        Returns
        -------
        list of dict
            Empty list if ray tracing fails or no valid sphere crossing is found.
            Otherwise a single-element list with keys:
              - 'label'       : str, 'D' or 'R'
              - 'elevation'   : float, elevation angle on R_sphere in degrees [-90, +90]
              - 'azimuth'     : float, azimuth angle in degrees [-180, +180]
              - 'tof'         : float, propagation time in nanoseconds
              - 'path_len'    : float, path length in meters
              - 'allowedUsed' : int, AraSim reflection mode flag
        """

        src = np.asarray(source_xyz, dtype=float)
        rec = np.asarray(station_center if station_center is not None else self.get_station_center(), dtype=float)
        R = float(R_sphere) # radius of the map

        z_src  = float(src[2])
        z_rec  = float(rec[2])

        dx     = float(src[0] - rec[0])
        dy     = float(src[1] - rec[1])
        r_src  = float(np.hypot(dx, dy))
        az_phi = float(np.arctan2(dy, dx))

        
        ROOT._rt_src.SetXYZ(r_src, 0.0, z_src)
        ROOT._rt_rec.SetXYZ(0.0,   0.0, z_rec)

        def _circle_intersect(p0, v, radius, t_max=np.inf):

            """
            Find the first intersection of p(t) = p0 + t*v with the
            station-centered circle x^2 + z^2 = radius^2 in the local
            propagation plane, where:
              - x >= 0 : radial distance from station, in meters
              - z      : vertical offset from station, in meters (positive upward)

            Parameterizes the path as p(t) = p0 + t*v and substitutes into
            the circle equation to get the quadratic a*t^2 + b*t + c = 0:
              a = |v|^2
              b = 2 * dot(p0, v)
              c = |p0|^2 - radius^2
            
            The quadratic has two roots t1 (smaller) and t2 (larger).
            t1 is the entry point into the circle, t2 is the exit point.
            We iterate smallest-first so we always return the nearest valid crossing.
            The two roots t1, t2 are the candidate intersections. Only roots
            with 0 <= t <= t_max  (lies on the segment (t_max=1) or forward ray (t_max=inf))
            and 
            x >= 0 are physically valid. (physically valid radial coordinate (radius cannot be negative))

            Used in two modes:
              - Segment (forward crossing):  t_max=1.0, v = p1 - p0
              - Ray (backward extension):    t_max=inf, v = local direction at p_far

            Parameters
            ----------
            p0 : tuple
                Start point (x, z) in meters.
            v : tuple
                Direction vector (vx, vz) in meters. For a segment, pass p1-p0.
            radius : float
                Circle radius in meters (= R_sphere).
            t_max : float, optional
                Upper bound on t. Default inf (semi-infinite ray).

            Returns
            -------
            tuple or None
                (x, z, t) of the first valid intersection in meters, or None.
            """
            
            x0, z0 = p0  # start point components
            vx, vz = v  # direction vector components
            a = vx*vx + vz*vz  # quadratic coefficient: |v|^2
            if a == 0.0:
                return None                         # zero-length direction vector, no intersection possible
            b    = 2.0*(x0*vx + z0*vz)              # quadratic coefficient: 2*dot(p0, v)
            c    = x0*x0 + z0*z0 - radius*radius    # quadratic coefficient: |p0|^2 - R^2
            disc = b*b - 4.0*a*c                    # discriminant: negative means no real intersection
            if disc < 0.0:
                return None                                 # ray misses the circle entirely
            sq = np.sqrt(disc)                               # square root of discriminant
            for t in ((-b - sq)/(2.0*a), (-b + sq)/(2.0*a)):   # t1 (near), t2 (far) roots, smallest first
                if 0.0 <= t <= t_max:                 # t must lie on segment (t_max=1) or forward ray (t_max=inf)
                    x, z = x0 + t*vx, z0 + t*vz       # intersection point
                    if x >= 0.0:                      # x >= 0 required: radial distance cannot be negative
                        return (x, z, t)              # return first valid crossing
            return None

        """
        findPaths from source to station: AraSim searches for all ray solutions in its 2D cylindrical propagation plane (x = horizontal radial distance
        from source, z = depth). Returns up to 2 solutions ordered by TOF:
          sol[0] -> direct path (D)   : shortest propagation time, travels directly
                                         through ice from source to station
          sol[1] -> reflected path (R): second-shortest, reflects off the ice surface
                                         before reaching the station
        """

        sol_cnt = ctypes.c_int()
        sol_err = ctypes.c_int()
        paths   = ROOT._rt_tf.findPaths(
            ROOT._rt_src, ROOT._rt_rec,
            float(frequency), ROOT.TMath.Pi() / 2,
            sol_cnt, sol_err,
            int(ROOT.RayTrace.SurfaceReflection), float(accuracy),
        )

        self._last_rt_npaths      = 0
        self._last_rt_paths       = []
        self._last_rt_fail_reason = None


        if sol_err.value != 0 or not paths or len(paths) == 0:
            self._last_rt_fail_reason = (
                f"findPaths error (serr={sol_err.value})"
                if sol_err.value != 0
                else "findPaths returned zero solutions"
            )
            return []


        sorted_paths = sorted(paths, key=lambda p: float(p.pathTime))

        self._last_rt_npaths = len(sorted_paths)
        self._last_rt_paths  = [
            {"tof": float(p.pathTime)*1e9, "path_len": float(p.pathLen),
             "launchAngle": float(p.launchAngle)}
            for p in sorted_paths
        ]



        # Label by TOF: shortest=D, second=R
        labeled = []
        if len(sorted_paths) >= 1: labeled.append((sorted_paths[0], "D"))
        if len(sorted_paths) >= 2: labeled.append((sorted_paths[1], "R"))

        if solution == 0:
            to_trace = [(p, l) for p, l in labeled if l == "D"]
        elif solution == 1:
            to_trace = [(p, l) for p, l in labeled if l == "R"]
            if not to_trace:
                # Only 1 solution from ray tracer — R simply does not exist
                #print(" [INFO] R solution not found: findPaths returned only 1 solution")
                self._last_rt_fail_reason = "R solution not found: only 1 ray-trace solution exists"
                return []
  

        allowed_used = int(ROOT.RayTrace.SurfaceReflection)

        # Straight-line fallback (used only if RT elevation is unphysical)

        _, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(rec, src)

        out = []
        for path, label in to_trace:
            sol_error = ctypes.c_int(0)
            ROOT.rt_doTrace_savePoints(
                z_src, float(path.launchAngle),
                ROOT.RayTrace.rayTargetRecord(z_rec, r_src),
                allowed_used, float(frequency), float(polarization), sol_error,
            )

            xs = np.asarray(list(ROOT.rt_get_x()), dtype=float)
            zs = np.asarray(list(ROOT.rt_get_z()), dtype=float)

            if sol_error.value != 0 or len(xs) < 2:
                print(f" [TRACE FAIL] label={label}  sol_error={sol_error.value}  npts={len(xs)}")
                continue

            # Receiver-centered coordinates in the local propagation plane
            xp       = r_src - xs       # radial distance from station, meters
            dz       = zs - z_rec       # vertical offset from station, meters
            d_to_rec = np.hypot(xp, dz) # distance from station along path, meters

            # Forward crossing: find where path crosses R_sphere outward from station
            # i_cross is the index where d steps from outside (d>=R) to inside (d<R)
            i_cross = next(
                (i for i in range(len(d_to_rec)-1, 0, -1)
                 if d_to_rec[i] <= R < d_to_rec[i-1]),
                None
            )

            if i_cross is not None:
                p_out = (float(xp[i_cross-1]), float(dz[i_cross-1])) # point just outside R_sphere
                p_in  = (float(xp[i_cross]),   float(dz[i_cross]))   # point just inside R_sphere
                inter = _circle_intersect( 
                    p_out, (p_in[0]-p_out[0], p_in[1]-p_out[1]), R, t_max=1.0)   # interpolate exact crossing on segment
                
            else:
                
                # Path never reaches R_sphere — source is inside R_sphere or path curves away.
                # Extend the ray from the outermost traced point using its local direction.
                imax          = int(np.argmax(d_to_rec))    # index of point farthest from station
                i_far, i_near = (0, 1) if imax == 0 else (imax, imax-1)     # far point and its neighbor for direction
                p_far         = (float(xp[i_far]),  float(dz[i_far]))       # outermost traced point
                p_near        = (float(xp[i_near]), float(dz[i_near]))      # neighbor point, defines local ray direction
                v_ext         = (p_far[0]-p_near[0], p_far[1]-p_near[1])    # local direction vector at p_far
                inter         = _circle_intersect(p_far, v_ext, R)        # extend ray until it hits R_sphere

            if inter is None:
                # Neither forward crossing nor backward extension reached R_sphere
                self._last_rt_fail_reason = (
                    f"no valid crossing/extension to R_sphere={R:.3f} m survived")
                continue

            xp_c, dz_c, _ = inter   # unpack crossing coordinates, discard t
            elev = np.degrees(np.arctan2(dz_c, xp_c)) # elevation angle at sphere crossing, degrees
            az   = np.degrees(az_phi)                       # azimuth angle

            if not (-90.0 <= elev <= 90.0):
                elev, az = elev_sl, az_sl   # unphysical elevation — fall back to straight-line

            out.append({
                "label"      : label,
                "elevation"  : float(elev),
                "azimuth"    : float(az),
                "tof"        : float(path.pathTime) * 1e9,
                "path_len"   : float(path.pathLen),
                "allowedUsed": allowed_used,
            })

        if not out:
            self._last_rt_fail_reason = (
                f"no valid crossing/extension to R_sphere={R:.3f} m survived")

        return out

    
    def get_raytraced_critical_angle(self, R_sphere):

        """
        Map the critical reception angle onto the reconstruction sphere via ray tracing.

        Integrates the critical-angle ray upward from the station through the ice
        using the Snell invariant n(z)*sin(theta) = n_air. The refractive index
        n(z) is evaluated at each step via constants.get_index_of_refraction(z)
        — no hardcoded ice model parameters.

        Two successful exit conditions:
          1. Ray crosses R_sphere within ice: returns the crossing elevation.
          2. Ray exits ice (z=0) before reaching R_sphere: returns the elevation
             of the surface-exit point. This occurs for large R_sphere (e.g. 300 m)
             where the ray reaches the surface before the sphere boundary.

        Falls back to get_critical_angle() (geometric, no ray bending) if:
          - sin(theta) >= 1 is encountered during integration, meaning the ray
            cannot propagate upward from that depth.
          - The sphere-crossing elevation is outside [-90, 90] degrees.
          - The integration loop exits without finding any crossing.

        Parameters
        ----------
        R_sphere : float
            Reconstruction sphere radius, in meters.

        Returns
        -------
        float
            Elevation angle in degrees where the critical-angle ray appears on
            the reconstruction sphere. Range: [0, 90].
            Returns get_critical_angle() on any failure.
        """

        z_rec  = float(self.get_station_center()[2])
        R = float(R_sphere)
        n_air  = constants.get_index_of_refraction(1e-3)
        #n_ice  = constants.get_index_of_refraction(z_rec)

        dz_step = 0.1  # integration step, meters 
        z_prev, r_prev, d_prev = z_rec, 0.0, 0.0

        while z_prev < 0.0:
            z_curr = z_prev + dz_step
            n_z    = constants.get_index_of_refraction(z_prev + 0.5*dz_step)
            sin_th = n_air / n_z

            if sin_th >= 1.0:
                return float(self.get_critical_angle())

            r_curr  = r_prev + (sin_th / np.sqrt(1.0 - sin_th**2)) * dz_step
            dz_curr = z_curr - z_rec
            d_curr  = np.hypot(r_curr, dz_curr)

 
            # Case 1: ray crosses R_sphere while still in ice
            if d_prev <= R < d_curr:
                frac     = (R - d_prev) / (d_curr - d_prev)          # fraction along step where d = R exactly
                r_cross  = r_prev + frac*(r_curr - r_prev)            # horizontal radial distance at crossing, meters
                dz_cross = (z_prev + frac*(z_curr - z_prev)) - z_rec  # vertical offset from station at crossing, meters
                elev     = float(np.degrees(np.arctan2(dz_cross, r_cross)))  # elevation angle at crossing, degrees
                if not (-90.0 <= elev <= 90.0):
                    return float(self.get_critical_angle())
                return elev

            # Case 2: ray exits ice surface (z=0) before reaching R_sphere
            if z_prev <= 0.0 <= z_curr:
                frac   = (0.0 - z_prev) / (z_curr - z_prev)  # fraction along step where z = 0 exactly
                r_surf = r_prev + frac*(r_curr - r_prev)      # horizontal radial distance at surface exit, meters
                elev   = float(np.degrees(np.arctan2(-z_rec, r_surf)))  # elevation using station depth as vertical leg, degrees
                if not (-90.0 <= elev <= 90.0):
                    return float(self.get_critical_angle())
                return elev

            z_prev, r_prev, d_prev = z_curr, r_curr, d_curr

        return float(self.get_critical_angle()) # return the default critical angle if crossing case 1 and 2 fail



    def get_known_landmarks(
        self,
        list_of_landmarks=None,
        R_map=None,
        list_of_cal_pulser_indices=None,
        spice_depth=None,
        solution=None,
    ):
        """
        Compute apparent directions of known landmarks on a reconstruction sphere.

        Each landmark is ray traced to the station and its direction is mapped onto
        the reconstruction sphere of radius R_map. The ray tracer (AraSim findPaths)
        always returns at least 1 solution for all landmarks and all stations tested.
        Solutions are ordered by propagation time (TOF):
          - Shortest TOF  -> Direct path (D)  : travels straight through ice
          - Second TOF    -> Reflected path (R): reflects off the ice surface

        Which solution is shown on the map depends on the map type:
          - direct maps   (e.g. pulser_*, distant_*_dir) : solution=0, shows D
          - reflected maps (e.g. distant_*_ref)           : solution=1, shows R

        Fallback rules:
          - D solution: shown on map unless the traced path never crosses R_map.
            This can happen if the source is closer to the station than R_map
            AND the backward extension of the path also fails to reach R_map
            (e.g. source inside R_map with a nearly vertical path). Falls back
            to straight-line (SL) direction in that case.

          - R solution: shown on map if findPaths returns 2 solutions AND the
            reflected path crosses R_map. Falls back to SL if:
              (a) findPaths returns only 1 solution — the reflected path
                  physically does not exist for this source-station geometry, OR
              (b) the reflected path never crosses R_map — same geometric
                  condition as the D fallback above.

        Parameters
        ----------
        list_of_landmarks : list of str, optional
            Landmark names to include. Supported values:
            'ICL', 'IC22S', 'SPT', 'IC1S', 'SPIce', 'WT', 'SPRESSO'.
            Pass ['all'] to include all of the above.
            Default: ['IC22S', 'ICL', 'SPRESSO'].
        R_map : float
            Reconstruction sphere radius corresponding to the map, in meters.
        list_of_cal_pulser_indices : list of int, optional
            Calpulser indices to include, e.g. [0, 1, 2, 3].
            Pass ['all'] for all four. Default: [1, 3].
        spice_depth : float, optional
            Depth of the SPIce pulser in meters. Required if 'SPIce' is in
            list_of_landmarks.
        solution : int
            Requested ray solution:
              - 0 : direct path (shortest propagation time), for direct maps
              - 1 : reflected path (second-shortest propagation time), for reflected maps

        Returns
        -------
        collect : dict
            Keys are landmark name strings. Values are [r, elevation, azimuth]
            where r is the straight-line distance in meters, elevation and
            azimuth are the ray-traced (or SL fallback) angles in degrees.
            Also includes:
              - 'critical_angle'    : geometric critical elevation angle, in degrees
              - 'critical_angle_rt' : ray-traced critical elevation angle mapped
                                      onto R_map, in degrees
        """
        

        if list_of_landmarks is None:
            list_of_landmarks = ["IC22S", "ICL", "SPRESSO"]
        elif list_of_landmarks == ["all"]:
            list_of_landmarks = ["ICL", "IC22S", "SPT", "IC1S", "WT", "SPRESSO"]

        if list_of_cal_pulser_indices is None:
            list_of_cal_pulser_indices = [1, 3]
        elif list_of_cal_pulser_indices == ["all"]:
            list_of_cal_pulser_indices = [0, 1, 2, 3]

        if spice_depth is not None and "SPIce" not in list_of_landmarks:
            list_of_landmarks.append("SPIce")

        collect = {}
        station_center = self.get_station_center()
        

        def _resolve(name, source_xyz, station_center, R_sphere, solution):

            """
            Resolve one landmark to its apparent direction on the reconstruction sphere.

            Computes the straight-line (SL) direction first — this is always valid
            and serves as the fallback. Then attempts ray tracing via
            get_raytraced_sphere_crossing_angles(). Falls back to SL with diagnostic
            output if ray tracing fails for any of the reasons documented in
            get_known_landmarks().


            Parameters
            ----------
            name : str
                Landmark label used for diagnostic printing and as the collect key.
            source_xyz : array-like
                Source position [x, y, z] in station-centric Cartesian coordinates,
                in meters.
            station_center : array-like
                Station center [x, y, z] in station-centric Cartesian coordinates,
                in meters. 

            R_sphere : float
                Reconstruction sphere radius, in meters.

            solution : int
                0 = direct path, 1 = reflected path. See get_known_landmarks().

            Returns
            -------
            list
                [r, elevation, azimuth] where r is the SL distance in meters,
                and elevation and azimuth are in degrees. Ray-traced if successful,
                straight-line fallback otherwise.

            """
            r_sl, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(
                station_center, source_xyz
            )


            rt_sols = self.get_raytraced_sphere_crossing_angles(
                source_xyz    = source_xyz,
                R_sphere      = R_sphere,
                solution      = solution,
                station_center= station_center,
            )

            if not rt_sols:
                return [r_sl, elev_sl, az_sl] # Fallback to SL solution
            
            return [r_sl, rt_sols[0]["elevation"], rt_sols[0]["azimuth"]] #always the requested solution, either D or R.


        
        # Local Calpulsers

        calpulser = self.get_local_CP(list_of_cal_pulser_indices)
        for cp in calpulser:
            collect[f"CP{cp}"] = _resolve(f"CP{cp}", calpulser[cp], station_center, R_map, solution)
        del calpulser
        
        # Distant pulsers

        for pulser in ["IC22S", "IC1S", "SPIce"]:
            if pulser in list_of_landmarks:
                this_pulser = self.get_distant_pulsers(pulser, spice_depth)
                collect[pulser] = _resolve(pulser, this_pulser, station_center, R_map, solution)

        # South Pole surface landmarks

        for known_loc in ["ICL", "WT", "SPRESSO", "SPT"]:
            if known_loc in list_of_landmarks:
                this_loc   = np.array(self.get_southpole_landmarks(known_loc))[0]
                st_centric = self.get_global_to_station_centric(this_loc[:3])
                if known_loc == "SPT":
                    st_centric[2] = 7.55984245e-02  # corrected SPT depth, meters, set to ICL depth 
                collect[known_loc] = _resolve(known_loc, st_centric, station_center, R_map, solution)

        
        # Critical angle (ray traced and fall back SL solution)

        collect["critical_angle_rt"] = self.get_raytraced_critical_angle(R_map)

        return collect
