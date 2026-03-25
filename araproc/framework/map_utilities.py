import ctypes
import numpy as np
import ROOT
import os
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



def load_raytrace_model(arasimsrc=None):
    """
    Load and initialize the AraSim ray-tracing model into ROOT's interpreter.

    Sets up the exponential refractive index model n(z) = A - B*exp(-C*|z|)
    with A=1.780, B=0.454 (= A - 1.326), C=0.0202 /m, and registers the
    rt_doTrace_savePoints function that saves all ray-path (x, z) points
    for reconstruction-sphere crossing calculations.

    Parameters
    ----------
    arasimsrc : str or None
        Path to the AraSim source directory containing RayTrace.h,
        RayTrace_IceModels.h, and Vector.h. If None, the path is taken
        from the ARA_SIM_DIR environment variable.
    """
    if arasimsrc is None:
        arasimsrc = os.getenv("ARA_SIM_DIR", None)

    if not arasimsrc:
        raise RuntimeError(
            "ARA_SIM_DIR environment variable not set. "
            "Please set it to the AraSim source directory."
        )

    arasimsrc = os.path.abspath(arasimsrc)

    # Load the AraSim ray-tracing headers into ROOT's interpreter so the
    # corresponding C++ classes can be used from Python.

    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/RayTrace_IceModels.h"')
    ROOT.gInterpreter.ProcessLine(f'#include "{arasimsrc}/Vector.h"')
    ROOT.gInterpreter.ProcessLine('#include <vector>')

    # Instantiate the attenuation model and the exponential refractive-index
    # model used by the ray tracer. The refractive-index parameters are the
    # same as those used in the PA analysis paper.

    ROOT.gInterpreter.ProcessLine(
        'auto _rt_atten = boost::shared_ptr<basicAttenuationModel>'
        '(new basicAttenuationModel);'
    )
    ROOT.gInterpreter.ProcessLine(
        'auto _rt_refr = boost::shared_ptr<exponentialRefractiveIndex>'
        '(new exponentialRefractiveIndex(1.3260,1.780,0.0202));'
    )

    # Create the AraSim ray-tracing engine from the refractive-index and
    # attenuation models above

    ROOT.gInterpreter.ProcessLine('RayTrace::TraceFinder _rt_tf(_rt_refr, _rt_atten);')

    # Declare source/receiver vectors used by the ray-tracing machinery.

    ROOT.gInterpreter.ProcessLine('Vector _rt_src; Vector _rt_rec;')

    # Declare helper C++ code inside ROOT.
    #
    # This block:
    #  1. defines persistent vectors to store the traced ray path
    #     (x, z, and RK step type),
    #  2. defines a callback functor that is called at each tracing step
    #     and saves the current ray position,
    #  3. wraps the AraSim doTrace call in a helper function
    #     (rt_doTrace_savePoints) that clears old path data, runs the trace,
    #     and records all intermediate path points, and
    #  4. exposes getter functions so Python can retrieve the saved path.
    #
    # We need this because for reconstruction-sphere mapping we do not only
    # need the final ray-trace solution; we also need the intermediate traced
    # path points to determine where the selected ray branch crosses the sphere.
    ROOT.gInterpreter.Declare(r'''
        #include <vector>

        // Traced x positions along the ray path.
        static std::vector<double> _rt_x;

        // Traced z positions along the ray path.
        static std::vector<double> _rt_z;

        // RK step classification saved at each traced point.
        static std::vector<int> _rt_stepType;

        // Callback invoked at each ray-tracing step to save the current
        // ray position and step type.
        struct RTPointSaver {
            void operator()(const RayTrace::fullRayPosition& p, RayTrace::RKStepType t) const {
                _rt_x.push_back(p.x);
                _rt_z.push_back(p.z);
                _rt_stepType.push_back((int)t);
            }
        };

        // Wrapper around AraSim doTrace that clears any previously saved
        // path information, runs the trace, and records all intermediate
        // path points through the callback above.
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

        // Accessors used from Python to retrieve the saved traced path.
        const std::vector<double>& rt_get_x(){ return _rt_x; }
        const std::vector<double>& rt_get_z(){ return _rt_z; }
        const std::vector<int>& rt_get_stepType(){ return _rt_stepType; }
    ''')

    # Confirm that the helper function was successfully declared in ROOT.
    if not hasattr(ROOT, "rt_doTrace_savePoints"):
        raise RuntimeError("Failed to declare rt_doTrace_savePoints in ROOT.")

#Load Model 

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
           SPT : South Pole Telescope (removed from landmark list due to depth bug in AraRoot)

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
        #elif landmark_type == "SPT":
        #    corner = self.geomTool.getSouthPoleTelescope()
        #    landmarks.append([corner[0], corner[1], corner[2]])
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




    def _circle_intersect(self, p0, v, radius, t_max=1.0):
        """
        Find the first intersection of a straight line segment (or ray) with the
        station-centered circle x^2 + z^2 = radius^2 in the local propagation plane.

        Here, p(t) = p0 + t*v parameterizes only the local straight segment used for
        interpolation between two successive sampled ray-trace points. It does NOT
        represent the full curved ray-trace path. The full path comes from AraSim;
        this helper is used only to determine where that discretely sampled path
        crosses the reconstruction circle more accurately than choosing the sampled
        point whose radius is merely closest to R_sphere.

        In the local propagation plane:
          - x >= 0 : radial distance from the station, in meters
          - z      : vertical offset from the station, in meters (positive upward)

        Substituting p(t) into x^2 + z^2 = radius^2 gives:
          a*t^2 + b*t + c = 0
        with
          a = |v|^2
          b = 2 * dot(p0, v)
          c = |p0|^2 - radius^2

        The two roots are checked in increasing order so the nearest valid crossing
        is returned first. Only roots with 0 <= t <= t_max are accepted.

        Parameters
        ----------
        p0 : tuple of float
            Start point (x, z), in meters.
        v : tuple of float
            Direction vector (vx, vz) in the local (x, z) plane. For a finite
            segment, pass p1 - p0.
        radius : float
            Circle radius, in meters (= R_sphere or R_map).
        t_max : float, optional
            Upper bound on t. In this file, t_max = 1.0 is used for a finite
            segment between successive sampled points.

        Returns
        -------
        tuple or None
            (x, z, t) for the first valid intersection, with x and z in meters,
            or None if no valid intersection exists.
        """
        x0, z0 = p0
        vx, vz = v

        a = vx * vx + vz * vz
        if a == 0.0:
            return None  # zero-length segment / direction vector

        b = 2.0 * (x0 * vx + z0 * vz)
        c = x0 * x0 + z0 * z0 - radius * radius
        disc = b * b - 4.0 * a * c
        if disc < 0.0:
            return None  # no real intersection with the circle

        sq = np.sqrt(disc)
        for t in ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)):
            if 0.0 <= t <= t_max:
                x = x0 + t * vx
                z = z0 + t * vz
                if x >= 0.0:  # radial distance from station must be non-negative
                    return (x, z, t)

        return None


    def _find_sphere_crossing_from_samples(
        self,
        xs,
        zs,
        z_ref,
        radius,
        reverse_search=False,
    ):
        """
        Find the intersection of sampled ray-trace points with the station-centered
        circle x^2 + z^2 = radius^2 in the local propagation plane.

        This helper is shared by:
          - get_raytraced_sphere_crossing_angles()
          - get_raytraced_critical_angle()

        It finds the first pair of consecutive sampled points that bracket the desired radius and
        then interpolates the crossing on that finite segment using
        self._circle_intersect().

        Parameters
        ----------
        xs : array-like
            Sampled horizontal coordinates in the local propagation plane, in meters.
        zs : array-like
            Sampled vertical coordinates in the local propagation plane, in meters.
        z_ref : float
            Reference station depth, in meters, used to convert z samples into
            station-centered vertical offsets dz = z - z_ref.
        radius : float
            Target circle / sphere radius, in meters.
        reverse_search : bool, optional
            If True, search from the end of the sampled path backward. This is used
            for source->station traces, where the desired crossing is the outermost
            one encountered before arriving at the station (for finding crossing point with the radius map for source
            within the radius map. 

            Example:
                For a source outside the map radius (e.g., r = 1000 m) and
                R_sphere = 300 m, the sampled distances along the ray may look like:
                    1000 → 800 → 500 → 350 → 280 → 100 → 0 (station)
                The crossing of R_sphere occurs between 350 → 280.
                We search backward to select this crossing before reaching
                the station.

            If False, search from the beginning forward. This is used for
            station->outward traces (for finding ray-traced critical angle)
            Example:
                For an outward trace starting at the station (r = 40m) with
                R_sphere = 300 m, the sampled distances may look like:
                    40 → 100 → 150 → 280 → 350 → 500
                The crossing of R_sphere occurs between 280 → 350.
                We search forward to select this first crossing after leaving
                the station.

        Returns
        -------
        tuple
            (inter, i_cross), where:
              - inter is (x_cross, z_cross, t) from self._circle_intersect(), or None
              - i_cross is the bracketing segment index, or None if no crossing exists
        """
        xp = np.asarray(xs, dtype=float)
        dz = np.asarray(zs, dtype=float) - float(z_ref)

        if len(xp) < 2:
            return None, None

        d = np.hypot(xp, dz)

        if reverse_search:
            # Samples run from source toward station: outside -> inside.
            # We want the last crossing before reaching the station.
            i_cross = next(
                (i for i in range(len(d) - 1, 0, -1) if d[i] <= radius < d[i - 1]),
                None,
            )
            if i_cross is None:
                return None, None

            p0 = (float(xp[i_cross - 1]), float(dz[i_cross - 1]))  # outside
            p1 = (float(xp[i_cross]),     float(dz[i_cross]))      # inside

        else:
            # Samples run outward from the station: inside -> outside.
            # We want the first crossing after leaving the station.
            i_cross = next(
                (i for i in range(1, len(d)) if d[i - 1] <= radius < d[i]),
                None,
            )
            if i_cross is None:
                return None, None

            p0 = (float(xp[i_cross - 1]), float(dz[i_cross - 1]))  # inside
            p1 = (float(xp[i_cross]),     float(dz[i_cross]))      # outside

        inter = self._circle_intersect(
            p0,
            (p1[0] - p0[0], p1[1] - p0[1]),
            float(radius),
            t_max=1.0,
        )

        return inter, i_cross


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
        The selected ray branch is intersected with the reconstruction sphere of
        radius R_sphere centered on the station. That sphere crossing defines the
        apparent direction (elevation, azimuth) that should be overlaid on the sky map.

        If the source->station traced samples do not themselves cross R_sphere
        (for example because the source lies inside R_sphere, or because the sampled
        branch reaches the station before crossing the sphere), the same branch is
        re-traced outward from the station using the branch receiptAngle. The sphere
        crossing is then determined from those outward ray-traced samples.

        If the mapped elevation lies outside the physical plotting range [-90, +90]
        degrees, the straight-line (SL) direction is used only as a plotting fallback
        and is flagged accordingly.

        Ray solutions are identified only by propagation-time ordering:
          - shortest TOF  -> label 'D' (direct)
          - second TOF    -> label 'R' (reflected / refracted)

        Parameters
        ----------
        source_xyz : array-like
            Source position [x, y, z] in station-centric Cartesian coordinates,
            in meters.
        R_sphere : float
            Reconstruction sphere radius, in meters. Must be positive.
        accuracy : float, optional
            Convergence tolerance passed to AraSim findPaths (dimensionless).
            Default is 0.001.
        solution : int, optional
            Requested solution type:
              - 0 : direct path (shortest propagation time)
              - 1 : reflected / refracted path (second-shortest propagation time)
        frequency : float, optional
            Signal frequency passed to the AraSim ray tracer, in GHz.
            Default is 0.5 GHz.
        polarization : float, optional
            Signal polarization passed to the AraSim ray tracer, in radians.
            Default is 0.0.
        station_center : array-like or None, optional
            Station center [x, y, z] in station-centric Cartesian coordinates,
            in meters. If None, the internally stored station center is used.

        Returns
        -------
        list of dict
            Each dict corresponds to one valid ray-traced solution and contains:
              - 'label'            : 'D' or 'R'
              - 'elevation'        : elevation angle on R_sphere, in degrees
              - 'azimuth'          : azimuth angle, in degrees
              - 'tof'              : propagation time, in nanoseconds
              - 'path_len'         : ray path length, in meters
              - 'allowedUsed'      : AraSim reflection-mode flag
              - 'used_sl_fallback' : bool, True only if plotting fell back to SL
        """
        if R_sphere is None:
            raise ValueError("R_sphere must be provided to get_raytraced_sphere_crossing_angles().")
        if R_sphere <= 0:
            raise ValueError(f"R_sphere must be positive, got {R_sphere}.")

        src = np.asarray(source_xyz, dtype=float)
        rec = np.asarray(
            station_center if station_center is not None else self.get_station_center(),
            dtype=float,
        )
        R = float(R_sphere)

        z_src = float(src[2])
        z_rec = float(rec[2])

        dx = float(src[0] - rec[0])
        dy = float(src[1] - rec[1])
        r_src = float(np.hypot(dx, dy))
        az_phi = float(np.arctan2(dy, dx))

        ROOT._rt_src.SetXYZ(r_src, 0.0, z_src)
        ROOT._rt_rec.SetXYZ(0.0, 0.0, z_rec)

        # findPaths searches in AraSim's 2D cylindrical propagation plane and
        # returns the available ray solutions. We later label them only by TOF:
        # shortest -> D, second-shortest -> R.
        sol_cnt = ctypes.c_int()
        sol_err = ctypes.c_int()
        paths = ROOT._rt_tf.findPaths(
            ROOT._rt_src,
            ROOT._rt_rec,
            float(frequency),
            ROOT.TMath.Pi() / 2,
            sol_cnt,
            sol_err,
            int(ROOT.RayTrace.SurfaceReflection),
            float(accuracy),
        )

        self._last_rt_npaths = 0
        self._last_rt_paths = []
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
        self._last_rt_paths = [
            {
                "tof": float(p.pathTime) * 1e9,
                "path_len": float(p.pathLen),
                "launchAngle": float(p.launchAngle),
            }
            for p in sorted_paths
        ]

        labeled = []
        if len(sorted_paths) >= 1:
            labeled.append((sorted_paths[0], "D"))
        if len(sorted_paths) >= 2:
            labeled.append((sorted_paths[1], "R"))

        if solution == 0:
            to_trace = [(p, l) for p, l in labeled if l == "D"]
        elif solution == 1:
            to_trace = [(p, l) for p, l in labeled if l == "R"]
            if not to_trace:
                self._last_rt_fail_reason = "R solution not found: only 1 ray-trace solution exists"
                return []
        else:
            raise ValueError(f"solution must be 0 or 1, got {solution}.")

        allowed_used = int(ROOT.RayTrace.SurfaceReflection)

        # Straight-line direction is used only as a last-resort plotting fallback
        # if the mapped RT elevation lies outside [-90, +90] degrees.
        _, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(rec, src)

        out = []
        for path, label in to_trace:
            used_sl_fallback = False

            # First trace: source -> station along the chosen TOF-ordered branch.
            sol_error = ctypes.c_int(0)
            ROOT.rt_doTrace_savePoints(
                z_src,
                float(path.launchAngle),
                ROOT.RayTrace.rayTargetRecord(z_rec, r_src),
                allowed_used,
                float(frequency),
                float(polarization),
                sol_error,
            )

            xs = np.asarray(list(ROOT.rt_get_x()), dtype=float)
            zs = np.asarray(list(ROOT.rt_get_z()), dtype=float)

            if sol_error.value != 0 or len(xs) < 2:
                print(f"[TRACE FAIL] label={label}  sol_error={sol_error.value}  npts={len(xs)}")
                continue

            # Convert source->station samples into receiver-centered coordinates.
            # xp_forward is radial distance from station in meters.
            xp_forward = r_src - xs

            # For source->station samples, search from the end backward because
            # the path terminates at the station and we want the outer sphere crossing.
            inter, i_cross = self._find_sphere_crossing_from_samples(
                xp_forward,
                zs,
                z_rec,
                R,
                reverse_search=True,
            )

            if inter is None:
                print(
                    f"[BACKTRACE INFO] label={label}  "
                    f"n_ray_solutions={self._last_rt_npaths}  "
                    f"forward_trace_found_no_crossing=True  "
                    f"R_sphere={R:.3f} m  "
                    f"source_r={r_src:.3f} m  "
                    f"source_z={z_src:.3f} m"
                )

                if not hasattr(path, "receiptAngle"):
                    self._last_rt_fail_reason = (
                        f"ray-traced solution exists for label={label}, but no forward crossing "
                        f"was found for R_sphere={R:.3f} m and receiptAngle is unavailable "
                        f"for backward tracing"
                    )
                    print(f"[BACKTRACE FAIL] {self._last_rt_fail_reason}")
                    continue

                # Re-trace the same branch outward from the station using AraSim,
                
                sol_error_back = ctypes.c_int(0)
                r_back_target = max(float(R) + 50.0, float(r_src) + 50.0)

                ROOT.rt_doTrace_savePoints(
                    z_rec,
                    float(path.receiptAngle),
                    ROOT.RayTrace.rayTargetRecord(float(z_src), r_back_target),
                    allowed_used,
                    float(frequency),
                    float(polarization),
                    sol_error_back,
                )

                xs_back = np.asarray(list(ROOT.rt_get_x()), dtype=float)
                zs_back = np.asarray(list(ROOT.rt_get_z()), dtype=float)

                if sol_error_back.value != 0 or len(xs_back) < 2:
                    self._last_rt_fail_reason = (
                        f"backward trace failed for label={label}: "
                        f"sol_error={sol_error_back.value}, npts={len(xs_back)}"
                    )
                    print(f"[BACKTRACE FAIL] {self._last_rt_fail_reason}")
                    continue

                inter, i_cross_back = self._find_sphere_crossing_from_samples(
                    xs_back,
                    zs_back,
                    z_rec,
                    R,
                    reverse_search=False,
                )

                if inter is not None:
                    print(
                        f"[BACKTRACE OK] label={label}  "
                        f"R_sphere={R:.3f} m  "
                        f"crossing_index={i_cross_back}  "
                        f"npts={len(xs_back)}"
                    )
                else:
                    self._last_rt_fail_reason = (
                        f"ray-traced solution found for label={label} "
                        f"(n_ray_solutions={self._last_rt_npaths}), but neither forward "
                        f"nor backward traced samples crossed R_sphere={R:.3f} m"
                    )
                    print(f"[BACKTRACE FAIL] {self._last_rt_fail_reason}")
                    continue

            xp_c, dz_c, _ = inter
            elev = float(np.degrees(np.arctan2(dz_c, xp_c)))
            az = float(np.degrees(az_phi))

            if not (-90.0 <= elev <= 90.0):
                self._last_rt_fail_reason = (
                    f"ray-traced elevation {elev:.3f} deg outside [-90, 90]; "
                    f"used straight-line fallback for plotting"
                )
                elev, az = elev_sl, az_sl
                used_sl_fallback = True

            out.append(
                {
                    "label": label,
                    "elevation": float(elev),
                    "azimuth": float(az),
                    "tof": float(path.pathTime) * 1e9,
                    "path_len": float(path.pathLen),
                    "allowedUsed": allowed_used,
                    "used_sl_fallback": used_sl_fallback,
                }
            )

        if not out and not self._last_rt_fail_reason:
            self._last_rt_fail_reason = f"no valid crossing to R_sphere={R:.3f} m survived"

        return out



    def get_raytraced_critical_angle(
        self,
        R_sphere,
        frequency=0.5,
        polarization=0.0,
    ):
        """
        Map the critical reception angle onto the reconstruction sphere using the
        AraSim ray tracer.

        The critical angle is first defined locally at the station via
        get_critical_angle(), which returns the geometric critical elevation
        in degrees. That local critical direction is then traced outward with
        AraSim, and the apparent critical angle on the sky map is defined by
        where that ray-traced trajectory crosses the reconstruction sphere of
        radius R_sphere centered on the station.

        The initial launch direction is obtained from get_critical_angle(),
        which computes the geometric critical angle at the station depth using
        the local index of refraction. This provides a first-order estimate of
        the direction from which a critically refracted ray would arrive at the
        station in a straight-line (locally uniform medium) approximation.

        This direction is then used to initialize the AraSim ray tracer at the
        station depth. The full ray-traced trajectory accounts for the depth-
        dependent refractive index profile n(z), allowing the curved ray path
        corresponding to that critical condition to be followed through the ice.

        The ray is launched from the station rather than from the surface. This
        relies on ray-path reversibility: the outward-traced trajectory is the
        time-reversed equivalent of a ray arriving at the station at the critical
        angle. Using the station depth as the starting point ensures a well-defined
        initial condition for the ray tracer and maintains consistency with how
        other ray-traced directions are constructed in this module.

        This function uses the same sampled sphere-crossing logic as
        get_raytraced_sphere_crossing_angles():
          1. Ray trace the trajectory outward from the station.
          2. Search the traced (x, z) samples for a crossing of R_sphere.
          3. Interpolate the crossing on the bracketing segment with
             self._circle_intersect() via the shared helper
             self._find_sphere_crossing_from_samples().

        Unlike the source-mapping case, no source->station ray-tracing step is
        required. The critical-angle trajectory is defined locally at the station,
        so it can be traced directly outward from the station. Therefore only a
        single outward ray trace is needed, and no fallback re-tracing is required.

        The ray is traced continuously through the ice-air interface; reaching
        z = 0 is not treated as a terminal condition. To ensure this, the outward
        target is placed slightly above the surface.

        Falls back to get_critical_angle() (geometric, no ray bending) only if
        no valid ray-traced crossing on R_sphere can be determined, or if the
        mapped ray-traced elevation is outside the physical plotting range
        [-90, +90] degrees.

        Parameters
        ----------
        R_sphere : float
            Reconstruction sphere radius, in meters. Must be positive.
        frequency : float, optional
            Signal frequency passed to the AraSim ray tracer, in GHz.
            Default is 0.5 GHz.
        polarization : float, optional
            Signal polarization passed to the AraSim ray tracer, in radians.
            Default is 0.0.

        Returns
        -------
        float
            Elevation angle, in degrees, where the ray-traced critical-angle
            trajectory appears on the reconstruction sphere. Returns
            get_critical_angle() only if no valid RT-based crossing can be obtained.
        """
        if R_sphere is None:
            raise ValueError("R_sphere must be provided to get_raytraced_critical_angle().")
        if R_sphere <= 0:
            raise ValueError(f"R_sphere must be positive, got {R_sphere}.")

        z_rec = float(self.get_station_center()[2])
        R = float(R_sphere)

        # Geometric critical angle (straight line, no bending) approximation evaluated at the station depth using the local
        # index of refraction that defines the direction from which a critically refracted ray would
        # arrive at the station.
        crit_elev_deg = float(self.get_critical_angle())

        # Convert the local elevation to the angle convention required by the
        # AraSim ray tracer. This direction is used as the initial launch angle
        # at the station for the full ray-traced propagation in n(z).
        crit_rt_angle = np.deg2rad(90.0 - crit_elev_deg)

        allowed_used = int(ROOT.RayTrace.SurfaceReflection)

        # Initialize the ray tracer at the station depth (z_rec). Although the
        # critical angle is defined for a ray approaching the surface from below,
        # we launch the ray outward from the station using the same direction.
        # By ray-path reversibility, this outward trajectory is equivalent to
        # the time-reversed incoming critical ray and therefore traces the same
        # physical path through the medium.
        #
        # The target is placed slightly above the surface to ensure that the ray
        # is allowed to continue through the ice-air interface rather than
        # effectively stopping at z = 0.
        r_target = float(R + 50.0)
        z_target = 0.050  # in meters; slightly above the ice surface

        sol_error = ctypes.c_int(0)
        ROOT.rt_doTrace_savePoints(
            z_rec,
            float(crit_rt_angle),
            ROOT.RayTrace.rayTargetRecord(z_target, r_target),
            allowed_used,
            float(frequency),
            float(polarization),
            sol_error,
        )

        xs = np.asarray(list(ROOT.rt_get_x()), dtype=float)
        zs = np.asarray(list(ROOT.rt_get_z()), dtype=float)

        if sol_error.value != 0 or len(xs) < 2:
            print(
                f"[CRITICAL ANGLE RT FAIL] forward trace failed: "
                f"sol_error={sol_error.value}, npts={len(xs)}"
            )
            return float(self.get_critical_angle())

        # Outward traces from the station are ordered station -> outward,
        # so we search forward for the first inside->outside crossing.
        inter, i_cross = self._find_sphere_crossing_from_samples(
            xs,
            zs,
            z_rec,
            R,
            reverse_search=False,
        )

        if inter is None:
            print(
                f"[CRITICAL ANGLE RT FAIL] no valid RT crossing to "
                f"R_sphere={R:.3f} m  critical_elevation={crit_elev_deg:.3f} deg"
            )
            return float(self.get_critical_angle())

        x_cross, z_cross, _ = inter
        elev = float(np.degrees(np.arctan2(z_cross, x_cross)))

        if not (-90.0 <= elev <= 90.0):
            print(
                f"[CRITICAL ANGLE RT WARN] mapped elevation {elev:.3f} deg outside [-90, 90]; "
                f"falling back to geometric critical angle with Straight Line approximation"
            )
            return float(self.get_critical_angle())

        return elev


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
        returns the available solutions, which are labeled only by propagation-time
        ordering:
          - shortest TOF  -> Direct path (D)
          - second TOF    -> Reflected / refracted path (R)

        Which solution is shown depends on the requested map / solution choice:
          - direct maps      : solution = 0, show D
          - reflected maps   : solution = 1, show R

        Fallback rules:
          - D solution: shown unless the ray-traced samples do not yield a valid
            crossing on R_map. In that case, the same branch is re-traced outward
            from the station using AraSim, and the crossing is determined from
            those traced samples. Falls back to straight-line (SL) only if no valid
            RT crossing survives, or if the mapped RT elevation is unphysical for
            plotting.

          - R solution: shown if findPaths returns at least two solutions and the
            requested reflected / refracted branch yields a valid RT crossing on
            R_map. Falls back to SL if:
              (a) only one solution exists, so the R branch is absent,
              (b) neither the source->station nor the station->outward traced
                  samples yield a valid crossing on R_map, or
              (c) the mapped RT elevation is unphysical for plotting.

        Parameters
        ----------
        list_of_landmarks : list of str, optional
            Landmark names to include. Supported values:
            'ICL', 'IC22S', 'SPT', 'IC1S', 'SPIce', 'WT', 'SPRESSO'.
            Pass ['all'] to include all of the above.
            Default is ['IC22S', 'ICL', 'SPRESSO'].
        R_map : float
            Reconstruction sphere radius corresponding to the map, in meters.
            Must be provided and must be positive.
        list_of_cal_pulser_indices : list of int, optional
            Calpulser indices to include, e.g. [0, 1, 2, 3].
            Pass ['all'] for all four. Default is [1, 3].
        spice_depth : float, optional
            Depth of the SPIce pulser, in meters. Required if 'SPIce' is included.
        solution : int
            Requested ray solution:
              - 0 : direct path (shortest propagation time)
              - 1 : reflected path (second-shortest propagation time)

        Returns
        -------
        collect : dict
            Keys are landmark-name strings. Values are [r, elevation, azimuth],
            where r is the straight-line distance in meters, and elevation and
            azimuth are the ray-traced (or SL fallback) angles in degrees.
            Also includes:
              - 'critical_angle_rt' : ray-traced critical elevation angle mapped
                                      onto R_map, in degrees
              - '_marker_status'    : plotting status for each landmark
        """
        if R_map is None:
            raise ValueError("R_map must be provided to get_known_landmarks().")
        if R_map <= 0:
            raise ValueError(f"R_map must be positive, got {R_map}.")

        if list_of_landmarks is None:
            list_of_landmarks = ["IC22S", "ICL", "SPRESSO"]
        elif list_of_landmarks == ["all"]:
            list_of_landmarks = ["ICL", "IC22S", "IC1S", "WT", "SPRESSO"] #"SPT" removed due to depth of SPT bug in AraRoot 

        if list_of_cal_pulser_indices is None:
            list_of_cal_pulser_indices = [1, 3]
        elif list_of_cal_pulser_indices == ["all"]:
            list_of_cal_pulser_indices = [0, 1, 2, 3]

        if spice_depth is not None and "SPIce" not in list_of_landmarks:
            list_of_landmarks.append("SPIce")

        collect = {}
        marker_status = {}
        station_center = self.get_station_center()

        def _resolve(name, source_xyz, station_center, R_sphere, solution):
            """
            Resolve one landmark to its apparent direction on the reconstruction sphere.

            Computes the straight-line (SL) direction first; this is always valid and
            serves as the fallback. Then attempts ray tracing via
            get_raytraced_sphere_crossing_angles(). Falls back to SL with diagnostic
            output if ray tracing fails for any of the reasons documented above.

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
                0 = direct path, 1 = reflected path.

            Returns
            -------
            list
                [r, elevation, azimuth], where r is the straight-line distance in
                meters, and elevation / azimuth are in degrees. Ray-traced if
                successful, straight-line fallback otherwise.
            """
            r_sl, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(
                station_center, source_xyz
            )

            rt_sols = self.get_raytraced_sphere_crossing_angles(
                source_xyz=source_xyz,
                R_sphere=R_sphere,
                solution=solution,
                station_center=station_center,
            )

            if not rt_sols:
                marker_status[name] = "sl_fallback"
                if hasattr(self, "_last_rt_fail_reason") and self._last_rt_fail_reason:
                    print(f"[RAYTRACE FAIL] {name}: {self._last_rt_fail_reason}")
                return [r_sl, elev_sl, az_sl]

            marker_status[name] = (
                "sl_fallback" if rt_sols[0].get("used_sl_fallback", False) else "raytraced"
            )

            return [r_sl, rt_sols[0]["elevation"], rt_sols[0]["azimuth"]]

        # Local calpulsers
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
        for known_loc in ["ICL", "WT", "SPRESSO"]: #"SPT" removed because of unknown depth, AraRoot SPT depth entry is incorrect 
            if known_loc in list_of_landmarks:
                this_loc = np.array(self.get_southpole_landmarks(known_loc))[0]
                st_centric = self.get_global_to_station_centric(this_loc[:3])
                #if known_loc == "SPT":
                #    st_centric[2] = 7.55984245e-02  # SPT is removed from surface landmark list due to depth bug in AraRoot
                collect[known_loc] = _resolve(known_loc, st_centric, station_center, R_map, solution)

        # Ray-traced critical angle
        collect["critical_angle_rt"] = self.get_raytraced_critical_angle(R_map)
        collect["_marker_status"] = marker_status

        return collect
