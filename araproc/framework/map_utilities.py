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



    def _circle_intersect(self, p0, v, radius, t_max=1.0):
        """
        Intersect the line segment p(t) = p0 + t*v with the circle x^2 + z^2 = radius^2
        in the local (rho, dz) plane, for t in [0, t_max].

        Parameters
        ----------
        p0 : (float, float)
            (x0, z0) = (rho0, dz0)
        v : (float, float)
            (vx, vz) direction from p0 to p1 (segment vector)
        radius : float
            Circle radius
        t_max : float
            Maximum allowed t (use 1.0 for segment endpoint)

        Returns
        -------
        (x, z, t) or None
            Intersection point and interpolation parameter t.
        """
        # Unpack the starting point (rho0, dz0) and segment direction (drho, ddz).
        x0, z0 = float(p0[0]), float(p0[1])
        vx, vz = float(v[0]), float(v[1])
        R = float(radius)

        # Substitute x(t)=x0+t*vx and z(t)=z0+t*vz into x^2+z^2=R^2:
        # (vx^2+vz^2) t^2 + 2(x0*vx+z0*vz) t + (x0^2+z0^2-R^2) = 0
        a = vx * vx + vz * vz
        if abs(a) < 1e-30: # Degenerate segment (zero length) => no unique intersection
            return None

        b = 2.0 * (x0 * vx + z0 * vz)
        c = x0 * x0 + z0 * z0 - R * R

        disc = b * b - 4.0 * a * c # find Discriminant
        if disc < 0.0: # Discriminant < 0 means no real intersection with the circle
            return None

        # Compute the two quadratic roots (potential entry/exit points)
        s = np.sqrt(disc)

        # Two solutions; choose the one that lies on the segment and is closest to p0 (smallest t >= 0)

        t_candidates = [(-b - s) / (2.0 * a), (-b + s) / (2.0 * a)]

        # Keep only roots that lie within the segment parameter range.
        t_valid = [t for t in t_candidates if 0.0 <= t <= float(t_max)]

        if not t_valid: # Intersections exist on the infinite line, but not on this segment portion
            return None

        t = min(t_valid) # Choose the first intersection along the segment (closest to p0)
        # Evaluate the intersection point p(t) = p0 + t*v
        x = x0 + t * vx
        z = z0 + t * vz
        return (float(x), float(z), float(t))


    def _find_sphere_crossing_from_samples(
        self,
        xs,
        zs,
        z_ref,
        radius,
        reverse_search=False,
        ):
        
        """
        Find the intersection of sampled ray-trace points with the reconstruction
        sphere (circle in the 2D propagation plane) centered at the station.

        Geometry
        --------
        Using local coordinates:
            rho = xs
            dz  = zs - z_ref

        The reconstruction circle is:
            rho^2 + dz^2 = radius^2
        i.e.
            rho^2 + (z - z_ref)^2 = radius_map^2   

        Parameters
        ----------
        xs : array-like
            Sampled horizontal coordinates rho (meters), station-centered in the propagation plane.
            For backward traces from the station: xs starts near 0 and increases outward.
            For forward traces (source to station path): station-centered rho is passed, 
            e.g. rho = r_src - x.
        zs : array-like
            Sampled absolute z coordinates (meters) in station-centric coordinates.
        z_ref : float
            Station z (meters). Used to compute dz = z - z_ref.
        radius : float
            Reconstruction Map radius (meters).
        reverse_search : bool, default False
            If True, pick the crossing closest to the far/source side (useful for forward traces).
            If False, pick the first crossing encountered moving outward from the station.

        Returns
        -------
        inter : tuple or None
            (rho_cross, dz_cross, t) where t in [0,1] is the interpolation parameter along the
            chosen segment. dz_cross is relative to z_ref.
        i_cross : int or None
            Segment index i for the segment (i-1 -> i) containing the intersection.
        """
        xp = np.asarray(xs, dtype=float)
        dz = np.asarray(zs, dtype=float) - float(z_ref)

        if len(xp) < 2:
            return None, None

        d = np.hypot(xp, dz)

        # Find ALL candidate crossing segments by sign change 
        crossing_indices = [
            i for i in range(1, len(d))
            if (d[i - 1] - float(radius)) * (d[i] - float(radius)) <= 0.0
        ]
        if not crossing_indices:
            return None, None

        i_cross = crossing_indices[-1] if reverse_search else crossing_indices[0]

        # Segment endpoints in local (rho, dz)
        p0 = (float(xp[i_cross - 1]), float(dz[i_cross - 1])) # outside
        p1 = (float(xp[i_cross]),     float(dz[i_cross]))     # inside


        inter = self._circle_intersect(
            p0,
            (p1[0] - p0[0], p1[1] - p0[1]),
            float(radius),
            t_max=1.0,
        )

        # If interpolation failed on the chosen segment, try other candidates
        if inter is None:
            for i in crossing_indices:
                p0 = (float(xp[i - 1]), float(dz[i - 1])) # inside
                p1 = (float(xp[i]),     float(dz[i]))     # outside
                inter = self._circle_intersect(
                    p0,
                    (p1[0] - p0[0], p1[1] - p0[1]),
                    float(radius),
                    t_max=1.0,
                )
                if inter is not None:
                    return inter, i
            return None, None

        return inter, i_cross

    def _find_sphere_crossing_exact(self, xs, zs, z_rec, R):
        """
        Exact segment-wise intersection with circle:
            x^2 + (z - z_rec)^2 = R^2

        Returns
        -------
        inter : tuple or None
            (x_cross, dz_cross)
        seg_idx : int or None
            segment index i for (i-1 -> i)
        """


        """
        Compute the intersection of a sampled ray path with the reconstruction circle
        (sphere cross-section) by solving the quadratic exactly on each sample segment.

        Geometry / units
        ---------------
        - xs is station-centered horizontal distance (rho) in meters.
        - zs is absolute z in meters (station-centric coordinates).
        - z_rec is the station center (receiver) z (circle center) in meters
        - R is the Map radius in meters. 

        The circle condition is:
            rho^2 + (z - z_rec)^2 = R^2
        Define dz = (z - z_rec), then:
            rho^2 + dz^2 = R^2

        Method:
        ---------
        For each segment between consecutive samples (i-1 -> i), the parametric form is:
            rho(t) = rho0 + t*(rho1 - rho0)
            dz(t)  = dz0  + t*(dz1  - dz0)
        for t in [0, 1]. We solve:  rho(t)^2 + dz(t)^2 = R^2
        and select the first valid intersection encountered in the provided ordering of samples.

        Parameters
        ----------
        xs, zs : array-like
            Sampled points along the ray in the 2D propagation plane:
              - xs: rho (m)
              - zs: absolute z (m)
            Typically ordered from station outward for backward traces.
        z_rec : float
            Station z coordinate (m), used to compute dz = z - z_rec.
        R : float
            Reconstruction radius (m).

        Returns
        -------
        inter : tuple or None
            (rho_cross, dz_cross) where:
              - rho_cross is in meters
              - dz_cross = z_cross - z_rec is in meters
            Returns None if no segment intersects the circle.
        seg_idx : int or None
            Segment index i such that the crossing lies on the segment from (i-1) to i.
        """
        xs = np.asarray(xs, dtype=float)
        zs = np.asarray(zs, dtype=float)
        if xs.size < 2:
            return None, None

        dzs = zs - float(z_rec) # Convert absolute z to station-relative dz for the circle equation.
        R = float(R)

        # Look segment-by-segment and solve for an intersection on each segment.
        for i in range(1, len(xs)):
            x0, dz0 = float(xs[i - 1]), float(dzs[i - 1])
            x1, dz1 = float(xs[i]),     float(dzs[i])

            # Segment direction vector (drho, ddz)
            dx = x1 - x0
            ddz = dz1 - dz0

            # Plug rho(t)=x0+t*dx and dz(t)=dz0+t*ddz into rho(t)^2 + dz(t)^2 = R^2:
            # -> a t^2 + b t + c = 0
            a = dx * dx + ddz * ddz
            if abs(a) < 1e-30: # Degenerate segment (two identical points); skip.
                continue

            b = 2.0 * (x0 * dx + dz0 * ddz)
            c = x0 * x0 + dz0 * dz0 - R * R

            disc = b * b - 4.0 * a * c # Discriminant determines whether the segment line intersects the circle.
            if disc < 0.0: # No real intersection for this segment
                continue

            # Two candidate roots (enter/exit); keep only those lying on the segment t in [0,1]
            s = np.sqrt(disc)
            t_candidates = [(-b - s) / (2.0 * a), (-b + s) / (2.0 * a)]
            t_valid = [t for t in t_candidates if 0.0 <= t <= 1.0]
            if not t_valid:
                continue

            t = min(t_valid) # Choose the earliest intersection along this segment.
            # Convert back to intersection point in (rho, dz).
            x = x0 + t * dx
            dz = dz0 + t * ddz
            return (float(x), float(dz)), int(i)

        # No segment intersected the circle.
        return None, None


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
        Ray-trace a source to station path using AraSim and return the apparent direction
        of the ray at its intersection with the reconstruction sphere of radius R_sphere.

        Coordinate conventions / units
        ------------------------------
        All distances are in meters and angles returned are in degrees.

        Inputs:
        - source_xyz: (x, y, z) in station-centric Cartesian coordinates (m).
          z is the absolute station-centric depth coordinate.
        - station_center: (x, y, z) in station-centric cartesian coordinates (m).
          If None, uses self.get_station_center().
        - R_sphere: reconstruction sphere radius (m), centered at the station.
 
        - rho is the station-centered horizontal distance in that plane (m), rho >= 0
        - dz = z - z_station (m)

        The reconstruction sphere in the (rho, z) plane is:
            rho^2 + (z - z_station)^2 = R_sphere^2
        i.e. rho^2 + dz^2 = R_sphere^2.

        How this works ?
        ---------------------------------
        1) Use AraSim TraceFinder::findPaths to obtain ray solutions between source and station.
           Solutions are ordered by time-of-flight; Hence:
             - solution=0 -> 'D' (shortest TOF, "direct")
             - solution=1 -> 'R' (second TOF, "refracted/alternate")
        2) Trace the chosen solution with rt_doTrace_savePoints to obtain sampled points.
        3) Determine where the sampled ray intersects the reconstruction circle:
           - Try a forward trace (source to station) converted to station-centered rho.
             For forward traces we select the "outer" crossing (reverse_search=True).
           - If no crossing is found (can happen when the source lies inside the Radius map of interest),
             do a backward trace outward from the station using a folded receiptAngle and
             take the first outward crossing (reverse_search=False).

        Returned direction definition
        -----------------------------
        The returned elevation is computed at the sphere-crossing point as:
            elevation = atan2(dz_cross, rho_cross)  [degrees]
        where dz_cross = (z_cross - z_station) and rho_cross is the station-centered horizontal
        distance at the crossing.

        Parameters
        ----------
        source_xyz : array-like
            Source coordinates [x, y, z] (m), station-centric.
        R_sphere : float
            Reconstruction sphere radius (m), must be > 0.
        accuracy : float, default 0.001
            Convergence tolerance passed to TraceFinder::findPaths.
        solution : int, default 0
            0 selects the shortest TOF solution ('D'); 1 selects the second ('R').
        frequency : float, default 0.5
            Frequency passed to the ray tracer (AraSim units; typically GHz in this interface).
        polarization : float, default 0.0
            Polarization argument passed to the ray tracer.
        station_center : array-like or None
            Station center [x, y, z] (m). If None, computed from geometry.

        Returns
        -------
        list of dict
            A list with 0 or 1 entries. On success, one dict is returned:
              - 'label' : 'D' or 'R'
              - 'elevation' : float (deg), elevation at the sphere crossing
              - 'azimuth' : float (deg), azimuth of source direction
              - 'tof' : float (ns), time of flight from AraSim pathTime
              - 'path_len' : float, path length from AraSim pathLen
              - 'allowedUsed' : int, reflection mode used
              - 'used_sl_fallback' : bool, True only if an unphysical angle forced SL fallback

            If no valid ray solution/crossing is found, returns [] and sets
            self._last_rt_fail_reason.
        """

        if R_sphere is None:
            raise ValueError("R_sphere must be provided to get_raytraced_sphere_crossing_angles().")
        if R_sphere <= 0:
            raise ValueError(f"R_sphere must be positive, got {R_sphere}.")

        # Source/receiver coordinates in station-centric cartesian frame (meters)

        src = np.asarray(source_xyz, dtype=float)
        rec = np.asarray(
            station_center if station_center is not None else self.get_station_center(),
            dtype=float,
        )
        # Extract sphere radius and absolute z (meters)
        R = float(R_sphere)
        z_src = float(src[2])
        z_rec = float(rec[2])
        # --- Convert 3D source position to 2D propagation-plane quantities ---
        # r_src: station-centered horizontal distance to the source (meters)
        # az_phi: azimuth angle of the source in station coordinates (radians)

        dx = float(src[0] - rec[0])
        dy = float(src[1] - rec[1])
        r_src = float(np.hypot(dx, dy))
        az_phi = float(np.arctan2(dy, dx))
        # --- Reset last-raytrace diagnostics (for debug) ---
        self._last_rt_npaths = 0
        self._last_rt_paths = []
        self._last_rt_fail_reason = None
        # --- Configure AraSim vectors in the ray-tracer propagation plane ---
        ROOT._rt_src.SetXYZ(r_src, 0.0, z_src)
        ROOT._rt_rec.SetXYZ(0.0, 0.0, z_rec)

        
        sol_cnt = ctypes.c_int()
        sol_err = ctypes.c_int()
        allowed_used = int(ROOT.RayTrace.SurfaceReflection)

        paths = ROOT._rt_tf.findPaths(
            ROOT._rt_src,
            ROOT._rt_rec,
            float(frequency),
            ROOT.TMath.Pi() / 2,
            sol_cnt,
            sol_err,
            allowed_used,
            float(accuracy),
        )
        # Handle "no solutions" or raytracer errors 
        if sol_err.value != 0 or not paths or len(paths) == 0:
            self._last_rt_fail_reason = (
                f"findPaths error (serr={sol_err.value})"
                if sol_err.value != 0
                else "findPaths returned zero solutions"
            )
            return []
        #Sort solutions by time-of-flight (TOF) and label them by order
        sorted_paths = sorted(paths, key=lambda p: float(p.pathTime))

        labeled = []
        if len(sorted_paths) >= 1:
            labeled.append((sorted_paths[0], "D")) # shortest TOF
        if len(sorted_paths) >= 2:
            labeled.append((sorted_paths[1], "R"))  # second-shortest TOF

        #Select which branch to trace based on 'solution' argument
        if solution == 0:
            to_trace = [(p, l) for p, l in labeled if l == "D"]
        elif solution == 1:
            to_trace = [(p, l) for p, l in labeled if l == "R"]
            if not to_trace:
                self._last_rt_fail_reason = "R solution not found: only 1 ray-trace solution exists"
                return []
        else:
            raise ValueError(f"solution must be 0 or 1, got {solution}.")

        # Straight-line (SL) direction computed once; used only as a fallback
        _, elev_sl, az_sl = self.get_relative_cartesian_to_spherical(rec, src)

        out = []
        for path, label in to_trace: 
            used_sl_fallback = False
            #Forward trace (source to station) and save all intermediate ray traced points
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
            #Retrieve the saved ray samples from the C++ helper vectors.
            xs = np.asarray(list(ROOT.rt_get_x()), dtype=float) # rho-like coordinate along the trace
            zs = np.asarray(list(ROOT.rt_get_z()), dtype=float) # absolute z samples

            if sol_error.value != 0 or len(xs) < 2:  # If the forward trace failed, skip this branch.
                self._last_rt_fail_reason = (
                    f"forward trace failed for label={label}: sol_error={sol_error.value}, npts={len(xs)}"
                )
                continue

            # Convert forward-trace rho coordinate into station-centered rho.
            # In this forward trace, x increases from source to station, so station-centered rho is (r_src - x).
            xp_forward = r_src - xs

            #Find sphere crossing from the forward samples
            # reverse_search=True chooses the crossing closer to the station side for forward-ordered points.
            inter, _ = self._find_sphere_crossing_from_samples(
                xp_forward,
                zs,
                z_rec,
                R,
                reverse_search=True,
            )

            # If no forward crossing: backward trace outward from station
            if inter is None:
                # Need receiptAngle to trace back outward along the same branch.
                if not hasattr(path, "receiptAngle"):
                    self._last_rt_fail_reason = (
                        f"no forward crossing for label={label} and receiptAngle unavailable for backward tracing"
                    )
                    continue

                # Fold the receipt angle into [0, pi/2] for a valid outward-launch angle.
                theta_receipt = float(path.receiptAngle)
                theta_back = theta_receipt if theta_receipt <= np.pi / 2 else (np.pi - theta_receipt)

                # Backward tracing distance: ensure it's far enough to cross the R_map circle

                # Choose a sufficiently distant backward-trace target so the station to outward ray is long enough
                # (and has enough sampled points) to reliably cross the reconstruction circle of radius R.
                # - 3*R: is choden to ensure we trace several sphere radii outward
                # - r_src+50 m: ensures we trace at least beyond the source horizontal range (plus margin) when R is small
                # - 100 m: absolute minimum to avoid too-short traces for very small R

                r_back_target = max(3.0 * R, r_src + 50.0, 100.0)
                z_back_target = z_rec

                
                # Backward trace (station to outward) and save all intermediate points.
                sol_error_back = ctypes.c_int(0)
                ROOT.rt_doTrace_savePoints(
                    z_rec,
                    float(theta_back),
                    ROOT.RayTrace.rayTargetRecord(float(z_back_target), float(r_back_target)),
                    allowed_used,
                    float(frequency),
                    float(polarization),
                    sol_error_back,
                )

                xs_back = np.asarray(list(ROOT.rt_get_x()), dtype=float)
                zs_back = np.asarray(list(ROOT.rt_get_z()), dtype=float)


                 # If backward trace failed, skip this branch
                if sol_error_back.value != 0 or len(xs_back) < 2:
                    self._last_rt_fail_reason = (
                        f"backward trace failed for label={label}: sol_error={sol_error_back.value}, npts={len(xs_back)}"
                    )
                    continue
                # Find the first outward crossing on the backward-trace samples.
                inter, _ = self._find_sphere_crossing_from_samples(
                    xs_back,
                    zs_back,
                    z_rec,
                    R,
                    reverse_search=False,
                )
                # If still no crossing, give up on this RT solution
                #record the failure reason and print it
                # so we can see why we fell back / skipped this branch.
                if inter is None:
                    self._last_rt_fail_reason = (
                        f"ray-traced solution found for label={label}, but neither forward nor backward samples crossed "
                        f"R_sphere={R:.3f} m"
                    )
                    continue
            # convert crossing point (rho_cross, dz_cross) into elevation-azimuth
            xp_c, dz_c, _ = inter
            elev = float(np.degrees(np.arctan2(dz_c, xp_c))) # elevation at sphere crossing (deg)
            az = float(np.degrees(az_phi))   # azimuth defined by 3D source position (deg)

            # Sanity check: if the mapped RT elevation is unphysical (outside [-90, 90]),
            # fall back to the straight-line direction. This should not happen, but keeping it as a safeguard.
            if not (-90.0 <= elev <= 90.0):
                elev, az = float(elev_sl), float(az_sl)
                used_sl_fallback = True
                self._last_rt_fail_reason = (
                    f"ray-traced elevation outside [-90,90]; used SL fallback for label={label}"
                )

            # Record output dictionary for this chosen solution branch
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

        # If everything failed, set a final diagnostic reason.
        if not out and not self._last_rt_fail_reason:
            self._last_rt_fail_reason = f"no valid crossing to R_sphere={R:.3f} m survived"

        return out


    def get_raytraced_critical_angle(self, R_sphere, frequency=0.5, polarization=0.0):

        """
        Compute the apparent critical-angle direction on the reconstruction sphere
        by ray tracing a near-critical ray launched from the station.

        Method : 
        ------------------------
        The function:
          1) computes the local geometric critical elevation at the station using
             the refractive indices (air vs ice),
          2) launches a ray from the station at that near-critical direction using AraSim,
          3) determines where that ray intersects the reconstruction sphere 
             of radius R_sphere centered at the station,
          4) returns the elevation angle at that sphere-crossing point.
          5) For "small" spheres (R_sphere less than the vertical distance to the turning
          point), the ray intersects the sphere on the upward leg; this is handled via
          exact segment-wise intersection in `_find_sphere_crossing_exact`.

        Coordinate conventions / units
        ------------------------------
        - All distances are in meters.
        - Returned angle is degrees.
        - z_rec = Station Center depth (z) coordinate (m).
        - rho = horizontal distance from station (m), rho >= 0
        - dz  = z - z_rec (m)
        The reconstruction circle is:  rho^2 + dz^2 = R_sphere^2

        Returned value
        --------------
        The returned elevation is:  elevation = atan2(dz_cross, rho_cross)   [degrees]
        where (rho_cross, dz_cross) is the intersection of the traced ray with the
        reconstruction circle.

        Parameters
        ----------
        R_sphere : float
            Reconstruction sphere radius (m), must be > 0.
        frequency : float, default 0.5
            Frequency passed to AraSim ray tracer (AraSim interface units; typically GHz).
        polarization : float, default 0.0
            Polarization argument passed to the ray tracer.

        Notes
        -----
        - This function returns a ray-traced critical-angle elevation.
        - If the AraSim trace fails (nonzero solver error / too few samples) or a turning point
          cannot be identified, the function falls back to the geometric critical elevation at
          the station (from `get_critical_angle()`). 
        - We always expect the ray-traced critical angle result; the fallback is retained as 
          a safeguard for unusual inputs.

        """


        if R_sphere is None or R_sphere <= 0:
            raise ValueError(f"R_sphere must be positive, got {R_sphere}")

        z_rec = float(self.get_station_center()[2])
        R = float(R_sphere)

        crit_elev_deg = float(self.get_critical_angle())

        # Convert elevation to AraSim launch angle convention (theta from +z axis, radians):
        # elevation = 90 - theta; so,  theta = 90 - elevation
        theta_launch  = np.deg2rad(90.0 - crit_elev_deg)

        # Choose a far-enough target so the traced path develops its turning point ---
        # Trace outward well beyond the sphere (10*R) and also enforce a minimum (2000 m)
        # to ensure robust sampling even for small R.

        r_target = max(10.0 * R, 2000.0)
        z_target = z_rec 

        #Run AraSim tracer from station outward and save all sampled points 
        sol_error = ctypes.c_int(0)
        ROOT.rt_doTrace_savePoints(
            z_rec,
            float(theta_launch), # launch angle (rad)
            ROOT.RayTrace.rayTargetRecord(z_target, r_target), # outward target (z, r)
            int(ROOT.RayTrace.SurfaceReflection),
            float(frequency),
            float(polarization),
            sol_error,
        )

        # Retrieve sampled ray points from the AraSim helper buffers.
        xs = np.asarray(list(ROOT.rt_get_x()), dtype=float)
        zs = np.asarray(list(ROOT.rt_get_z()), dtype=float)

        # Fallback: if trace failed or produced too few points, return geometric critical angle
        # this has not been observed, but we keep it as a safeguard
        if sol_error.value != 0 or len(xs) < 2:
            print("[CRIT] Trace failed -> geometric fallback")
            return crit_elev_deg

        # Identify the turning point (maximum z along the trajectory)
        # Turning point is detected by slope sign change: dz/ds goes from + to -.
        turning_idx = None
        for i in range(1, len(xs) - 1):
            if (zs[i] - zs[i-1]) >= 0 and (zs[i+1] - zs[i]) <= 0:
                turning_idx = i
                break

        # If no turning point is found, we cannot define the critical ray reliably fallback.
        # this has not been observed, but we keep it as a safeguard
        if turning_idx is None:
            print("[CRIT] No turning point found -> geometric fallback")
            return crit_elev_deg
        # Compute vertical offset of turning point relative to station 
        z_turn = zs[turning_idx]
        dz_h   = z_turn - z_rec # dz at the turning point (m)

        # Determine intersection with reconstruction circle rho^2 + dz^2 = R^2
        if abs(dz_h) < R:
            # If the turning point lies inside the circle vertically, the circle is intersected
            # on the upward-going leg at dz=dz_h; solve rho from rho = sqrt(R^2 - dz^2).
            xp_c  = np.sqrt(R**2 - dz_h**2)
            dz_c  = dz_h
        else:
            # Otherwise, find the exact segment-wise intersection using samples up to the turning point.
            inter, _ = self._find_sphere_crossing_exact(
                xs[:turning_idx+1],
                zs[:turning_idx+1],
                z_rec,
                R,
            )
            # If no intersection is found (unlikely), fall back to geometric critical angle.
            if inter is None:
                return crit_elev_deg
            xp_c, dz_c = inter

        # Convert the crossing point (rho_cross, dz_cross) to elevation angle (degrees)
        elev = float(np.degrees(np.arctan2(dz_c, xp_c)))
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
        Compute apparent directions of known landmarks on a reconstruction sphere of chosen radius R_map

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
            crossing on R_map (source is inside R_map). 
            In that case, the same branch is re-traced outward
            from the station using AraSim, and the crossing is determined from
            those traced samples. Falls back to straight-line (SL) only if no valid
            RT crossing survives, or if the mapped RT elevation is unphysical for
            plotting. We always find a D solution for all landmarks in all stations. 

          - R solution: shown if findPaths returns at least two solutions and the
            requested reflected / refracted branch yields a valid RT crossing on
            R_map. Falls back to SL if only one solution exists, so the R branch is absent (that is the 
            case for surface landmarks where we get SL fallback since AraSim returned only 1 solution 
            for surface landmarks)
              

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
            list_of_landmarks = ["ICL", "IC22S", "IC1S", "WT", "SPRESSO", "SPT"]
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
                Landmark label used for the collect key.
            source_xyz : array-like
                Source position [x, y, z] in station-centric Cartesian coordinates,
                in meters.
            station_center : array-like
                Station center [x, y, z] in station-centric Cartesian coordinates,
                in meters.
            R_sphere : float
                Reconstruction sphere radius, in meters.
            solution : int
                0 = direct path, 1 = second-shortest TOF path ('R')

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
        
        for known_loc in ["ICL", "WT", "SPRESSO", "SPT"]:
            if known_loc in list_of_landmarks:
                this_loc = np.array(self.get_southpole_landmarks(known_loc))[0]
                st_centric = self.get_global_to_station_centric(this_loc)
                collect[known_loc] = _resolve(known_loc, st_centric, station_center, R_map, solution)

        # Ray-traced critical angle
        
        critical_angle_rt = self.get_raytraced_critical_angle(R_map)
        collect["critical_angle_rt"] = critical_angle_rt
        collect["_marker_status"] = marker_status

        return collect




