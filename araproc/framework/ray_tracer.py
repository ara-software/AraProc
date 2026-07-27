import os
import threading
import ctypes
import math

import ROOT


_INITIALIZATION_LOCK = threading.Lock()
_RAYTRACE_MODEL_LOADED = False


def load_raytrace_model(arasimsrc=None):
    """
    Load and initialize the AraSim ray-tracing model in ROOT's interpreter.

    The refractive-index model is:

        n(z) = A - B * exp(-C * abs(z))

    with A=1.780, B=0.454, and C=0.0202 /m.

    Parameters
    ----------
    arasimsrc : str or os.PathLike or None
        AraSim source directory containing ``RayTrace.h``,
        ``RayTrace_IceModels.h``, and ``Vector.h``. If None, use the
        ``ARA_SIM_DIR`` environment variable.

    Raises
    ------
    RuntimeError
        If no AraSim source directory is configured, a required header is
        missing, or ROOT fails to declare the helper functions.
    """
    global _RAYTRACE_MODEL_LOADED

    if _RAYTRACE_MODEL_LOADED:
        return

    with _INITIALIZATION_LOCK:
        if _RAYTRACE_MODEL_LOADED:
            return

        if arasimsrc is None:
            arasimsrc = os.getenv("ARA_SIM_DIR")

        if not arasimsrc:
            raise RuntimeError(
                "ARA_SIM_DIR environment variable not set. "
                "Please set it to the AraSim source directory or pass "
                "arasimsrc explicitly to load_raytrace_model()."
            )

        arasimsrc = os.path.abspath(os.fspath(arasimsrc))

        required_headers = (
            "RayTrace.h",
            "RayTrace_IceModels.h",
            "Vector.h",
        )
        missing_headers = [
            name
            for name in required_headers
            if not os.path.isfile(os.path.join(arasimsrc, name))
        ]

        if missing_headers:
            raise RuntimeError(
                f"AraSim source directory {arasimsrc!r} is missing: "
                f"{', '.join(missing_headers)}"
            )

        # Load the AraSim ray-tracing headers into ROOT's interpreter so the
        # corresponding C++ classes can be used from Python.
        ROOT.gInterpreter.ProcessLine(
            f'#include "{os.path.join(arasimsrc, "RayTrace.h")}"'
        )
        ROOT.gInterpreter.ProcessLine(
            f'#include "{os.path.join(arasimsrc, "RayTrace_IceModels.h")}"'
        )
        ROOT.gInterpreter.ProcessLine(
            f'#include "{os.path.join(arasimsrc, "Vector.h")}"'
        )
        ROOT.gInterpreter.ProcessLine("#include <vector>")

        # Instantiate the attenuation model and the exponential refractive-index
        # model used by the ray tracer. The refractive-index parameters are the
        # same as those used in the PA analysis paper.
        ROOT.gInterpreter.ProcessLine(
            "auto _rt_atten = boost::shared_ptr<basicAttenuationModel>"
            "(new basicAttenuationModel);"
        )
        ROOT.gInterpreter.ProcessLine(
            "auto _rt_refr = boost::shared_ptr<exponentialRefractiveIndex>"
            "(new exponentialRefractiveIndex(1.3260,1.780,0.0202));"
        )

        # Create the AraSim ray-tracing engine from the refractive-index and
        # attenuation models above.
        ROOT.gInterpreter.ProcessLine(
            "RayTrace::TraceFinder _rt_tf(_rt_refr, _rt_atten);"
        )

        # Declare source/receiver vectors used by the ray-tracing machinery.
        ROOT.gInterpreter.ProcessLine("Vector _rt_src; Vector _rt_rec;")

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
        ROOT.gInterpreter.Declare(
            r"""
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
                void operator()(
                    const RayTrace::fullRayPosition& p,
                    RayTrace::RKStepType t
                ) const {
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

                return _rt_tf.doTrace<
                    RayTrace::fullRayPosition,
                    RTPointSaver
                >(
                    src_depth,
                    launch_theta,
                    target,
                    allowedReflections,
                    frequency,
                    polarization,
                    sol_error,
                    cb
                );
            }

            // Accessors used from Python to retrieve the saved traced path.
            const std::vector<double>& rt_get_x() {
                return _rt_x;
            }

            const std::vector<double>& rt_get_z() {
                return _rt_z;
            }

            const std::vector<int>& rt_get_stepType() {
                return _rt_stepType;
            }
            """
        )

        if not hasattr(ROOT, "rt_doTrace_savePoints"):
            raise RuntimeError(
                "Failed to declare rt_doTrace_savePoints in ROOT."
            )

        _RAYTRACE_MODEL_LOADED = True


def raytrace_model_is_loaded():
    """
    Return whether this module successfully initialized the ray tracer.

    Returns
    -------
    bool
        True after ``load_raytrace_model()`` completes successfully.
    """
    return _RAYTRACE_MODEL_LOADED


def get_trace_finder():
    """
    Return the initialized AraSim ``RayTrace::TraceFinder``.

    Returns
    -------
    ROOT.RayTrace.TraceFinder
        The ROOT proxy for the initialized trace finder.
    """
    load_raytrace_model()
    return ROOT._rt_tf


def get_source_vector():
    """
    Return the reusable ray-tracing source vector.
    """
    load_raytrace_model()
    return ROOT._rt_src


def get_receiver_vector():
    """
    Return the reusable ray-tracing receiver vector.
    """
    load_raytrace_model()
    return ROOT._rt_rec


def get_saved_trace_points():
    """
    Return copies of the points saved by the most recent trace.

    Returns
    -------
    tuple of list
        ``(x, z, step_type)``.
    """
    load_raytrace_model()

    return (
        list(ROOT.rt_get_x()),
        list(ROOT.rt_get_z()),
        list(ROOT.rt_get_stepType()),
    )

def find_ray_paths(
    source_xyz,
    receiver_xyz,
    frequency=0.5,
    polarization=None,
    accuracy=0.001,
    allowed_reflections=None,
):
    """
    Find AraSim ray-trace solutions between two Cartesian positions.

    The refractive-index model is horizontally uniform, so the 3D source and
    receiver positions are reduced to AraSim's 2D propagation plane.

    Parameters
    ----------
    source_xyz : array-like
        Source ``[x, y, z]`` in meters.
    receiver_xyz : array-like
        Receiver ``[x, y, z]`` in meters, in the same coordinate system as
        ``source_xyz``.
    frequency : float, default 0.5
        Frequency passed to AraSim.
    polarization : float or None, default None
        Polarization angle passed to ``TraceFinder.findPaths``. If None, use
        ``pi / 2``, matching the existing map-utilities implementation.
    accuracy : float, default 0.001
        AraSim path-finding accuracy.
    allowed_reflections : int or None, default None
        AraSim reflection mode. If None, use
        ``ROOT.RayTrace.SurfaceReflection``.

    Returns
    -------
    list
        AraSim ``TraceRecord`` objects sorted by increasing path time.
        The first entry is the direct solution and the second entry is the
        reflected/refracted solution, matching ``map_utilities.py``.

        An empty list is returned if AraSim reports an error or finds no path.
    """
    load_raytrace_model()

    if len(source_xyz) != 3:
        raise ValueError(
            f"source_xyz must contain exactly three coordinates; "
            f"got {source_xyz!r}."
        )

    if len(receiver_xyz) != 3:
        raise ValueError(
            f"receiver_xyz must contain exactly three coordinates; "
            f"got {receiver_xyz!r}."
        )

    source_x, source_y, source_z = map(float, source_xyz)
    receiver_x, receiver_y, receiver_z = map(float, receiver_xyz)

    horizontal_distance = math.hypot(
        source_x - receiver_x,
        source_y - receiver_y,
    )

    ROOT._rt_src.SetXYZ(horizontal_distance, 0.0, source_z)
    ROOT._rt_rec.SetXYZ(0.0, 0.0, receiver_z)

    solution_count = ctypes.c_int(0)
    solution_error = ctypes.c_int(0)

    if polarization is None:
        polarization = ROOT.TMath.Pi() / 2.0

    if allowed_reflections is None:
        allowed_reflections = int(ROOT.RayTrace.SurfaceReflection)

    paths = ROOT._rt_tf.findPaths(
        ROOT._rt_src,
        ROOT._rt_rec,
        float(frequency),
        float(polarization),
        solution_count,
        solution_error,
        int(allowed_reflections),
        float(accuracy),
    )

    if solution_error.value != 0 or not paths:
        return []

    return sorted(paths, key=lambda path: float(path.pathTime))

__all__ = [
    "find_ray_paths",
    "get_receiver_vector",
    "get_saved_trace_points",
    "get_source_vector",
    "get_trace_finder",
    "load_raytrace_model",
    "raytrace_model_is_loaded",
]
