import os
import warnings

# Get the directory of ray tracing tables from an environment variable,
#  with a fallback or error if not set
ray_trace_tables_dir = os.getenv('RAY_TRACE_TABLES', None)

# If the environment variable is not set, use a default
if not ray_trace_tables_dir:
    ray_trace_tables_dir = "/data/ana/ARA/processing/support/raytrace_timing_tables/"
    warnings.warn(
        "The 'RAY_TRACE_TABLES' environment variable is not set and will use default: "
        + ray_trace_tables_dir
    )
# If the environment variable is set, make sure it's an absolute path
else: 
    ray_trace_tables_dir = os.path.abspath(ray_trace_tables_dir)

# Check if the directory exists
if not os.path.isdir(ray_trace_tables_dir):
    raise FileNotFoundError(
        f"The directory '{ray_trace_tables_dir}' does not exist. Please ensure the path is correct."
    )