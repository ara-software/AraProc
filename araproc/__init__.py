import logging
import ROOT
import os
import subprocess

# set the git hash
try:
    path = os.path.abspath(os.path.dirname(__file__))

    # Walk up the directory tree to find .git
    while path != "/" and not os.path.isdir(os.path.join(path, ".git")):
        path = os.path.dirname(path)

    if os.path.isdir(os.path.join(path, ".git")):
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=path,
            stderr=subprocess.DEVNULL
        ).strip().decode("utf-8")
except Exception:
    raise ValueError("Could not determine git hash")

# libAraEvent
try:
    success_flag_libAraEvent = ROOT.gSystem.Load(os.environ.get("ARA_ROOT_LIB_DIR")+"/libAraEvent.so")
except: # noqa: E722
    logging.critical("Loading libAraEvent.so failed. Stop all work!")
    raise ImportError("Failed to load libAraEvent.so")
if success_flag_libAraEvent not in [0, 1]:
        raise ImportError("Failed to load libAraEvent.so")

# libAra (this is AraSim)
try:
    success_flag_libAra = ROOT.gSystem.Load(os.environ.get("ARA_SIM_LIB_DIR")+"/libAra.so")
except: # noqa: E722
    logging.critical("Loading libAra.so failed. Stop all work!")
    raise ImportError("Failed to load libAra.so")
if success_flag_libAra not in [0, 1]:
        raise ImportError("Failed to load libAra.so")


# libAraCorrelator
try:
    success_flag_libAraCorrelator = ROOT.gSystem.Load(os.environ.get("ARA_ROOT_LIB_DIR")+"/libAraCorrelator.so")
except: # noqa: E722
    logging.critical("Loading libAraCorrelator.so failed. Stop all work!")
    raise ImportError("Failed to load libAraCorrelator.so")
if success_flag_libAraCorrelator not in [0, 1]:
        raise ImportError("Failed to load libAraCorrelator.so")

# FFTtools
try:
    success_flag_libFFT = ROOT.gSystem.Load(os.environ.get("ARA_DEPS_INSTALL_DIR")+"/lib/libRootFftwWrapper.so")
except: # noqa: E722
    logging.critical("Loading libRootFftwWrapper.so failed. Stop all work!")
    raise ImportError("Failed to load libRootFftwWrapper.so")
if success_flag_libFFT not in [0, 1]:
    raise ImportError("Failed to load libRootFftwWrapper.so")

success_flag_FFTinclude = ROOT.gInterpreter.Declare(f'#include "{os.environ.get("ARA_DEPS_INSTALL_DIR")}/include/FFTtools.h"')
if not success_flag_FFTinclude:
    raise ImportError("Including FFTtools.h failed. Stop all work!")
