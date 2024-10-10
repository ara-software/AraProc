import logging
import ROOT
import os

# libAraEvent
try:
    success_flag_libAraEvent = ROOT.gSystem.Load(os.environ.get("ARA_UTIL_INSTALL_DIR")+"/lib/libAraEvent.so")
except: # noqa: E722
    logging.critical("Loading libAraEvent.so failed. Stop all work!")
    raise ImportError("Failed to load libAraEvent.so")
if success_flag_libAraEvent not in [0, 1]:
        raise ImportError("Failed to load libAraEvent.so")

# libAra (this is AraSim)
try:
    success_flag_libAra = ROOT.gSystem.Load(os.environ.get("ARA_UTIL_INSTALL_DIR")+"/lib/libAra.so")
except: # noqa: E722
    logging.critical("Loading libAra.so failed. Stop all work!")
    raise ImportError("Failed to load libAra.so")
if success_flag_libAra not in [0, 1]:
        raise ImportError("Failed to load libAra.so")


# libAraCorrelator
try:
    success_flag_libAraCorrelator = ROOT.gSystem.Load(os.environ.get("ARA_UTIL_INSTALL_DIR")+"/lib/libAraCorrelator.so")
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

success_flag_FFTinclude = ROOT.gInterpreter.Declare(f'#include "s{os.environ.get("ARA_DEPS_INSTALL_DIR")}/include/FFTtools.h"')
if not success_flag_FFTinclude:
    raise ImportError("Including FFTtools.h failed. Stop all work!")
