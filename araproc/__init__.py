import logging
import ROOT
import os

try:
    ROOT.gSystem.Load(os.environ.get("ARA_UTIL_INSTALL_DIR")+"/lib/libAraEvent.so")
    ROOT.gSystem.Load(os.environ.get("ARA_UTIL_INSTALL_DIR")+"/lib/libAraCorrelator.so")
    ROOT.gSystem.Load(os.environ.get("ARA_DEPS_INSTALL_DIR")+"/lib/libRootFftwWrapper.so")
    ROOT.gInterpreter.ProcessLine(f'#include "{os.environ.get("ARA_DEPS_INSTALL_DIR")}/include/FFTtools.h"')
except:
    logging.critical("Loading libAraEvent.so failed. Stop all work!")
    raise
