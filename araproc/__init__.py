import logging
import ROOT
import os

try:
    ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR')+"/lib/libAraEvent.so")
except:
    logging.critical("Loading libAraEvent.so failed. Stop all work!")
    raise
