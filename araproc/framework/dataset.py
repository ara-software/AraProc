import logging
import os
import ROOT

ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR')+"/lib/libAraEvent.so")

class AraDataset:

    def __init__(self, 
                 path_to_data_file : str,
                 path_to_pedestal_file : str,
                 ):
          
        self.path_to_root_file = None
        self.path_to_pedestal_file = None
        self._ttree = None
        self._root_tfile = None
        
        # check if they gave us strings
        if not isinstance(path_to_data_file, str):
            raise TypeError("Path to root data file must be a string")

        if not isinstance(path_to_pedestal_file, str):
            raise TypeError("Path to pedestal file must be a string")

        # check if the file exists
        if not os.path.exists(path_to_data_file):
            raise FileNotFoundError(f"Data root File ({path_to_data_file}) not found")

        if not os.path.exists(path_to_pedestal_file):
            raise FileNotFoundError(f"Pedestal File ({path_to_pedestal_file}) not found")

        # check they gave us files and not directories
        if not os.path.isfile(path_to_data_file):
            raise ValueError(f"{path_to_data_file} looks like a directory, not a file")

        # check they gave us files and not directories
        if not os.path.isfile(path_to_pedestal_file):
            raise ValueError(f"{path_to_pedestal_file} looks like a directory, not a file")
        
        self.path_to_data_file = path_to_data_file
        self.path_to_pedestal_file = path_to_pedestal_file
        self.open_tfile()

    def open_tfile(self):
        try:
            self._root_tfile = ROOT.TFile(self.path_to_data_file, "READ")
            logging.info(f"Successfully opened {self.path_to_data_file}")
        except:
            logging.critical(f"Opening {self.path_to_data_file} failed")
            raise



    