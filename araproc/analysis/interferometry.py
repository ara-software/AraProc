import numpy as np
import os
import ROOT


class RayTraceCorrelatorWrapper:

    """
    A wrapper class around the AraRoot RayTraceCorrelator

    This things has the responsibility to manage interferometry
    via the AraRoot RayTraceCorrelator.
    
    ...

    Attributes
    ----------
    station_id : int
       The station ID
    correlators : dictionary
        A dictionary.
        Keys = radii
        Values = the ray trace correlator instance

    """

    def __init__(self,
                 station_id : int,
                 num_channels : int =  16 # by default 16
                 ):
        if station_id not in [1, 2, 3, 4, 5]:
            raise KeyError(f"Station {station_id} is not supported")
        if not isinstance(num_channels, int):
            raise ValueError(f"{num_channels} needs to be an int")
        if (num_channels < 0) or not np.isfinite(num_channels):
            raise ValueError(f"{num_channels} is invalid (negative, not finite, etc.)")

        self.station_id = station_id
        self.num_channels = num_channels

        self.correlators = {}


    def __sanitize_inputs(self,
                           path_to_dir_file,
                           path_to_ref_file,
                           ):
        
        # check if they gave us strings
        if not isinstance(path_to_dir_file, str):
            raise TypeError("Path to dir file must be a string")
        if not isinstance(path_to_ref_file, str):
            raise TypeError("Path to ref file must be a string")

        # check if the file exists
        if not os.path.exists(path_to_dir_file):
            raise FileNotFoundError(f"Direct root File ({path_to_dir_file}) not found")
        if not os.path.exists(path_to_ref_file):
            raise FileNotFoundError(f"Ref root File ({path_to_ref_file}) not found")

        # check they gave us files and not directories
        if not os.path.isfile(path_to_dir_file):
            raise ValueError(f"{path_to_dir_file} looks like a directory, not a file")
        # check they gave us files and not directories
        if not os.path.isfile(path_to_ref_file):
            raise ValueError(f"{path_to_ref_file} looks like a directory, not a file")
        
        # now that we're sure these are clean, we can assign the variables
        self.path_to_data_file = path_to_dir_file
        self.path_to_pedestal_file = path_to_ref_file


    def add_rtc(self, 
                radius : float,
                path_to_dir_file : str,
                path_to_ref_file : str,
                ref_name = None, 
                angular_size : float  = 1., 
                ):

        """
        Parameters
        ----------
        ref_name : None
            The refernece name you want used as the key in the self.correlators dictionary.
            If you don't put anything in, the key will be set to the radius.
        radius : float
            The radius for this correlator in meters.
        path_to_dir_file : str
            the full path to the direct ray tracing file
        path_to_pedestal_file : str
            the full path to the reflected file
        angular_side : float
            Size of the angular grid in degrees.
            Default is 1 degrees.

        """
        
        # open the file, establish the tree, and its properties
        self.__sanitize_inputs(path_to_dir_file, path_to_ref_file)

        if ref_name is None:
            ref_name = radius

        if (radius < 0) or not np.isfinite(radius):
            raise ValueError(f"{radius} is invalid")

        if (angular_size < 0) or not np.isfinite(angular_size):
            raise ValueError(f"{angular_size} is invalid")

        the_correlator = ROOT.RayTraceCorrelator(self.station_id,
                                                 self.num_channels,
                                                 radius,
                                                 angular_size,
                                                 path_to_dir_file,
                                                 path_to_ref_file
                                                 )
        ROOT.SetOwnership(the_correlator, True) # take ownership
        the_correlator.LoadTables()
        self.correlators[ref_name] = the_correlator

