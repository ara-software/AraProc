import ctypes
import ROOT

def get_corr_map_peak(the_map  = None):
        
        """
        Little utility helper function

        Parameters
        ----------
        the_map : ROOT.TH2D
            A 2D map you want the peak found
        
        
        Returns
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

def calculate_avg_antenna_xyz(station_id, num_channels):
    """
    Calculates the average x, y, and z coordinates of antennas for a given station.

    Parameters
    ----------
    station_id : int
        ID of the station from which to retrieve antenna information.
    num_channels : int
        Number of channels (antennas) at the station.

    Returns
    -------
    avg_x : float
        Average x-coordinate of the antennas at the specified station.
    avg_y : float
        Average y-coordinate of the antennas at the specified station.
    avg_z : float
        Average z-coordinate of the antennas at the specified station.
    """
    geom_tool = ROOT.AraGeomTool.Instance()
    
    # Initialize sums for x, y, and z coordinates
    sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
    
    # Sum up x, y, and z coordinates for all antennas
    for ant in range(num_channels):
        ant_info = geom_tool.getStationInfo(station_id).getAntennaInfo(ant).antLocation
        sum_x += ant_info[0]
        sum_y += ant_info[1]
        sum_z += ant_info[2]
    
    # Calculate averages
    avg_x = sum_x / num_channels
    avg_y = sum_y / num_channels
    avg_z = sum_z / num_channels
    
    return avg_x, avg_y, avg_z
