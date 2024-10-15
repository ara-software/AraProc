import ROOT
import numpy as np
import math
import os

# HELP HERE: constants that we might want to store elsewhere, or can I retrieve them somehow?
channel_num = 16

def get_surface_correlation_ratio(corr_map, station_id, z_thresh):
    """
    Calculates the ratio of the maximum surface correlation value within a specific
    theta range to the overall peak correlation value. The theta range is determined
    based on the average antenna z-coordinate and a threshold.

    Parameters
    ----------
    corr_map : dict
        A dictionary containing:
        - 'map': ROOT.TH2D histogram with correlation values mapped by theta and phi.
        - 'corr': The overall peak correlation value.
        - 'radius': The correlation radius for the station.
    station_id : int
        The station ID for which to calculate the correlation ratio.
    z_thresh : float
        The threshold value to define the lower bound of the theta range.

    Returns
    -------
    surf_corr_ratio : float
        The ratio of the maximum surface correlation value within the theta range 
        to the overall peak correlation.
    max_surf_theta : float
        The theta angle corresponding to the maximum surface correlation value.
    max_surf_phi : float
        The phi angle corresponding to the maximum surface correlation value.
    """
    
    # Extract the radius from corr_map
    corr_radius = float(corr_map.get("radius", None))
    if corr_radius is None:
        raise ValueError("Radius not found in corr_map.")

    # Initialize AraGeomTool instance to retrieve antenna information
    geomTool = ROOT.AraGeomTool.Instance()

    # Collect z-coordinates of all antennas to calculate the average z-coordinate
    z = []
    for ant in range(channel_num):
        # Retrieve each antenna's location for the given station
        ant_location = geomTool.getStationInfo(station_id).getAntennaInfo(ant).antLocation
        z.append(ant_location[2])  # Append only the z-coordinate

    # Calculate the average z-coordinate of the antennas
    avg_z = sum(z) / len(z)

    # Check if the radius is large enough to visualize the surface
    if corr_radius < abs(avg_z):
        raise RuntimeError("Surface not visible: radius is smaller than average z coordinate of the antennas.")

    # Calculate the angles at the surface and z threshold relative to the horizontal
    theta_at_surface = math.degrees(math.asin(abs(avg_z) / corr_radius))
    theta_at_z_thresh = math.degrees(math.asin((abs(avg_z) - z_thresh) / corr_radius))

    # Access the ROOT.TH2D histogram from 'map' in corr_map
    map = corr_map.get("map", None)
    if map is None:
        raise ValueError("The 'map' key was not found in corr_map.")

    # Initialize lists to store correlation values within the specified theta range
    surf_corr = []
    tot_corr_idx = []

    # Define the theta range based on theta_at_surface and theta_at_z_thresh
    theta_min = min(theta_at_surface, theta_at_z_thresh)
    theta_max = max(theta_at_surface, theta_at_z_thresh)

    # Iterate over each bin in the map to filter by theta range
    nbins_x = map.GetNbinsX()
    nbins_y = map.GetNbinsY()
    for binx in range(1, nbins_x + 1):
        for biny in range(1, nbins_y + 1):
            # Get theta and phi for the current bin and retrieve its correlation value
            theta = map.GetYaxis().GetBinCenter(biny)
            phi = map.GetXaxis().GetBinCenter(binx)
            corr_value = map.GetBinContent(binx, biny)

            # Store correlation values and bin indices if within the theta range
            if theta_min <= theta <= theta_max:
                surf_corr.append(corr_value)
                tot_corr_idx.append((binx, biny))

    # Ensure that there are correlation values within the theta range
    if not surf_corr:
        raise RuntimeError("No correlation values found in the specified theta range!")

    # Calculate the maximum correlation within the specified theta range
    max_surf_corr = np.nanmax(surf_corr)
    max_surf_corr_index = np.nanargmax(surf_corr)

    # Identify bin location of the maximum correlation value
    max_surf_bin = tot_corr_idx[max_surf_corr_index]
    max_surf_bin_x, max_surf_bin_y = max_surf_bin
    max_surf_theta = map.GetYaxis().GetBinCenter(max_surf_bin_y)
    max_surf_phi = map.GetXaxis().GetBinCenter(max_surf_bin_x)

    # Access the overall peak correlation value from 'corr' in corr_map
    corr_value = corr_map.get("corr", None)
    if corr_value is None:
        raise ValueError("The 'corr' key was not found in corr_map.")

    # Calculate the ratio of the max surface correlation to the peak correlation
    surf_corr_ratio = max_surf_corr / corr_value

    return surf_corr_ratio, max_surf_theta, max_surf_phi
