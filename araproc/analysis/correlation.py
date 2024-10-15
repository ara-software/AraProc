import ROOT
import numpy as np
import math
import os

# constants that we might want to store elsewhere or can I retrieve it somehow?
channel_num = 16

def get_surface_correlation_ratio(corr_map, station_id, z_thresh):
    # Extract the radius from corr_map
    corr_radius = float(corr_map.get("radius", None))
    if corr_radius is None:
        raise ValueError("Radius not found in corr_map.")

    # Now extract the z distance of the station (antenna z- coords avg) from the surface

    # create a geomtool
    geomTool = ROOT.AraGeomTool.Instance()

    # Lists to store coordinates
    z = []

    # Iterate over all antennas
    for ant in range(channel_num):
        # Retrieve the location of each antenna
        ant_location = geomTool.getStationInfo(station_id).getAntennaInfo(ant).antLocation
    
        # Append the z coordinates to the list
        z.append(ant_location[2])

    # Calculate the average of z coordinates
    avg_z = sum(z) / len(z)

    # DELETE LATER Print avg_z
    print(f"Average z: {avg_z}")

    # Check if the radius is smaller than avg_z
    if corr_radius < abs(avg_z):
        raise RuntimeError("Surface not visible: radius is smaller than average z coordinate of the antennas.")

    # Calculate angles with respect to the horizontal
    theta_at_surface = math.degrees(math.asin(abs(avg_z) / corr_radius))
    theta_at_z_thresh = math.degrees(math.asin((abs(avg_z)-z_thresh)/ corr_radius))

    print("theta_at_surface", theta_at_surface)
    print("theta_at_z_thersh", theta_at_z_thresh)

    # Access the ROOT.TH2D histogram from the 'map' key in corr_map
    map = corr_map.get("map", None)
    if map is None:
        raise ValueError("The 'map' key was not found in corr_map.")

    # Initialize lists to store correlation values and their corresponding bin indices
    surf_corr = []
    tot_corr_idx = []

    # Define theta range
    theta_min = min(theta_at_surface, theta_at_z_thresh)
    theta_max = max(theta_at_surface, theta_at_z_thresh)

    # Loop over all bins in the map
    nbins_x = map.GetNbinsX()
    nbins_y = map.GetNbinsY()

    for binx in range(1, nbins_x + 1):
        for biny in range(1, nbins_y + 1):
            # Get the theta and phi for the current bin
            theta = map.GetYaxis().GetBinCenter(biny)
            phi = map.GetXaxis().GetBinCenter(binx)
            corr_value = map.GetBinContent(binx, biny)

            # Check if theta is within the desired range
            if theta_min <= theta <= theta_max:
                # Store the correlation value and the current bin (binx, biny)
                surf_corr.append(corr_value)
                tot_corr_idx.append((binx, biny))

    # Raise an error if no correlations are found in the specified range
    if not surf_corr:
        raise RuntimeError("No correlation values found in the specified theta range!")

    # Find the maximum correlation in the filtered list
    max_surf_corr = np.nanmax(surf_corr)
    max_surf_corr_index = np.nanargmax(surf_corr)

    # Get the bin indices corresponding to the maximum correlation
    max_surf_bin = tot_corr_idx[max_surf_corr_index]
    max_surf_bin_x, max_surf_bin_y = max_surf_bin

    # Retrieve the corresponding theta and phi values for the maximum correlation
    max_surf_theta = map.GetYaxis().GetBinCenter(max_surf_bin_y)
    max_surf_phi = map.GetXaxis().GetBinCenter(max_surf_bin_x)

    # Print results
    print(f"Maximum correlation value within theta range: {max_surf_corr}")
    print(f"Bin indices with max correlation: (x: {max_surf_bin_x}, y: {max_surf_bin_y})")
    print(f"Corresponding theta: {max_surf_theta}")
    print(f"Corresponding phi: {max_surf_phi}")
    print(f"Theta range: ({theta_min}, {theta_max})")

    # Access the "corr" value from corr_map
    corr_value = corr_map.get("corr", None)
    if corr_value is None:
        raise ValueError("The 'corr' key was not found in corr_map.")

    print("mac corr value: ", corr_value)

    surf_corr_ratio = max_surf_corr / corr_value

    return surf_corr_ratio
