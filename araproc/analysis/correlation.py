import ctypes

def get_max_corr_and_coords(corr_map):
    """
    Finds the maximum correlation value in a single map and extracts the 
    corresponding theta and phi.

    Parameters
    ----------
    corr_map : ROOT.TH2D
        A 2D histogram containing correlation values, where the axes represent
        theta and phi angles.

    Returns
    -------
    max_corr : float
        The maximum correlation value in the map.
    max_theta : float
        The theta coordinate corresponding to the maximum correlation value.
    max_phi : float
        The phi coordinate corresponding to the maximum correlation value.
    """
    # Find the bin with the maximum correlation
    max_bin = corr_map.GetMaximumBin()
    binx = ctypes.c_int()
    biny = ctypes.c_int()
    binz = ctypes.c_int() # Unused but needed
    corr_map.GetBinXYZ(max_bin, binx, biny, binz)

    # Get the correlation value and corresponding theta, phi
    max_corr = corr_map.GetBinContent(max_bin)
    max_theta = corr_map.GetYaxis().GetBinCenter(biny.value)
    max_phi = corr_map.GetXaxis().GetBinCenter(binx.value)

    return max_corr, max_theta, max_phi