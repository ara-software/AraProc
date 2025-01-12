import numpy as np

# Physical constants

speed_of_light = 0.2998 # m/ns

# AraProc global constants

valid_station_ids = [100, 2, 3, 4, 5]

rf_channels_ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
vpol_channel_ids = [0, 1, 2, 3, 4, 5, 6, 7]
hpol_channel_ids = [8, 9, 10, 11, 12, 13, 14, 15]

num_string = 4
num_dda = 4
rf_chans_per_dda = 8
num_antennas_per_dda = 4
num_blocks_per_dda = 512
num_samples_per_block = 64
num_samples_per_dda = 32768
num_polarization = 2
num_trigger_type = 3
num_rf_channels = len(rf_channels_ids)

# ice model
def get_index_of_refraction(z): 
    """
    Computes the index of refraction for a depth z.

    Parameters
    ----------
    z : float
    Height relative to ice surface in meters (below surface is negative).
    
    Returns
    -------
    n : float
    Index of refraction
    """

    # in-air
    if z > 0:
        return 1

    #in-ice

    # ice model parameters from PA paper
    A = 1.780
    B = 0.454
    C = 0.0202

    n = A - B*np.exp(-C*np.abs(z))

    return n


