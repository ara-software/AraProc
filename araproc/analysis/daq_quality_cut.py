import numpy as np
import logging
import yaml
from araproc.framework import constants
import importlib.resources as pkg_resources
import araproc.framework.config_files as config_files

# Set up logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

def process_event_info(useful_event, station_id):

    """
    Process a useful_event to extract event information, determine trigger type, and calculate block length.

    Parameters
    ----------
    useful_event : object
        The useful event object containing event information.
    station_id : int
        The station ID used to determine pulser time.

    Returns
    -------
    blk_len : np.ndarray
        Calculated block length for the event.
    trig_type : int
        Trigger type (0: RF trigger, 1: Pulser trigger, 2: Soft trigger).
    block_nums : np.ndarray
        Block numbers for the event.
    channel_mask : np.ndarray
        Channel mask for the event.
    """

    # Extract block information
    num_blocks = len(useful_event.blockVec)
    block_nums = np.empty(num_blocks)  # Create an empty array for block numbers
    channel_mask = np.empty(num_blocks, dtype=int)  # Create an empty array for channel masks
    for n in range(num_blocks):
        this_block_id = useful_event.blockVec[n].getBlock()
        block_nums[n] = int(this_block_id)
        channel_mask[n] = useful_event.blockVec[n].channelMask
    # Extract event-level information
    read_win = useful_event.numReadoutBlocks  # Number of readout blocks

    if useful_event.isCalpulserEvent():
       trig_type = 1
    elif useful_event.isRFTrigger():
       trig_type = 0
    elif useful_event.isSoftwareTrigger():
       trig_type = 2
    elif station_id == 5: 
        # Station 5 force readouts from PA don't have a trigger type. Log them as an RF event
        trig_type = 0
    else:
       logger.error("__________________ Unknown trigger type for the event _______________________")  # error message
    num_ddas = constants.num_dda 
    blk_len = read_win // num_ddas

    return num_ddas,blk_len, trig_type, block_nums, channel_mask


def get_daq_structure_errors(num_ddas,blk_len_sort, trig_sort, irs_block_number, channel_mask):

    """
    Detect DAQ structure errors for a single event based on provided event data.

    Parameters
    ----------
    num_ddas : int
        Total number of Digitizing DAughter board.
    blk_len_sort : int
        The block lengths for the event.
    trig_sort : int
        trigger information for the event.
    irs_block_number : array-like
        The IRS block number for the event.
    channel_mask : array-like
        The channel mask for the event.

    Returns
    -------
    bool
        True if there are DAQ structure errors, False otherwise.
    """

    num_chs = constants.rf_chans_per_dda   # Number of channels
    num_blocks_per_dda = constants.num_blocks_per_dda
    bi_ch_mask = 1 << np.arange(num_chs, dtype=int)  # Bit mask for channels
    dda_ch = np.arange(num_ddas, dtype=int)  # DDA channels
    dda_idx = (channel_mask & 0x300) >> 8  # DDA index from channel mask
    max_blk_diff = -(num_blocks_per_dda - 1)  # Maximum block difference (wrap-around at 512)

    daq_st_err = np.full(5, 0, dtype=int)  # Error array for this single event

    # Block index for this event
    blk_idx_evt = np.asarray(irs_block_number, dtype=int)

    # Check for block length errors (remainder when divided by number of DDAs)
    daq_st_err[0] = len(blk_idx_evt) % num_ddas
    if daq_st_err[0] != 0:
        return True

    # Reshape block indices into groups of DDAs for comparison
    blk_idx_reshape = np.reshape(blk_idx_evt, (-1, num_ddas))
    daq_st_err[1] = int(np.any(blk_idx_reshape != blk_idx_reshape[:, 0][:, np.newaxis]))

    for dda in range(num_ddas):
        blk_idx_dda = blk_idx_reshape[:, dda]
        first_block_idx = blk_idx_dda[0]
        last_block_idx = blk_idx_dda[-1]
        block_diff = len(blk_idx_dda) - 1

        # Check block continuity, accounting for wrap-around at 512
        if first_block_idx + block_diff != last_block_idx:
            if num_blocks_per_dda - first_block_idx + last_block_idx != block_diff:
                daq_st_err[2] += 1

        # Check block index increments
        blk_diff = np.diff(blk_idx_dda).astype(int)
        incre_flag = np.any(np.logical_and(blk_diff != 1, blk_diff != max_blk_diff))
        if incre_flag:
            daq_st_err[2] += 1

    # Check DDA index consistency
    dda_idx_evt = np.asarray(dda_idx, dtype=int)
    dda_idx_reshape = np.reshape(dda_idx_evt, (-1, num_ddas))
    daq_st_err[3] = int(np.any(dda_idx_reshape != dda_ch[np.newaxis, :]))

    # Check channel mask consistency
    ch_mask_evt = np.asarray(channel_mask, dtype=int)
    ch_mask_reshape = np.repeat(ch_mask_evt[:, np.newaxis], num_chs, axis=1)
    ch_mask_bit = ch_mask_reshape & bi_ch_mask[np.newaxis, :]
    daq_st_err[4] = int(np.any(ch_mask_bit != bi_ch_mask[np.newaxis, :]))

    # Return True if there are any DAQ structure errors
    return np.sum(daq_st_err) > 0


def get_read_win_limit(st,config):

    """
    Detect readout window errors for a single event based on event data.

    Parameters
    ----------
    st : int
        Station id number.
    config : int
        Configuration Number
    Returns
    -------
    rf_readout_limit: int
        RF trig readout window limit.
    soft_readout_limit: int
        Soft trig readout window limit.
    """

    file = pkg_resources.open_text(config_files,"analysis_configs.yaml")
    file_content = yaml.safe_load(file)
    #Get the readout limits for different triggers
    readout_limits = file_content[f"station{st}"][f"config{config}"]["readout_limits"]
    rf_readout_limit = readout_limits["rf_readout_limit"]
    soft_readout_limit = readout_limits["soft_readout_limit"]
    file.close()
    return rf_readout_limit, soft_readout_limit

def get_readout_window_errors(blk_len_sort, trig_sort, channel_mask, st,config):

    """
    Detect readout window errors for a single event based on event data.

    Parameters
    ----------
    blk_len_sort : int
        The sorted block length for the event.
    trig_sort : int
        The trigger information for the event.
    channel_mask : array-like
        The channel mask for the event.
    config : int
        Configuration Number
    st : int
        The station ID for the event.

    Returns
    -------
    bool
        True if there are readout window errors, False otherwise.
    """

    # Get readout window limits for RF and software triggers
    rf_read_win_len, soft_read_win_len = get_read_win_limit(st,config)

    # Define boolean conditions for various readout errors
    single_read_bools = blk_len_sort < 2
    rf_cal_read_bools = blk_len_sort < rf_read_win_len
    rf_read_bools = np.logical_and(rf_cal_read_bools, trig_sort == 0)
    cal_read_bools = np.logical_and(rf_cal_read_bools, trig_sort == 1)
    soft_read_bools = np.logical_and(blk_len_sort != soft_read_win_len, trig_sort == 2)

    # Initialize the readout window error array for the single event (4 error types)
    read_win_err = np.full(4, 0, dtype=int)

    # Populate the error array based on the conditions
    read_win_err[0] = int(single_read_bools)  # Single block error
    read_win_err[1] = int(rf_read_bools)      # Bad RF readout window error
    read_win_err[2] = int(cal_read_bools)     # Bad calibration readout window error
    read_win_err[3] = int(soft_read_bools)    # Bad software readout window error

    # Return True if there are any readout window errors
    return np.sum(read_win_err) > 0


def check_daq_quality(useful_event, station_id, config):

    """
    Main function to extract DAQ structure and readout window errors for a single event.

    Parameters
    ----------
    useful_event : object
        The useful event object containing event information.
    station_id : int
        The station ID for the event.
    config : int
        Configuration Number

    Returns
    -------
    combined_errors : bool
        True if there are DAQ structure errors or readout window errors, False otherwise.
        
    """

    # Extract necessary information from useful_event using the merged function
    num_ddas,blk_len, trig_type, irs_block_number, channel_mask = process_event_info(useful_event, station_id)

    # Get DAQ structure errors
    daq_errors = get_daq_structure_errors(num_ddas,blk_len, trig_type, irs_block_number, channel_mask)

    # Get readout window errors
    readout_errors = get_readout_window_errors(blk_len, trig_type, channel_mask,station_id,config)

    combined_errors = daq_errors or readout_errors

    return combined_errors

