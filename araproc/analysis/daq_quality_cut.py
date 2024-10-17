import numpy as np
import os
import uproot
import ROOT
import logging

# Set up logging configuration
logging.basicConfig(level=logging.INFO,  # Change to DEBUG for more detailed logs
                    format='%(message)s')

logger = logging.getLogger(__name__)

# Load ARA libraries
ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR') + "/lib/libAraEvent.so")
ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR') + "/lib/libRootFftwWrapper.so.3.0.1")

def get_sub_info(data):

    """
    Extract relevant event information from a ROOT file.
    Returns several arrays with sorted event details and metadata.
    """

    file = uproot.open(data)
    evtTree = file['eventTree']
    st_arr = np.asarray(evtTree['event/RawAraStationEvent/RawAraGenericHeader/stationId'], dtype=int)
    station_id = st_arr[0]
    num_evts = len(st_arr)
    run = int(np.asarray(evtTree['run'], dtype=int)[0])
    evt_num = np.asarray(evtTree['event/eventNumber'], dtype=int)
    unix_time = np.asarray(evtTree['event/unixTime'], dtype=int)
    read_win = np.asarray(evtTree['event/numReadoutBlocks'], dtype=int)
    irs_block_number = np.asarray(evtTree['event/blockVec/blockVec.irsBlockNumber']) & 0x1ff
    channel_mask = np.asarray(evtTree['event/blockVec/blockVec.channelMask'])
    time_stamp = np.asarray(evtTree['event/timeStamp'], dtype=int)
    trigger_info = np.asarray(evtTree['event/triggerInfo[4]'], dtype=int)[:, 2]

    # Sort event data
    evt_sort_idx = np.argsort(evt_num)
    evt_num_sort = evt_num[evt_sort_idx]
    unix_time_sort = unix_time[evt_sort_idx]
    trig_type = get_trig_type(time_stamp, trigger_info, station_id)
    trig_type_sort = trig_type[evt_sort_idx]

    blk_len = get_block_length(read_win, trig_type)
    blk_len_sort = blk_len[evt_sort_idx]

    return (blk_len_sort, trig_type_sort, evt_num_sort, unix_time_sort, 
            station_id, run, num_evts, evt_num, unix_time, read_win, 
            irs_block_number, channel_mask, time_stamp)

def get_trig_type(time_stamp, trigger_info, station_id):

    """
    Determine the trigger type based on time stamps and trigger information.
    Returns an array indicating the type of trigger for each event.
    """

    pulserTime = np.array([254, 245, 245, 400, 400], dtype=int)
    trig_type = np.full(len(time_stamp), 0, dtype=int)

    trig_type[np.abs(time_stamp - pulserTime[station_id - 1]) < 1e4] = 1
    trig_type[trigger_info == 1] = 2

    return trig_type

def get_block_length(read_win, trig_type):

    """
    Calculate the block length for different trigger types.
    Returns the basic block length.
    """

    num_ddas = 4
    blk_len = read_win // num_ddas

    rf_blk_len = np.copy(blk_len).astype(float)
    rf_blk_len[trig_type != 0] = np.nan
    cal_blk_len = np.copy(blk_len).astype(float)
    cal_blk_len[trig_type != 1] = np.nan
    soft_blk_len = np.copy(blk_len).astype(float)
    soft_blk_len[trig_type != 2] = np.nan

    return blk_len

def get_daq_structure_errors(blk_len_sort, trig_sort, irs_block_number, channel_mask, run, num_evts):

    """
    Detect DAQ structure errors based on provided event data.
    Returns an array indicating the types of DAQ errors for each event.
    """

    num_ddas = 4  # Update based on your context
    num_chs = 8   # Update based on your context
    bi_ch_mask = 1 << np.arange(num_chs, dtype=int)
    dda_ch = np.arange(num_ddas, dtype=int)
    dda_idx = (channel_mask & 0x300) >> 8
    max_blk_diff = -(512 - 1)  # Update based on your context
    daq_st_err = np.full((num_evts, 5), 0, dtype=int)

    for evt in range(num_evts):
        blk_idx_evt = np.asarray(irs_block_number[evt], dtype=int)
        daq_st_err[evt, 0] = len(blk_idx_evt) % num_ddas
        if daq_st_err[evt, 0] != 0:
            continue

        blk_idx_reshape = np.reshape(blk_idx_evt, (-1, num_ddas))
        daq_st_err[evt, 1] = int(np.any(blk_idx_reshape != blk_idx_reshape[:, 0][:, np.newaxis]))

        for dda in range(num_ddas):
            blk_idx_dda = blk_idx_reshape[:, dda]
            first_block_idx = blk_idx_dda[0]
            last_block_idx = blk_idx_dda[-1]
            block_diff = len(blk_idx_dda) - 1

            if first_block_idx + block_diff != last_block_idx:
                if 512 - first_block_idx + last_block_idx != block_diff:
                    daq_st_err[evt, 2] += 1

            blk_diff = np.diff(blk_idx_dda).astype(int)
            incre_flag = np.any(np.logical_and(blk_diff != 1, blk_diff != max_blk_diff))
            if incre_flag:
                daq_st_err[evt, 2] += 1

        dda_idx_evt = np.asarray(dda_idx[evt], dtype=int)
        dda_idx_reshape = np.reshape(dda_idx_evt, (-1, num_ddas))
        daq_st_err[evt, 3] = int(np.any(dda_idx_reshape != dda_ch[np.newaxis, :]))

        ch_mask_evt = np.asarray(channel_mask[evt], dtype=int)
        ch_mask_reshape = np.repeat(ch_mask_evt[:, np.newaxis], num_chs, axis=1)
        ch_mask_bit = ch_mask_reshape & bi_ch_mask[np.newaxis, :]
        daq_st_err[evt, 4] = int(np.any(ch_mask_bit != bi_ch_mask[np.newaxis, :]))

        if sum(daq_st_err[evt, :]) > 0:
            logger.warning('Bad block index: Run %d, Event %d, Errors: %s, Total: %d', run, evt, daq_st_err[evt, :], sum(daq_st_err[evt, :]))

    return daq_st_err


def get_read_win_limit(st, run):

    """
    Determine the readout window limits based on station and run number.
    Returns the RF and soft readout limits.
    """

    if st == 1:             ############ Once you insert original limits for your station please remove this comment ########
        if run < 4029:
            rf_readout_limit = 20
        elif run > 4028 and run < 9749:
            rf_readout_limit = 26
        elif run > 9748:
            rf_readout_limit = 28
    if st == 2:              ############ Once you insert original limits for your station please remove this comment ########
        if run < 4029:
            rf_readout_limit = 20
        elif run > 4028 and run < 9749:
            rf_readout_limit = 26
        elif run > 9748:
            rf_readout_limit = 28
    elif st == 3:             ############ Once you insert original limits for your station please remove this comment ########
        if run < 3104:
            rf_readout_limit = 20
        elif run > 3103 and run < 10001:
            rf_readout_limit = 26
        elif run > 10000:
            rf_readout_limit = 28
    elif st == 4:
        if run < 3600:
            rf_readout_limit = 28
        elif run > 3600:
            rf_readout_limit = 32
    elif st == 5:              ############ Once you insert original limits for your station please remove this comment ########
        if run < 3600:
            rf_readout_limit = 28
        elif run > 3600:
            rf_readout_limit = 32

    if st == 1:                 ############ Once you insert original limits for your station please remove this comment ########
        if run < 9505:
            soft_readout_limit = 8
        else:
            soft_readout_limit = 12
    if st == 2:                  ############ Once you insert original limits for your station please remove this comment ########
        if run < 9505:
            soft_readout_limit = 8
        else:
            soft_readout_limit = 12
    elif st == 3:                ############ Once you insert original limits for your station please remove this comment ########
        if run < 10001:
            soft_readout_limit = 8
        else:
            soft_readout_limit = 12
    elif st == 4:
        soft_readout_limit = 10
    elif st == 5:                ############ Once you insert original limits for your station please remove this comment ########
        soft_readout_limit = 10

    return rf_readout_limit, soft_readout_limit

def get_readout_window_errors(blk_len_sort, trig_sort, channel_mask, run, st, num_evts):

    """
    Detect readout window errors based on the event data.
    Returns an array indicating readout window errors for each event.
    """

    rf_read_win_len, soft_read_win_len = get_read_win_limit(st, run)
    single_read_bools = blk_len_sort < 2
    rf_cal_read_bools = blk_len_sort < rf_read_win_len
    rf_read_bools = np.logical_and(rf_cal_read_bools, trig_sort == 0)
    cal_read_bools = np.logical_and(rf_cal_read_bools, trig_sort == 1)
    soft_read_bools = np.logical_and(blk_len_sort != soft_read_win_len, trig_sort == 2)

    bad_single_evts = np.where(single_read_bools)[0]
    bad_rf_evts = np.where(rf_read_bools)[0]
    bad_cal_evts = np.where(cal_read_bools)[0]
    bad_soft_evts = np.where(soft_read_bools)[0]

    read_win_err = np.full((num_evts, 4), 0, dtype=int)
    read_win_err[:, 0] = np.in1d(np.arange(num_evts), bad_single_evts).astype(int)  # single block
    read_win_err[:, 1] = np.in1d(np.arange(num_evts), bad_rf_evts).astype(int)  # bad rf readout window
    read_win_err[:, 2] = np.in1d(np.arange(num_evts), bad_cal_evts).astype(int)  # bad cal readout window
    read_win_err[:, 3] = np.in1d(np.arange(num_evts), bad_soft_evts).astype(int)  # bad soft readout window

    return read_win_err

def get_all_errors(file_path):

    """
    Main function to analyze a ROOT file and return DAQ and readout window errors.
    """

    blk_len_sort, trig_sort, evt_sort, unix_sort, st, run, num_evts, evt_num, unix_time, read_win, irs_block_number, channel_mask, time_stamp = get_sub_info(file_path)

    daq_errors = get_daq_structure_errors(blk_len_sort, trig_sort, irs_block_number, channel_mask, run, num_evts)
    readout_errors = get_readout_window_errors(blk_len_sort, trig_sort, channel_mask, run, st, num_evts)

    return daq_errors, readout_errors

