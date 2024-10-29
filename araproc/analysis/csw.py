import numpy as np
from scipy.interpolate import Akima1DInterpolator
from scipy.signal import correlate
import ROOT

from araproc.analysis import standard_reco as sr
from araproc.framework.dataset import AnalysisDataset
from araproc.framework import waveform_utilities as wu

def _get_abspeak(xs, ys):
    """
    Given x and y values, find the peak (of the absolute value) of the y array.
    Return the corresponding value in the `xs` array followed by the peak in 
      the `ys` array.
    """
    max_idx = np.nanargmax( np.abs(ys) )
    return xs[max_idx], ys[max_idx]

def _get_peak(xs, ys):
    """
    Given x and y values, find the peak of the y array.
    Return the corresponding value in the `xs` array followed by the peak in 
      the `ys` array.
    """
    max_idx = np.nanargmax( ys )
    return xs[max_idx], ys[max_idx]

def _get_channels_to_csw(excluded_channels, data, polarization):
    """
    Return an array of RF channel IDs that will be used in the CSW based on the 
      user-provided `excluded_channels` and the polarization of the antennas.
    """

    # Determine which channels should be used in the csw
    if excluded_channels is not None: 
        if not isinstance(excluded_channels, np.ndarray):
            excluded_channels = np.array(excluded_channels)
    else: 
        excluded_channels = np.array([])
    
    channels_to_csw = []
    for ch in data.rf_channel_indices: 
        if (ch not in excluded_channels) and (data.rf_channel_polarizations[ch] == polarization):
            channels_to_csw.append(ch)

    return channels_to_csw[:]

def _get_arrival_delays_reco(
    data, reco, reco_results, channels_to_csw, reference_ch, 
    which_distance, solution
):
    """
    Extract the arrival time of the signals from the reconstruction then 
      calculate the time delays between each channel and the reference channel
    """

    # dv.plot_skymap(
    #     reco_results[reco_key]["map"],
    #     ouput_file_path=f"./station_{data.station_id}_run_{data.run_number}_event_{e}_map.png")

    # Get all the arrival delays
    arrival_times = {ch: 0.0 for ch in channels_to_csw }
    for c, ch_ID in enumerate(channels_to_csw):
        arrival_times[ch_ID] = reco.lookup_arrival_time(
            ch_ID, reco_results['theta'], reco_results['phi'], 
            which_distance=which_distance, solution=solution
        )
    
    # Center all arrival delays around a reference channel
    reference_arrival_time = arrival_times[reference_ch]
    arrival_delays = {ch: 0.0 for ch in channels_to_csw }
    for ch_ID in arrival_delays.keys():
        arrival_delays[ch_ID] = arrival_times[ch_ID] - reference_arrival_time

    return arrival_delays

def _get_arrival_delays_xcorr(
    wavepacket, reco, channels_to_csw, reference_ch):
    """
    Use the cross correlation results already calculated by the reconstruction 
      algorithm to determine arrival delays
    """
    
    arrival_delays = {ch: 0. for ch in channels_to_csw }
    for ch_ID in arrival_delays.keys():
        if ch_ID == reference_ch: continue
        cross_corr = reco._StandardReco__get_correlation_function(ch_ID, reference_ch, wavepacket)
        arrival_delays[ch_ID] = _get_peak(cross_corr.GetX(), cross_corr.GetY())[0]

    return arrival_delays

def _get_arrival_delays_AraRoot_xcorr(
    waveform_bundle, reco, pol, reference_ch, reco_delays, wf_type, zoom_window=20):
    """
    Determine the arrival delays from each channel by finding the time of
      max correlation between each channel and the reference channel in a 
      +/- `zoom_window` nanosecond window around the expected arrival delay of a 
      signal given the results from reconstructing the event.
    """

    # 
    pairs = reco.pairs_v if pol==0 else reco.pairs_h
    channels_in_pairs = np.unique([[ch1, ch2] for (idx, (ch1, ch2)) in pairs])

    wf_map = ROOT.std.map("int", "TGraph*")()
    ROOT.SetOwnership(wf_map, True)
    for chan_i in waveform_bundle.keys():
        wf_map[chan_i] = waveform_bundle[chan_i]

    cross_corr_unhilberted = reco.rtc_wrapper.correlators["nearby"].GetCorrFunctions(
        pairs, wf_map, False)

    delays = {}
    for ch_ID in channels_in_pairs:

        if ch_ID == reference_ch: 
            delay = 0
        else: 
            pair_idx = reco.get_pair_index(ch_ID, reference_ch, pairs)
            xcorr = cross_corr_unhilberted[pair_idx]
            zoomed_indices = np.where(
                (np.abs( np.asarray(xcorr.GetX()) - reco_delays[ch_ID] )) < zoom_window
            )[0]
            if len(zoomed_indices) == 0: 
                if wf_type != "soft": 
                    print(
                        f"Touble calculating csw for {wf_type} event with "
                        f"channels {ch_ID} and {reference_ch}")
                delay = _get_peak(
                    np.asarray(xcorr.GetX()), 
                    np.asarray(xcorr.GetY())
                )[0] 
            else: 
                print(
                    ch_ID, len(zoomed_indices), xcorr.GetN(), 
                    round(reco_delays[ch_ID], 1), 
                    np.nanargmax( np.asarray( xcorr.GetY() ) ),
                    np.round(_get_peak(
                        np.asarray(xcorr.GetX())[zoomed_indices], 
                        np.asarray(xcorr.GetY())[zoomed_indices]
                    )[0], 1)
                )
                print(ch_ID, xcorr.GetX()[zoomed_indices[0]], xcorr.GetX()[zoomed_indices[1]],
                    xcorr.GetY()[zoomed_indices[0]], xcorr.GetY()[zoomed_indices[-1]])
                delay = _get_peak(
                    np.asarray(xcorr.GetX())[zoomed_indices], 
                    np.asarray(xcorr.GetY())[zoomed_indices]
                )[0] 

        # AraRoot always compared channels with smaller IDs to channels with 
        #   larger IDs but we always want to compare to the reference channel.
        #   Correct for this if the current ch_ID is larger than the ref. ch_ID
        if ch_ID > reference_ch: 
            delay *= -1
        
        delays[ch_ID] = delay

    return delays

def _get_arrival_delays(waveforms, times, channels_to_csw, reference_ch):

    arrival_delays = {ch_ID: None for ch_ID in channels_to_csw}

    volts_ref_nan_idx = np.isnan(waveforms[reference_ch])
    volts_ref_nonan = np.copy( waveforms[reference_ch] )
    volts_ref_nonan[volts_ref_nan_idx] = 0.0

    for ch_ID in channels_to_csw: 

        volts_nan_idx = np.isnan(waveforms[ch_ID])
        volts_nonan = np.copy(waveforms[ch_ID])
        volts_nonan[volts_nan_idx] = 0.0

        xcorr = correlate(volts_nonan, volts_ref_nonan, method='direct')
        xcorr_idxpeak, xcorr_vpeak = _get_peak(np.arange(len(xcorr)), xcorr)
        arrival_delay_bin_shift = xcorr_idxpeak - len(xcorr)//2
        arrival_delays_time = (
            np.sign(arrival_delay_bin_shift) * 
            ( times[abs(arrival_delay_bin_shift)] - times[0] ))
        arrival_delays[ch_ID] = arrival_delay_bin_shift

    return arrival_delays

def _get_catchall_time_array(waveforms, channels_to_csw, arrival_delays):
    """
    For a list of antennas with waveforms that will be shifted by the provided 
      `arrival_delays`, build a time array from the earliest waveform time
      to the latest waveform time, using the same binning as the waveforms.
    """

    # Get a time array that stretches from the earliest (arrival delay shifted)
    #   timestamp to the latest. 
    earliest_wf_time = 123456
    latest_wf_time = -123456
    for c, ch_ID in enumerate(channels_to_csw):
        waveform_times = np.asarray(waveforms[ch_ID].GetX())
        if (waveform_times[0] + arrival_delays[ch_ID]) < earliest_wf_time:
            earliest_wf_time = (waveform_times[0] + arrival_delays[ch_ID])
        if (waveform_times[-1] + arrival_delays[ch_ID]) > latest_wf_time:
            latest_wf_time = (waveform_times[-1] + arrival_delays[ch_ID])
        del waveform_times
    waveform_dt = waveforms[ch_ID].GetX()[1] - waveforms[ch_ID].GetX()[0]
    big_times = np.arange(
        earliest_wf_time, latest_wf_time+waveform_dt, waveform_dt)
    
    return big_times

def _shift_wf(xs, ys, x_shift, new_xs):
    """
    Shift a waveform with `volts` and `times` so that the waveform now exists
        along the axis `new_times` using an Akima interpolation..
    
    Parameters
    ----------
    volts : np.ndarray
        The original voltage array for this waveform.
    times : np.ndarray
        The original time array for this waveform.
    new_times : np.ndarray
        The requested new time array for this waveform

    Returns
    -------
    new_volts : np.ndarray
        The new voltage array corresponding to the requested `new_time` array.
    """

    # Interpolate the new wf
    wf_interpolation = Akima1DInterpolator(xs+x_shift, ys)
    new_ys = wf_interpolation(new_xs)

    return new_ys

def _get_channel_with_max_V(waveforms, waveforms_to_analyze):
    """
    Out of a given waveforms for different channels, 
      find and return the one with the greatest maximum voltage
    """

    max_voltage = {'voltage': -1, 'channel':-123456}
    for ch_ID in waveforms_to_analyze:
        this_max_voltage = np.max(waveforms[ch_ID].GetY())
        if this_max_voltage > max_voltage['voltage']:
            max_voltage['voltage'] = this_max_voltage
            max_voltage['channel'] = ch_ID 

    return max_voltage['channel']

def _get_shifts(waveforms, reference_ch):
    """
    For each waveform, calculate the number of bins you need to roll by.

    This function will fail if any waveforms passed to it are entirely 
      composed of `nan`s

    Parameters
    ----------
    waveforms : dict
        Dictionary where keys are channel IDs and values are the corresponding
            waveform as TGraphs
    n_sols : int
        Number of ray solutions to solve for

    Returns
    -------
    roll_shifts : np.ndarray
        The amount to roll each waveform by for each ray solution
        Dimensions: ( n_channels, n_sols )
    """

    ant_first_nonnans = {
        ch_ID: np.where( ~np.isnan(wf) )[0][0] 
        for ch_ID, wf in waveforms.items()
    }

    # Determine all shifts for rolling based on the earliest waveform
    roll_shifts = {
        ch_ID: this_ant_first_nonnan - ant_first_nonnans[reference_ch]
        for ch_ID, this_ant_first_nonnan in ant_first_nonnans.items()
    }

    return roll_shifts

def _roll_wf(wf, shift, times):
    """
    Roll the provided waveform. 

    More Info
    ---------
    It is assumed that the waveform has some number of `nan` values leading 
      up to and following the actual waveform. There is some trickiness to 
      properly roll the function in this scenario so I have an example where 
      `.` represents a `nan` value, the  `T` is the sample in the waveforms 
      that have been aligned and should remain aligned, and the other letters 
      are random waveform values.
    We want non-nan samples of this waveform
        A: ........acTeunw.....
    to be rolled to match this waveform
        B: ......rseuTca.......
    We need to roll by 2 samples but if we roll the whole waveform A 
      as is, it will look like this
        A: ......acTeunw.......
    which we don't want because the peak `T` gets misaligned. So, first, the 
      non-nan parts of waveform A are rolled forward by 2 samples to get this
        A: ........nwacTeu.....
    The peaks are now misaligned so we roll the entire waveform backwards
      now by 2 samples to get this
        A: ......nwacTeu.......
    So the timing and location of the peak should now match with waveform B
        B: ......rseuTca.......
    """

    # The interpolator uses nan's for extrapolated values that we don't want 
    #   to roll as a part of the entire waveform. 
    # Isolate the first and first-from-last non-nan value 
    #   then roll just the waveform forward by the requested `shift` and save
    idx_non_pad_values = ~np.isnan(wf)
    start_idx = np.where( 
        idx_non_pad_values.any(), idx_non_pad_values.argmax(), -123456)
    end_idx = np.where(
        idx_non_pad_values.any(),
        len(wf) - np.flip(idx_non_pad_values).argmax() - 1,
        -123456
    )
    rolled_nonnan_data = np.roll(
        wf[ start_idx : end_idx+1 ], shift
    )
    wf[ start_idx : end_idx+1 ] = rolled_nonnan_data

    # The peak is now offset from where if should be by `+shift`
    # Roll the entire waveform by `-shift` to move the peak back
    wf = np.roll(wf, -shift)

    return wf

def get_csw_reco(e,
    data : AnalysisDataset, 
    useful_event, 
    solution : int ,
    polarization : int,
    excluded_channels = None,
    which_distance : str = "distant",
    return_rolled_wfs : bool = False
):
    """
    From a given `useful_event`, return the coherently summed waveform (CSW)
    """

    # Pull wavepacket from the dataset
    wavepacket = data.get_wavepacket(useful_event=useful_event, which_traces='filtered')

    # Make a list of channels to analyze
    channels_to_csw = _get_channels_to_csw(excluded_channels, data, polarization)

    # Get reconstruction results
    reco = sr.StandardReco(data.station_id, excluded_channels=excluded_channels)
    reco_results = reco.do_standard_reco(wavepacket)
    pol_int_to_str = {0: "v", 1: "h"}
    distance_to_key = {'distant': 'distant', 'nearby': 'pulser'}
    reco_key = f"{distance_to_key[which_distance]}_{pol_int_to_str[polarization]}"

    # Calculate the arrival delays in each channel relative to some 
    #   reference channel
    reference_ch = _get_channel_with_max_V(wavepacket['waveforms'], channels_to_csw)
    arrival_delays = _get_arrival_delays_reco(
        data, reco, reco_results[reco_key], channels_to_csw, reference_ch, 
        which_distance, solution)
    
    # Shift each waveform by the arrival delay then interpolate each waveform 
    #   so it exists within an array that stretches from the earliest wf 
    #   time to the last wf time. Waveform values are set to `nan` at all 
    #   times for which it does not exist. These will be trimmed later.
    big_times = _get_catchall_time_array(wavepacket['waveforms'], channels_to_csw, arrival_delays)
    interpolated_wfs = {
        ch_ID: np.full(len(big_times), -123456.0) for ch_ID in channels_to_csw}
    for ch_ID in channels_to_csw:
        interpolated_wfs[ch_ID] = _shift_wf(
            *wu.tgraph_to_arrays(wavepacket['waveforms'][ch_ID]), 
            -arrival_delays[ch_ID], big_times
        )

    # Roll the non-nan part of the waveforms around the preselected reference 
    #   channel & build the CSW
    roll_shifts = _get_shifts(interpolated_wfs, reference_ch)
    rolled_wfs = {ch_ID: np.full(len(big_times), -123456.0) for ch_ID in channels_to_csw}
    csw = np.zeros((1,len(big_times)))
    for c, ch_ID in enumerate(channels_to_csw):
        rolled_wfs[ch_ID] = _roll_wf(interpolated_wfs[ch_ID], roll_shifts[ch_ID], big_times)
        # rolled_wfs[ch_ID] = interpolated_wfs[ch_ID] ## Uncomment to use unrolled wfs in csw
        csw = np.nansum( np.dstack( (csw[0], rolled_wfs[ch_ID]) ), axis=2) 

    # Un-nest the csw. csw.shape was (1,len(big_times)) but is now len(big_times)
    csw = csw[0]

    # Trim the CSW of leading and trailing nan values. Trim time array accordingly.
    idx_non_pad_values = csw != 0.0
    start_idx = np.where( 
        idx_non_pad_values.any(), idx_non_pad_values.argmax(), -123456)
    end_idx = np.where(
        idx_non_pad_values.any(),
        len(csw) - np.flip(idx_non_pad_values).argmax() - 1,
        -123456
    )
    final_csw = csw[ start_idx : end_idx+1 ]
    final_times = big_times[ start_idx : end_idx+1 ]

    # Return requested data
    if return_rolled_wfs: 
        return final_times, final_csw, big_times, rolled_wfs
    else:
        return final_times, final_csw
    
def _trim_array(times, values, trim, trim_method='peak'):

    if trim == 0: 
        return times, values
    
    value_max_idx = np.nanargmax(values)

    front_trim = trim // 2
    back_trim = int( np.ceil( trim / 2 ) )
    if value_max_idx - front_trim < 0:
        # The peak of the signal is within the front region of the waveform 
        #   about to be trimmed. Switch to trimming off the back entirely.
        front_trim = 0
        back_trim = trim
    elif back_trim - value_max_idx < 0:
        # The peak of the signal is within the back region of the waveform
        #   about to be trimmed. Switch to trimming off the front entirely.
        front_trim = trim
        back_trim = 0

    if back_trim == 0:
        return times[front_trim:], values[front_trim:]
    else: 
        return times[front_trim:-back_trim], values[front_trim:-back_trim]

def get_csw(e,
    data : AnalysisDataset, 
    useful_event, 
    solution : int ,
    polarization : int,
    excluded_channels = None,
    which_distance : str = "distant",
    return_rolled_wfs : bool = False
):
    """
    From a given `useful_event`, return the coherently summed waveform (CSW)
    """

    # Pull wavepacket from the dataset
    wavepacket = data.get_wavepacket(
        useful_event=useful_event, which_traces='filtered')

    # Make a list of channels to analyze
    channels_to_csw = _get_channels_to_csw(excluded_channels, data, polarization)

    # In case some channels have different lengths than others, choose the
    #   smallest waveform size for the length of the csw
    csw_length = 123456
    for ch_ID in channels_to_csw:
        if wavepacket['waveforms'][ch_ID].GetN() < csw_length: 
            csw_length = wavepacket['waveforms'][ch_ID].GetN() 

    # Open up the reconstruction for this event and polarization
    reco = sr.StandardReco(data.station_id, excluded_channels=excluded_channels)
    reco_results = reco.do_standard_reco(wavepacket)
    pol_int_to_str = {0: "v", 1: "h"}
    distance_to_key = {'distant': 'distant', 'nearby': 'pulser'}
    reco_key = f"{distance_to_key[which_distance]}_{pol_int_to_str[polarization]}"

    # Uncomment these lines to plot the skymap for this event's reconstruction
    # from araproc.framework import data_visualization as dv
    # dv.plot_skymap(
    #     reco_results[reco_key]["map"],
    #     ouput_file_path=f"./u-station_{data.station_id}_run_{data.run_number}_event_{e}_map.png")

    # Get arrival delays relative to the reference channel (the one with the
    #   highest voltage response)
    reference_ch = _get_channel_with_max_V(wavepacket['waveforms'], channels_to_csw)
    arrival_delays_reco = _get_arrival_delays_reco(
        data, reco, reco_results[reco_key], channels_to_csw, reference_ch, 
        which_distance, solution)
    if useful_event.isSoftwareTrigger(): wf_type = 'soft'
    elif useful_event.isCalpulserEvent(): wf_type = 'cp'
    else: wf_type = 'rf'
    arrival_delays = _get_arrival_delays_AraRoot_xcorr(
        wavepacket['waveforms'], reco, polarization, reference_ch, arrival_delays_reco, wf_type)
    
    # Prepare the final CSW waveform time and voltage arrays
    csw_values = np.zeros((1, csw_length))
    csw_times = _trim_array(
        np.asarray(wavepacket['waveforms'][reference_ch].GetX()), 
        np.asarray(wavepacket['waveforms'][reference_ch].GetY()), 
        wavepacket['waveforms'][reference_ch].GetN()-csw_length
    )[0]
    dt = csw_times[1] - csw_times[0]
     
    # 
    rolled_wfs = {ch_ID: np.full(len(csw_times), -123456.0) for ch_ID in channels_to_csw}
    for ch_ID in channels_to_csw: 

        # Load this channels voltage and time arrays. Shift time by arrival delay
        arrival_delay = arrival_delays[ch_ID]
        values = np.asarray(wavepacket['waveforms'][ch_ID].GetY())
        times = np.asarray(wavepacket['waveforms'][ch_ID].GetX()) - arrival_delay

        # Determine if the time binning of this channel matches the time binning
        #   of the CSW. Ex: If the CSW has times of [0.9, 1.4, 1.9] then a channel
        #   with a time array of [1.6, 2.1, 2.6] needs to be rebinned/shifted 
        #   backwards by 0.2 but a channel with the time array [3.4, 3.9, 4.4 ] 
        #   doesn't because each bin starts and ends on the same fraction of a 
        #   nanosecond that the CSW does. 
        rebinning_shift = round( (times[0] - csw_times[0]) % dt, 4) 
        if rebinning_shift != 0: 
            waveform_model = Akima1DInterpolator(times, values)
            times = times - rebinning_shift
            values = waveform_model(times)

        # Trim the array to match the CSW length but try not to remove the 
        #   maximal waveform point
        # TODO could combine this logic with the rolling
        times, values = _trim_array(times, values, len(values)-csw_length)

        # Now roll the waveform so that the start and end times of the waveform
        #   line up exactly with the CSW
        roll_shift = (times[0] - csw_times[0])/dt
        # Check that roll_shift is close to an integer
        rolled_wfs[ch_ID] = np.roll( values, int(roll_shift) )
        rolled_times = np.linspace(
            times[0] - roll_shift*(times[1] - times[0]),
            times[-1] - roll_shift*(times[1] - times[0]),
            len(times)
        )

        # Add this channel's waveform to the CSW
        csw_values = np.nansum( np.dstack( (csw_values[0], rolled_wfs[ch_ID]) ), axis=2) 
        

    # arrival_delays = _get_arrival_delays_xcorr(
    #     wavepacket, reco, channels_to_csw, reference_ch)
    
    # # Shift each waveform by the arrival delay then interpolate each waveform 
    # #   so it exists within an array that stretches from the earliest wf 
    # #   time to the last wf time. Waveform values are set to `nan` at all 
    # #   times for which it does not exist. These will be trimmed later.
    # big_times = _get_catchall_time_array(wavepacket['waveforms'], channels_to_csw, arrival_delays)
    # interpolated_wfs = {
    #     ch_ID: np.full(len(big_times), -123456.0) for ch_ID in channels_to_csw}
    # for ch_ID in channels_to_csw:
    #     interpolated_wfs[ch_ID] = _shift_wf(
    #         *wu.tgraph_to_arrays(wavepacket['waveforms'][ch_ID]), 
    #         -arrival_delays[ch_ID], big_times
    #     )

    # # Roll the non-nan part of the waveforms around the preselected reference 
    # #   channel & build the CSW
    # roll_shifts = _get_shifts(interpolated_wfs, reference_ch)
    # rolled_wfs = {ch_ID: np.full(len(big_times), -123456.0) for ch_ID in channels_to_csw}
    # csw = np.zeros((1,len(big_times)))
    # for c, ch_ID in enumerate(channels_to_csw):
    #     rolled_wfs[ch_ID] = _roll_wf(interpolated_wfs[ch_ID], roll_shifts[ch_ID], big_times)
    #     # rolled_wfs[ch_ID] = interpolated_wfs[ch_ID] ## Uncomment to use unrolled wfs in csw
    #     csw = np.nansum( np.dstack( (csw[0], rolled_wfs[ch_ID]) ), axis=2) 

    # Un-nest the csw. csw.shape was (1,len(big_times)) but is now len(big_times)
    csw_values = csw_values[0]

    # # Trim the CSW of leading and trailing nan values. Trim time array accordingly.
    # idx_non_pad_values = csw != 0.0
    # start_idx = np.where( 
    #     idx_non_pad_values.any(), idx_non_pad_values.argmax(), -123456)
    # end_idx = np.where(
    #     idx_non_pad_values.any(),
    #     len(csw) - np.flip(idx_non_pad_values).argmax() - 1,
    #     -123456
    # )
    # final_csw = csw[ start_idx : end_idx+1 ]
    # final_times = big_times[ start_idx : end_idx+1 ]

    # Return requested data
    if return_rolled_wfs: 
        return csw_times, csw_values, np.copy(csw_times), rolled_wfs
    else:
        return csw_times, csw_values

def get_csw_snr(csw : np.ndarray):
    """
    Calculate the signal-to-noise ratio of the provided coherently summed 
        waveform (CSW). Quarters the CSW and uses the mean of the two quarters 
        with the lowest RMS as the RMS.

    Parameters
    ----------
    csw : np.ndarray
        Coherently Summed Waveform
    
    Returns
    -------
    csw_snr : float
        The signal-to-noise ratio of the coherently summed waveform. 
    """

    # Calculate peak-to-peak
    p2p = ( np.max(csw) - np.min(csw) )

    # Use the average of the two quarters of the CSW with the lowest RMS
    #   for the RMS value in the SNR calculation
    quartered_csw = np.array_split(csw, 4)
    rms_values = [ np.sqrt(np.mean( quarter**2 )) for quarter in quartered_csw]
    sorted_rms = sorted(rms_values)
    rms = np.mean( sorted_rms[:2] )

    # Calculate SNR
    csw_snr = p2p / 2 / rms

    return csw_snr

def plot_csw(csw_times, csw, plot_dir):
    import matplotlib.pyplot as plt

    # plot CSW
    fig, ax = plt.subplots()
    ax.plot(csw_times, csw)
    ax.set_xlabel("Time")
    ax.set_ylabel("Voltage")
    ax.set_title("CSW")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/CSW.png", dpi=400)

    return

def plot_rolled_wfs_overlayed(
    times, rolled_wfs, save_name, excluded_channels, data, polarization, 
    xlims=None, ylims=None, plot_csw=True):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    # Make the 4x4 of rolled plots
    fig, ax = plt.subplots()
    ax.set_title(f"Waveforms for Entry CSW", size=16)
    ax.set_xlabel("Time", size=14, labelpad=15.0)
    ax.set_ylabel("Voltage", size=14, labelpad=23.0)

    if plot_csw:
        csw = np.zeros((1,len(times)))
        for ch_ID, wf in rolled_wfs.items():
            csw = np.nansum( np.dstack( (csw[0], wf) ), axis=2) 
        ax.fill_between(
            times, np.zeros(len(times)), csw[0], label='CSW',
            color='k', alpha=0.05, lw=0)

    # Create 4x4 subplots and plot 16 waveforms
    antenna_number = 0
    channels_to_csw = _get_channels_to_csw(excluded_channels, data, polarization)

    cmap = plt.get_cmap('gist_rainbow')
    for i, (channel, wfs) in enumerate(rolled_wfs.items()): 
        ax.plot(
            times, wfs, label=channels_to_csw[i], alpha=0.4,
            color=cmap( i/len(rolled_wfs) )
        )
        ax.scatter(*_get_peak(times, wfs), color=cmap( i/len(rolled_wfs)))
    if xlims != None: 
        ax.set_xlim(xlims)
    if ylims != None: 
        ax.set_ylim(ylims)
    ax.legend(loc=1)
            
    # Tighten and save figure
    plt.tight_layout()
    plt.savefig(f"{save_name}", dpi=300)

    return

def plot_rolled_wfs(
    times, rolled_wfs, save_name, excluded_channels, data, polarization, xlims=None, ylims=None):
    import matplotlib.pyplot as plt

    # Make the 4x4 of rolled plots
    fig = plt.figure(figsize=(10,16))
    plt.suptitle(f"Waveforms for Entry CSW", size=16)

    # Create the major subplot so we can have universal x and y labels
    ax_major = fig.add_subplot(111) 
    ax_major.spines['top'].set_color('none') # Turn off 
    ax_major.spines['bottom'].set_color('none')
    ax_major.spines['left'].set_color('none')
    ax_major.spines['right'].set_color('none')
    ax_major.tick_params(
        labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_major.set_xlabel("Time", size=14, labelpad=15.0)
    ax_major.set_ylabel("Voltage", size=14, labelpad=23.0)

    # Create 4x4 subplots and plot 16 waveforms
    gs = fig.add_gridspec(8, 2, hspace=0, wspace=0)
    axs = gs.subplots(sharex='col', sharey='row')
    antenna_number = 0
    channels_to_csw = _get_channels_to_csw(excluded_channels, data, polarization)
    def get_axs_coordinates(antenna_number): 
        """
        Returns
        -------
        horizontal coordinate : int
        vertical coordinate : int
        """
        return antenna_number // 8, antenna_number % 8
    for c, (channel, wfs) in enumerate(rolled_wfs.items()): 
        if antenna_number >= 16: 
            print("Too many antennas to plot")
            break
        axs_horizontal, axs_vertical = get_axs_coordinates(channels_to_csw[c])
        axs[axs_vertical][axs_horizontal].plot(times, wfs, label=channels_to_csw[c])
        axs[axs_vertical][axs_horizontal].legend(loc=2)
        antenna_number += 1
    if xlims != None: 
        for i in range(8):
            for j in range(2):
                axs[i][j].set_xlim(xlims)
    if ylims != None: 
        for i in range(8):
            for j in range(2):
                axs[i][j].set_ylim(ylims)
            
    # Tighten and save figure
    plt.tight_layout()
    plt.savefig(f"{save_name}", dpi=300)

    return

# # A1
# run = "019353"
# data = AnalysisDataset(
#     station_id = 100,
#     path_to_data_file=f"/data/exp/ARA/2020/L1/10pct/ARA01/1030/run{run}/event{run}.root",
#     path_to_pedestal_file=f"/data/ana/ARA/ARA01/ped_full/ped_full_qualities_A1_R{int(run)}.dat",
#     is_simulation = False
# )
# xlims = (290,235)
# event = 7

# # A2
# run = "019290"
# data = AnalysisDataset(
#     station_id = 2,
#     path_to_data_file=f"/data/exp/ARA/2021/unblinded/L1/ARA02/0302/run{run}/event{run}.root",
#     path_to_pedestal_file=f"/data/ana/ARA/ARA02/ped_full/ped_full_values_A2_R{int(run)}.dat",
#     is_simulation = False
# )
# xlims=(90,135)
# event = 5

# # A3
# run = "018686"
# data = AnalysisDataset(
#     station_id = 3,
#     path_to_data_file=f"/data/exp/ARA/2020/L1/10pct/ARA03/0828/run{run}/event{run}.root",
#     path_to_pedestal_file=f"/data/ana/ARA/ARA03/ped_full/ped_full_values_A3_R{int(run)}.dat",
#     is_simulation = False
# )
# xlims=(60,120)
# event = 3

# # A4
# run = "004214"
# data = AnalysisDataset(
#     station_id = 4,
#     path_to_data_file=f"/data/exp/ARA/2018/L1/10pct/ARA04/0530/run{run}/event{run}.root",
#     path_to_pedestal_file=f"/data/ana/ARA/ARA04/ped_full_from_Martin/ped_full_values_A4_run{run}.dat",
#     is_simulation = False
# )
# xlims=(-10,45)
# event = 64

# A5
run = "005626"
data = AnalysisDataset(
    station_id = 5,
    path_to_data_file=f"/data/exp/ARA/2019/blinded/L1/ARA05/0701/run{run}/event{run}.root",
    path_to_pedestal_file=f"/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run{run}.dat",
    is_simulation = False
)
xlims=(125,160)
event = 54

useful_event = data.get_useful_event(event)

csw_times, csw, rolled_wf_times, rolled_wfs = get_csw(
    event,
    data, 
    useful_event,
    0,
    0, 
    which_distance="nearby",
    return_rolled_wfs=True, 
    excluded_channels=data.excluded_channels
    # excluded_channels=np.concatenate((data.excluded_channels,[2,3,4,5,6,7]))
)

print("\nCSW_SNR", round(get_csw_snr(csw),1))
print()
file_id = f"A{data.station_id}_r{run}_e{event}"
# plot_csw(csw_times, csw, "./")
# plot_rolled_wfs(
#     rolled_wf_times, rolled_wfs, f"./u-wfs_for_csw_{file_id}.png", data.excluded_channels, data, 0,
#     ylims=(-800, 1100))
# plot_rolled_wfs_overlayed(
#     rolled_wf_times, rolled_wfs, f"./u-wfs_for_csw_{file_id}_overlayed_with_csw.png", 
#     data.excluded_channels, data, 0, xlims=xlims 
#     )
# plot_rolled_wfs_overlayed(
#     rolled_wf_times, rolled_wfs, f"./u-wfs_for_csw_{file_id}_overlayed.png", 
#     data.excluded_channels, data, 0,
#     # xlims=(125,160), 
#     # xlims=(290, 340),
#     plot_csw=False)