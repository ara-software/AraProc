import ROOT
import numpy as np

from araproc.framework import waveform_utilities as wfu

def get_vpp(waveform, use_local = True, time_window = 20):

    """ 
    Calculates the peak-to-peak voltage of a signal using local (Default) and global extrema.

    Parameters
    ----------
    waveform: ROOT.TGraph, tuple, or list
        A TGraph or a tuple/list containing two np.ndarrays: (time, trace).
    use_local: bool, optional (Default = True)
        If True, returns the peak-to-peak voltage based on local extrema. 
        If False, returns the peak-to-peak voltage based on global extrema.
    time_window: float, optional (Default = 20 ns)
        Time window in nanoseconds to calculate the peak-to-peak voltage.

    Returns
    -------
    vpp : float
        Peak-to-peak voltage in same units as trace.
    """

    if(isinstance(waveform, ROOT.TGraph)):
      time, trace = wfu.tgraph_to_arrays(waveform)
      if len(trace) != len(time):
        raise ValueError("Trace and time must have the same length.")

    elif isinstance(waveform, (list, tuple)):
      if len(waveform) != 2:
        raise ValueError("Waveform tuple or list must have length 2!")
      time, trace = waveform
      if (not isinstance(time, np.ndarray)) or (not isinstance(trace, np.ndarray)):
        raise ValueError("Waveform must be tuple or list of np.ndarrays!")        
      if len(trace) != len(time):
        raise ValueError("Trace and time must have the same length.")

    else:
      raise Exception("Unsupported data type in snr.get_vpp. Abort")

    # global extrema case
    if not use_local:
      vpp = trace.max() - trace.min()
      return vpp

    # if there aren't at least two extrema
    if trace.max() == trace.min():
      return 0.0

    # find points at least half as large as the global min/max
    upper_peak_idx = trace > 0.5*trace.max() 
    lower_peak_idx = trace < 0.5*trace.min() 

    # limit to just these points
    upper_peak_times = time[upper_peak_idx]
    upper_peak_voltages = trace[upper_peak_idx]
    
    lower_peak_times = time[lower_peak_idx]
    lower_peak_voltages = trace[lower_peak_idx]

    # find differences between all extrema pairs
    all_dt = upper_peak_times[:, np.newaxis] - lower_peak_times 
    all_vpp = upper_peak_voltages[:, np.newaxis] - lower_peak_voltages 

    # limit to extrema within time_window of each other
    mask = np.abs(all_dt) <= time_window
    if not mask.any():
      return 0.0
    
    valid_vpp = all_vpp[mask]

    # finally select the max vpp within a time_window
    vpp = np.max(valid_vpp)

    return vpp

def get_min_segmented_rms(waveform, nsegs=8):

    """
    Calculates the RMS of a voltage trace by calculating the RMS of segments
    voltage trace and returning the minimum value.

    Parameters
    ----------
    waveform: TGraph or np.ndarray
        A TGraph or np.ndarray of the waveform voltage.
    nsegs: int
        Number of segments to break waveform into.

    Returns
    -------
    rms : float
        The minimum RMS among all segments in same units as trace.
    """

    if(isinstance(waveform, ROOT.TGraph)):
      _, trace = wfu.tgraph_to_arrays(waveform)
    elif(isinstance(waveform, np.ndarray)):
      trace = np.copy(waveform)
      
      # don't be bothered by some unused dimensions
      trace = np.squeeze(trace)

      if(trace.ndim != 1):
        raise Exception("Trace is not 1d in snr.get_min_segmented_rms. Abort")

    else:
      raise Exception("Unsupported data type in snr.get_min_segmented_rms. Abort")

    rms = 1e100
    traceLen = len(trace)
    segLen = traceLen // nsegs
    if(segLen < 2):
      raise Exception("Number of segments cannot be more than number of points in trace. Abort.")

    segRem = traceLen % nsegs

    for i in range(nsegs):

      start = i*segLen
      if(i < segRem):
        start += i

      end = start + segLen
      if(i < segRem):
        end += 1

      thisRms = np.sqrt(np.mean(trace[start:end+1]**2))

      if(thisRms < rms):
        rms = thisRms

    return rms

def get_jackknife_rms(waveform, nSamples=50):

    """
    Calculates the RMS of a voltage trace via a sort-of jackknife technique,
    which attempts to identify the parts of the trace which do not contain 
    signal with the goal of using as much of the trace to estimate the RMS 
    as possible. The RMS is then calculated using all parts of the trace which
    are determined to not be contain signal (ie outlier voltages). The trace
    is segmented into equal-sample segments to avoid sensitivity to the trace
    length.

    Parameters
    ----------
    waveform: TGraph or np.ndarray
        A TGraph or np.ndarray of the waveform voltage.
    nSamples: int
        Number of samples for each segment.

    Returns
    -------
    rms : float
        The jackknifed RMS in same units as trace.
    """
    
    if(isinstance(waveform, ROOT.TGraph)):
      _, trace = wfu.tgraph_to_arrays(waveform)
    elif(isinstance(waveform, np.ndarray)):
      trace = np.copy(waveform)
      
      # don't be bothered by some unused dimensions
      trace = np.squeeze(trace)

      if(trace.ndim != 1):
        raise Exception("Trace is not 1d in snr.get_jackknife_rms. Abort")

    else:
      raise Exception("Unsupported data type in snr.get_jackknife_rms. Abort")

    # first split the trace into segments of nSamples 
    # then calculate the RMS on the trace _without_ that segment (ie the "all-but RMS")
    allRms = []
    traceLen = len(trace)
    trace2 = np.square(trace)
    trace2CumSum = np.cumsum(trace2)
    trace2Sum = trace2CumSum[-1] # sum of squares of full trace
    
    start_idx = np.arange(0, traceLen, nSamples)
    end_idx = np.clip(start_idx + nSamples, 0, traceLen)
    segLens = end_idx - start_idx     
    nSegs = len(segLens)

    segTrace2Sum = trace2CumSum[end_idx-1] - trace2CumSum[start_idx-1]*(start_idx > 0) # sum of squares in each segment (note: if start == 0, we don't subtract anything)
    allRms = np.sqrt((trace2Sum-segTrace2Sum)/(traceLen - segLens)) # all-but rms is calculated from sum of full-trace squares _minus_ sum of segment-only squares
 
    # segments that have outlier (ie nonthermal) voltages, will have
    # an all-but RMS that is significantly smaller than others
    # here we attempt to identify those outlier segments by detecting
    # all-but RMS values which are smaller than the mean of the other
    # all-but RMS values by more the 1 standard deviation of the other all-but RMS values 
    # non-outlier segments are collected to calculate the final RMS
    allRms = np.asarray(allRms)
    subMean = (allRms.sum()-allRms)/(nSegs-1.) # mean of all other all-but RMS values
    subStd = np.sqrt((np.sum(np.square(allRms-subMean)) - np.square(allRms-subMean))/(nSegs-1.)) # std of all other all-but RMS values
    
    mask = (allRms >  subMean - subStd) # if all-but RMS isn't an outlier mark this segment to be included
    mask = np.repeat(mask, nSamples)[:traceLen] # convert the outlier mask array from being per-segment to per-sample
    
    # calculate RMS from non-outlier segments
    rms = np.sqrt(np.mean(trace2[mask]))   
 
    return rms

def get_snr(waveform):

    """
    Calculates SNR of a single voltage trace.

    Parameters
    ----------
    waveform: TGraph or np.ndarray
        A TGraph or np.ndarray of the waveform voltage.

    Returns
    -------
    snr : float
        The SNR of the waveform.
    """

    vpp = get_vpp(waveform)
    rms = get_jackknife_rms(waveform)
    if(rms == 0.0):
      return 0

    snr = vpp/rms/2.0

    return snr

def collect_snrs(wave_bundle, excluded_channels=[]):

    """
    Collects the SNRs of each included channel.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    snrs : list of floats 
        The snrs of all included channels
    """
    
    chans = list(wave_bundle.keys())

    snrs = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      waveform = wave_bundle[chan]
      snr = get_snr(waveform)
      snrs.append(snr)

    return snrs

def get_avg_snr(wave_bundle, excluded_channels=[]):

    """
    Calculates channel-wise averaged SNR.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    avg_snr : float
        The average SNR.
    """

    snrs = collect_snrs(wave_bundle, excluded_channels)
    avg_snr = np.mean(snrs)

    return avg_snr

def get_third_highest_snr(wave_bundle, excluded_channels=[]):

    """
    Calculates the 3rd highest SNR.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    third_snr : float
        The third highest SNR.
    """

    snrs = collect_snrs(wave_bundle, excluded_channels)

    # check that we have enough snrs to have a third
    if len(snrs) < 3:
      raise Exception("Not enough channels included to have a 3rd highest SNR. Abort.")

    sorted_snrs = np.sort(snrs) # sort ascending
    sorted_snrs = sorted_snrs[::-1] # switch to descending order

    third_snr = sorted_snrs[2] # get the 3rd highest 

    return third_snr
    
def get_snr_ratio(wave_bundle, excluded_channels=[]):

    """
    Calculates the ratio of the 1st-to-2nd highest SNR.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    snr_ratio : float
        The SNR ratio.
    """

    snrs = collect_snrs(wave_bundle, excluded_channels)

    sorted_snrs = np.sort(snrs) # sort ascending
    sorted_snrs = sorted_snrs[::-1] # switch to descending order

    snr_ratio = sorted_snrs[0] / sorted_snrs[1] # get the 3rd highest 

    return snr_ratio 

def get_avg_rms(wave_bundle, excluded_channels=[]):

    """
    Calculates channel-wise averaged RMS.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    avg_rms : float
        The average RMS.
    """
    
    chans = list(wave_bundle.keys())

    rms = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      waveform = wave_bundle[chan]
      thisRms = get_jackknife_rms(waveform)
      rms.append(thisRms)

    avg_rms = np.mean(rms)

    return avg_rms


