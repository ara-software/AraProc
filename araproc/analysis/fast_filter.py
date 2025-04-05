import numpy as np
import ROOT

from araproc.framework import waveform_utilities as wfu
from araproc.framework import constants as const
from araproc.analysis import snr

def get_peak_time(waveform):

    """
    Finds the time of the absolute maximum of the trace.

    Parameters
    ----------
    waveform: ROOT.TGraph, tuple, or list
        A TGraph or a tuple/list containing two np.ndarrays: (time, trace).
  
    Returns
    -------
    tPeak : float
        Time of the peak.
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
      raise Exception("Unsupported data type in fast_filter.get_peak_time. Abort")

    tPeak, _ = snr.get_windowed_vpp(time, trace)

    return tPeak

def get_wavefront_rms(wave_bundle, station_id, excluded_channels=[]):

    """
    Calculates the RMS of the real part of channel pair-wise calculated cosZenith.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    station_id : int
        Station id
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    wavefront_rms :  float
        The wavefront rms
    """

    all_chans = list(wave_bundle.keys())
    
    geom_tool = ROOT.AraGeomTool.Instance()
    ROOT.SetOwnership(geom_tool, False)
    station_info = geom_tool.getStationInfo(station_id)

    chan_pos = {}
    chan_peak = {}
    for chan in all_chans:
        chan_pos[chan] = np.asarray(station_info.getAntennaInfo(chan).antLocation) 
        chan_peak[chan] = get_peak_time(wave_bundle[chan]) 

    # find cosTheta for each included channel pair
    cosTheta = []
    for i in range(len(all_chans)):
        ch1 = all_chans[i]
        if ch1 in excluded_channels:
            continue

        for j in range(i+1, len(all_chans)):
            ch2 = all_chans[j]
            if ch2 in excluded_channels:
                continue

            dt = chan_peak[ch1] - chan_peak[ch2]
            dr = chan_pos[ch1] - chan_pos[ch2]
            rho = np.sqrt(dr[0]*dr[0] + dr[1]*dr[1])
            dz = dr[2]
            l = np.sqrt(rho*rho + dz*dz)
        
            zavg = (chan_pos[ch1][2]+chan_pos[ch2][2])/2.
            n = const.get_index_of_refraction(zavg)
            v = const.speed_of_light/n # m/ns

            beta = np.arctan2(rho, dz) # zenith angle of line connecting chan2 to chan1
            # calculate real part of cosTheta
            thisCosTheta = -np.cos(beta) * v*dt/l 
            if(abs(v*dt/l) < 1):
                thisCosTheta += np.sin(beta)*np.sqrt(1.-(v*dt/l)**2)
            cosTheta.append(thisCosTheta)
    cosTheta = np.asarray(cosTheta)

    # calculate rms
    wavefront_rms = np.mean(cosTheta**2) - np.mean(cosTheta)**2

    return wavefront_rms

def get_avg_spacetime_interval(wave_bundle, station_id, excluded_channels=[]):

    """
    Calculates the average spacetime interval, s2 = dr^2 - v^2*dt^2, between peaks in channel-pairs.
    Positive values indicate a causal relationship between the signals in each channel.

    Parameters
    ----------
    wave_bundle: dict of TGraphs or np.ndarrays
        Dictionary of waveform TGraphs or np.ndarrays to be averaged.
    station_id : int
        Station id
    excluded_channels: list
        List of dictionary keys to exclude from average.

    Returns
    -------
    avg_s2 :  float
        The average spacetime interval (since we are using a (-,+,+,+) signature this is >= 0 for causal signals)
    """

    all_chans = list(wave_bundle.keys())

    geom_tool = ROOT.AraGeomTool.Instance()
    ROOT.SetOwnership(geom_tool, False)
    station_info = geom_tool.getStationInfo(station_id)

    chan_pos = {}
    chan_peak = {}
    for chan in all_chans:
        chan_pos[chan] = np.asarray(station_info.getAntennaInfo(chan).antLocation) 
        chan_peak[chan] = get_peak_time(wave_bundle[chan]) 

    # find s2 for each included channel pair
    s2 = []
    for i in range(len(all_chans)):
        ch1 = all_chans[i]
        if ch1 in excluded_channels:
            continue

        for j in range(i+1, len(all_chans)):
            ch2 = all_chans[j]
            if ch2 in excluded_channels:
                continue

            dt = chan_peak[ch1] - chan_peak[ch2]
            dr = chan_pos[ch1] - chan_pos[ch2]
            l = np.sqrt(np.sum(dr**2))
        
            zavg = (chan_pos[ch1][2]+chan_pos[ch2][2])/2.
            n = const.get_index_of_refraction(zavg)
            v = const.speed_of_light/n # m/ns
           
            s2.append(l**2 - v**2*dt**2)

    s2 = np.asarray(s2)
    avg_s2 = np.mean(s2)

    return avg_s2 




