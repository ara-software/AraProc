import ROOT
import numpy as np

from araproc.framework import waveform_utilities as wfu

def get_vpp(waveform):

    """
    Calculates the peak-to-peak voltage of a voltage trace.

    Parameters
    ----------
    waveform: TGraph or np.ndarray
        A TGraphi or np.ndarray of the waveform voltage.

    Returns
    -------
    vpp : float
        Peak-to-peak voltage in same units as trace.
    """ 

    if(isinstance(waveform, ROOT.TGraph)):
      _, trace = wfu.tgraph_to_arrays(waveform)
    elif(isinstance(waveform, np.ndarray)):
      trace = waveform
    else:
      raise Exception("Unsupported data type in snr.get_vpp. Abort")

    vMax = trace.max()
    vMin = trace.min()
    vpp = vMax - vMin

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
      trace = waveform
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
    rms = get_min_segmented_rms(waveform)
    if(rms == 0.0):
      return 0

    snr = vpp/rms/2.0

    return snr

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

    chans = list(wave_bundle.keys())

    avg_snr = []
    for chan in chans:
      if(chan in excluded_channels):
        continue

      waveform = wave_bundle[chan]
      snr = get_snr(waveform)
      avg_snr.append(snr)

    avg_snr = np.mean(avg_snr)

    return avg_snr
