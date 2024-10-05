import numpy as np

from araproc.framework import waveform_utilities as wfu

def get_vpp(waveform):

    """
    Calculates the peak-to-peak voltage of a voltage trace.

    Parameters
    ----------
    waveform: TGraph
        A numpy array of floats containing the voltages for the trace.

    Returns
    -------
    vpp : float
        Peak-to-peak voltage in same units as trace.
    """ 

    _, trace = wfu.tgraph_to_arrays(waveform)

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
    waveform: TGraph
        A numpy array of floats containing the voltages for the trace.
    nsegs: int
        Number of segments to break waveform into.

    Returns
    -------
    rms : float
        The minimum RMS among all segments in same units as trace.
    """

    _, trace = wfu.tgraph_to_arrays(waveform)

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
