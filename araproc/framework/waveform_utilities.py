from array import array
import numpy as np
import ROOT

def tgraph_to_arrays(the_graph):
    xvals = np.asarray(the_graph.GetX())
    yvals = np.asarray(the_graph.GetY())
    return xvals, yvals

def arrays_to_tgraph(xvals, yvals):
    if len(xvals) != len(yvals):
        raise Exception("Lenfth of x and y arrays is mismatched. Abort.")
    return ROOT.TGraph(len(xvals), xvals, yvals)

def get_freqs(times):
    # really simple function to get the frequencies
    nsamples = len(times)
    dt = times[1] - times[0]
    return np.fft.rfftfreq(nsamples, dt)

def get_dt_and_sampling_rate(times):
    dt = times[1] - times[0]
    sampling_rate = 1/dt
    return dt, sampling_rate

def time2freq(times, volts):
    """
    Performs forward FFT with correct normalization that conserves the power
    Note the returned values are *complex*

    Parameters
    ----------
    times: np.array
        time samples
    volts: np.array
        voltage samples to be FFT'd

    Returns
    -------
    freqs: np.array
        the frequencies at which the fft was evaluated
    the_fft
        the complex FFT of the voltages
        I have named it the_fft specifically to discourage people accidentally
        misinterpreting this as the spectrum.
        To get the traditional power spectrum, you should do
        np.abs(the_fft).

    """
    
    dt, sampling_rate = get_dt_and_sampling_rate(times)

    # an additional sqrt(2) is added because negative frequencies are omitted
    the_fft = np.fft.rfft(volts, axis=-1) / sampling_rate * 2 ** 0.5
    
    freqs = get_freqs(times) # get the frequencies
    
    return freqs, the_fft

def freq2time(times, spectrum):
    """
    Performs backward FFT with correct normalization that conserves the power.
    Please note that it is the user's job to supply the time axis.
    Python can't reconstruct that for you.


    Parameters
    ----------
    times: np.array
        time samples
    spectrum: complex np array
        the complex numpy array containing the fft

        
    Returns
    -------
    times: np.array
        time samples
    volts: np.array
        voltage samples after the ifft

    """
    
    dt, sampling_rate = get_dt_and_sampling_rate(times)
    
    # you need to know the number of samples in the time domain 
    # (relevant if time trace has an odd number of samples)
    n_samples_t_domain = len(times)

    volts = np.fft.irfft(spectrum, axis=-1, n=n_samples_t_domain) * sampling_rate / 2 ** 0.5
    return times, volts
