import numpy as np 
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import chi2
from araproc.framework import waveform_utilities as wfu

# Define the erf-linear model for curve fitting
def erf_linear(x, A, B):
    """
    Defines the erf-linear model for curve fitting.

    Parameters
    ----------
    x : np.ndarray
        Independent variable, typically the time fraction or index.
    A : float
        Amplitude scaling factor for the error function.
    B : float
        Scaling factor controlling the width of the error function.

    Returns
    -------
    np.ndarray
        The result of applying the erf-linear model to `x`.
    """
    return (A * erf(x / B) + x) / (A * erf(1 / B) + 1)


# Function to calculate impulsivity-related variables
def calculate_impulsivity_measures(channel_wf,channel_time):
    """
    Calculates impulsivity and other statistical measures from a waveform's voltage and time arrays.

    Parameters
    ----------
    channel_wf : np.ndarray
        Array containing the voltage values of the waveform.
    channel_time : np.ndarray
        Array containing the corresponding time values of the waveform.

    Returns
    -------
    result : dict
        A dictionary containing impulsivity measures and various fit statistics:
        - 'impulsivity' : float
            The calculated impulsivity of the waveform.
        - 'slope' : float
            Slope of the linear regression on the CDF.
        - 'intercept' : float
            Intercept of the linear regression on the CDF.
        - 'ks' : float
            Kolmogorov-Smirnov statistic for the difference between fitted and actual CDF.
        - 'r_value' : float
            Correlation coefficient of the linear fit.
        - 'p_value' : float
            p-value of the linear regression.
        - 'std_err' : float
            Standard error of the regression slope.
        - 'impLinChi2' : float
            Chi-square value for the linear fit.
        - 'impErfLinChi2' : float
            Chi-square value for the erf-linear fit.
        - 'impSig' : float
            Significance of erf-linear fit over the linear fit.
        - 'impErfA' : float
            Fitted parameter A for the erf-linear model.
        - 'impErfB' : float
            Fitted parameter B for the erf-linear model.
    """
    result = {}

    # Hilbert transform to get the envelope of the waveform
    hilbert_envelope = wfu.get_hilbert_envelope(channel_wf)   
    # Find the index of the maximum value in the Hilbert envelope
    hill_max_idx = np.argmax(hilbert_envelope)
    hill_max = hilbert_envelope[hill_max_idx]

    # Sorting based on closeness to the maximum index
    closeness = np.abs(np.arange(len(channel_wf)) - hill_max_idx)
    clo_sort_idx = np.argsort(closeness)

    # Sort the Hilbert envelope by closeness to the maximum
    sorted_waveform = hilbert_envelope[clo_sort_idx]

    # Cumulative distribution function (CDF) calculation
    cdf = np.cumsum(sorted_waveform)
    cdf /= np.max(cdf)

    # Linear regression to get slope, intercept, and other statistics
    slope, intercept, r_value, p_value, std_err = linregress(np.arange(len(cdf)), cdf)
    cdf_fit = slope * np.arange(len(cdf)) + intercept
    t_frac = np.linspace(0, 1, len(cdf))

    # Calculate Kolmogorov-Smirnov statistic
    ks = np.max(np.abs(cdf_fit - cdf))

    # Perform erf-linear fit on the CDF
    popt, _ = curve_fit(erf_linear, t_frac, cdf, p0=[intercept / slope, 1e-2], bounds=([0, 1e-6], [3 * intercept / slope, 0.5]))
    A_fit, B_fit = popt
    cdf_erf_fit = erf_linear(t_frac, A_fit, B_fit)

    # Calculate impulsivity
    impulsivity = 2 * np.mean(cdf) - 1

    # Calculate linear chi2 (difference between linear fit and ideal x=y line)
    chi2_linear = np.sum((cdf - cdf_fit) ** 2)

    # Calculate erf-linear chi2 (difference between erf-linear fit and data)
    chi2_erf_linear = np.sum((cdf - cdf_erf_fit) ** 2)

    # Calculate significance of erf-linear fit over linear fit using Wilks' theorem
    dChi2 = chi2_linear - chi2_erf_linear
    impSig = np.sign(dChi2)*np.sqrt(chi2.ppf(chi2.cdf(abs(dChi2), 2), 1))
    impSig = min(10, impSig) # bound the significance

    # Store the results
    result['impulsivity'] = impulsivity
    result['slope'] = slope
    result['intercept'] = intercept
    result['ks'] = ks  # Kolmogorov-Smirnov statistic
    result['r_value'] = r_value
    result['p_value'] = p_value
    result['std_err'] = std_err
    result['impLinChi2'] = chi2_linear
    result['impErfLinChi2'] = chi2_erf_linear
    result['impSig'] = impSig
    result['impErfA'] = A_fit
    result['impErfB'] = B_fit

    return result

