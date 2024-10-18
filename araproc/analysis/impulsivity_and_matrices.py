import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.special import erf

# Define the erf-linear model for curve fitting
def erf_linear(x, A, B):
    return (A * erf(x / B) + x) / (A * erf(1 / B) + 1)

# Generate a Gaussian noise waveform (600 ns)
np.random.seed(42)
sampling_interval = 0.5  # ns
time_array = np.arange(0, 600, sampling_interval)
voltage_array = np.random.normal(0, 1, len(time_array))

# Function to calculate impulsivity-related variables
def calculate_impulsivity_measures(voltage_array,time_array):
    result = {}

    # Hilbert transform to get the envelope of the waveform
    hilbert_envelope = np.abs(hilbert(voltage_array))

    # Find the index of the maximum value in the Hilbert envelope
    hill_max_idx = np.argmax(hilbert_envelope)
    hill_max = hilbert_envelope[hill_max_idx]

    # Sorting based on closeness to the maximum index
    closeness = np.abs(np.arange(len(voltage_array)) - hill_max_idx)
    clo_sort_idx = np.argsort(closeness)

    # Sort the Hilbert envelope by closeness to the maximum
    sorted_waveform = hilbert_envelope[clo_sort_idx]

    # Cumulative distribution function (CDF) calculation
    cdf = np.cumsum(sorted_waveform)
    cdf /= np.max(cdf)

    # Linear regression to get slope, intercept, and other statistics
    slope, intercept, r_value, p_value, std_err = linregress(np.arange(len(cdf)), cdf)
    cdf_fit = slope * np.arange(len(cdf)) + intercept
    t_frac = np.linspace(0,1,len(cdf))

    # Calculate Kolmogorov-Smirnov statistic 
    ks = np.max(np.abs(cdf_fit - cdf))
    # Perform erf-linear fit on the CDF
    popt, _ = curve_fit(erf_linear, t_frac, cdf, p0=[intercept / slope, 0.01], bounds=(0, [3 * intercept / slope, 0.5]))
    A_fit, B_fit = popt
    cdf_erf_fit = erf_linear(t_frac, A_fit, B_fit)

    # Calculate impulsivity
    impulsivity = 2 * np.mean(cdf) -1 # np.sum(cdf[1:] * (np.arange(1, len(cdf)) - np.arange(0, len(cdf) - 1))) - 1
    print(impulsivity)
    # Calculate linear chi2 (difference between linear fit and ideal x=y line)
    chi2_linear = np.sum((cdf - np.arange(len(cdf)) / len(cdf)) ** 2)

    # Calculate erf-linear chi2 (difference between erf-linear fit and data)
    chi2_erf_linear = np.sum((cdf - cdf_erf_fit) ** 2)

    # Store the results
    result['impulsivity'] = impulsivity
    result['hill_max'] = hill_max
    result['slope'] = slope
    result['intercept'] = intercept
    result['ks'] = ks  # Kolmogorov-Smirnov statistic
    result['r_value'] = r_value
    result['p_value'] = p_value
    result['std_err'] = std_err
    result['impLinChi2'] = chi2_linear
    result['impErfLinChi2'] = chi2_erf_linear
    result['impErfA'] = A_fit
    result['impErfB'] = B_fit

    # Plotting all results at the end

    ################ These plotting scripts are for debugging #################
    ################# Once  this function pass the reviews abd corrections, these plot scripts will be removed #########
    plt.figure(figsize=(15, 10))
    plt.subplot(3, 1, 1)
    plt.plot(time_array, voltage_array, label='Waveform')
    plt.plot(time_array, hilbert_envelope, label='Hilbert Envelope')
    plt.title('Waveform and Hilbert Envelope')
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude')
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(np.arange(len(cdf)), cdf, label='CDF')
    plt.plot(np.arange(len(cdf)), cdf_fit, label='Linear Fit', linestyle='--')
    plt.axhline(np.mean(cdf),label = 'mean cdf', c  = 'black', linestyle='--')
    plt.title('CDF with Linear Fit')
    plt.xlabel('Samples')
    plt.ylabel('CDF Value')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(np.arange(len(cdf)), cdf, label='CDF')
    plt.plot(np.arange(len(cdf)), cdf_erf_fit, label='Erf-Linear Fit', linestyle='--')
    plt.axhline(np.mean(cdf),label = 'mean cdf', c  = 'black', linestyle='--')
    plt.title('CDF with Erf-Linear Fit')
    plt.xlabel('Samples')
    plt.ylabel('CDF Value')
    plt.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

    return result

# Call the function and print results

