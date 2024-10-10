import numpy as np
from scipy.ndimage import uniform_filter1d

class RPRCalculator:
    def __init__(self, dt = 1/3.2, sum_win=25,use_debug = False):
        self.dt = dt
        self.sum_win = sum_win
        self.sum_win_idx = int(np.round(self.sum_win / self.dt))
        self.use_debug = use_debug

    def get_max_info(self):
        """Find the maximum bin, voltage, and time for the current channel."""
        self.max_bin = np.nanargmax(self.pad_v)
        self.max_val = self.pad_v[self.max_bin]
        self.max_time = self.pad_t[self.max_bin]

        if self.use_debug:
            print(f"Max value: {self.max_val}, Max time: {self.max_time}")

    def get_sqrt_volt_sum_wf(self):
        """Apply a rolling mean (smoothing) to the squared waveform."""
        self.pad_v = np.sqrt(uniform_filter1d(self.pad_v, size=self.sum_win_idx, mode='constant'))

    def get_mean_sigma_in_no_max(self):
        """Calculate mean and standard deviation excluding the window around the signal peak."""
        bin_8 = self.pad_num // 8  # Define a window size based on pad_num
        pad_no_max = self.pad_v.copy()  # Create a copy of the waveform
        front_idx = max(0,int(self.max_bin - bin_8))
        
        if front_idx < 0:
            front_idx = 0  # Ensure the index doesn't go negative
        
        pad_no_max[front_idx:self.max_bin + bin_8 + 1] = np.nan  # Exclude the signal window around the peak
        
        self.pad_mean = np.nanmean(pad_no_max)  # Compute mean excluding the peak window
        self.pad_sigma = np.nanstd(pad_no_max)  # Compute standard deviation excluding the peak window

        if self.use_debug:
            print(f"Mean: {self.pad_mean}, Sigma: {self.pad_sigma}")

    def get_ch_sliding_v2_snr_uw(self):
        """Compute the RPR (SNR-like) value."""
        # Smoothing the waveform (apply rolling mean)

        self.get_sqrt_volt_sum_wf()
        # Get max bin, value, and corresponding time
        self.get_max_info()

        # Calculate mean and sigma after excluding the peak window
        self.get_mean_sigma_in_no_max()

        # Calculate SNR-like (RPR) value
        self.rpr_arr = (self.max_val - self.pad_mean) / self.pad_sigma
        
        # Handle channels where sigma is non-positive
        nega_sigma_idx = self.pad_sigma <= 0
        if nega_sigma_idx:
           print(f'negative std detected, which may be unphysical')
           self.rpr_arr = -9999 ## Assign some unphysical value
        
        if self.use_debug:
            print(f"RPR: {self.rpr_arr}, Max Val: {self.max_val}, Mean: {self.pad_mean}, Sigma: {self.pad_sigma}")

    def run_rpr_calculation(self, pad_v=None, pad_t=None, pad_num=None, snr=None):
        """Prepare input data, then compute the RPR value."""

        self.pad_v = pad_v ** 2  # Square the voltage (squaring for SNR-like calculation)
        self.pad_v[np.isnan(self.pad_v)] = 0  # Replace NaNs with zeros
        self.pad_t = pad_t
        self.pad_num = pad_num

        # Compute RPR (SNR-like) value
        self.get_ch_sliding_v2_snr_uw()

        if not self.use_debug:
           # Clear unnecessary variables if not debugging
           del self.pad_v, self.pad_t, self.pad_num, self.max_bin, self.max_val, self.max_time, self.pad_mean, self.pad_sigma



