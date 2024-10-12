import ctypes

def get_corr_map_peak(the_map  = None):
        
        """
        Little utility helper function

        Parameters
        ----------
        the_map : ROOT.TH2D
            A 2D map you want the peak found
        
        
        Returns
        peak_corr : float
            The peak in the correlation  map
        peak_phi : int
            The theta value of the peak of the correlation map
        peak_phi : int
            The phi value of the peak of the correlation map
        """
        
        if the_map is None:
            raise Exception("No TH2D map was passed")
    
        _peakZ = ctypes.c_int()
        _peakTheta = ctypes.c_int()
        _peakPhi = ctypes.c_int()
        the_map.GetMaximumBin(_peakPhi, _peakTheta, _peakZ)

        peak_corr = the_map.GetMaximum()
        peak_phi = the_map.GetXaxis().GetBinCenter(_peakPhi.value)
        peak_theta = the_map.GetYaxis().GetBinCenter(_peakTheta.value)
        return peak_corr, peak_phi, peak_theta
