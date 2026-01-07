import math
import numpy as np
import ROOT

from araproc.framework import constants as const


class AraVertexReco:

    """
    AraVertex-based vertex reconstruction wrapper.
    Assumes: wavepacket["waveforms"] contains calibrated, interpolated TGraph objects
    """

    def __init__(self, station_id: int, excluded_channels = np.array([])):

        if station_id not in const.valid_station_ids:
            raise KeyError(f"Station {station_id} is not supported")

        if not isinstance(excluded_channels, np.ndarray) or excluded_channels.ndim != 1:
            raise ValueError("excluded_channels must be a 1D numpy array")

        self.station_id = station_id
        self.num_channels = const.num_rf_channels
        self.excluded_channels = [int(ch) for ch in excluded_channels]

        # load ARA geomtool
        self.araGeom = ROOT.AraGeomTool.Instance()
        ROOT.SetOwnership(self.araGeom, True)

        # Compute station COG (center of gravity) or station center
        antenna_average = np.zeros(3)
        for i in range(16): # 16 RF channels
            for ii in range(3): # 3 axes (x, y, z)
                antenna_average[ii] += (self.araGeom.getStationInfo(station_id).getAntennaInfo(i).antLocation[ii])
        antenna_average /= 16.0

        # load AraVertex + handler
        self.Reco = ROOT.AraVertex()
        self.Reco.SetCOG(*antenna_average)

        self.RecoHandler = ROOT.AraRecoHandler()
        self.chanLocation = self.RecoHandler.getVectorOfChanLocations(self.araGeom, station_id)

        # supported polarizations
        self.reco_pol = {"aravertex_v", "aravertex_h"}

    def do_aravertex_reco(self, wavepacket, snr_threshold = 5.0):

        """
        Run AraVertex reconstruction.

        Parameters
        ----------
        wavepacket : dict{
                "event": int,
                "waveforms": {ch : TGraph},
                "trace_type": str}

        snr_threshold : float
            Hit-finding SNR threshold

        Returns
        -------
        reco_results : dict{
              "aravertex_v": {...},
              "aravertex_h": {...}}
        """

        reco_results = {}
        waveform_bundle = wavepacket["waveforms"]

        for pol_key, pol in [("aravertex_v", 0), ("aravertex_h", 1)]:

            if pol_key not in self.reco_pol:
                raise RuntimeError(f"Unsupported polarization key: {pol_key}")

            # reset internal AraVertex state
            self.Reco.clear()

            # build waveform vector in channel order
            waveforms = []
            for ch in range(self.num_channels):
                if ch not in waveform_bundle:
                    raise KeyError(f"Missing waveform for channel {ch}")
                waveforms.append(waveform_bundle[ch])

            # polarization-specific channel exclusions
            if pol == 0:
                # VPol reco: exclude HPol channels
                excluded_channels_pol = self.excluded_channels + const.hpol_channel_ids
            else:
                # HPol reco: exclude VPol channels
                excluded_channels_pol = self.excluded_channels + const.vpol_channel_ids

            # identify hits
            self.RecoHandler.identifyHitsPrepToVertex(self.chanLocation, self.Reco, self.station_id, pol, excluded_channels_pol, waveforms, snr_threshold)

            # run vertexing
            reco_out = self.Reco.doPairFitSpherical()

            # handle failures -- set to an unphysical value that won't be cut by the zenith cut
            if not (math.isfinite(reco_out.theta) and math.isfinite(reco_out.phi)):
                reco_results[pol_key] = {
                    "valid": False,
                    "theta": -1234.,
                    "phi": -1234.,
                    "R": np.nan}
                continue

            # convert to ARA conventions
            theta = 90.0 - reco_out.theta * ROOT.TMath.RadToDeg() # in degree
            phi = reco_out.phi * ROOT.TMath.RadToDeg() # in degree
            R = reco_out.R # in meter

            # check whether the "not enough hits" value was returned by AraVertex -- set to an unphysical value that won't be cut by the zenith cut
            # "not enough hits" value is theta = 90.0 (elevation) and phi = 0.0
            if abs(theta - 90.0) < 1e-3 and phi == 0.0:
                reco_results[pol_key] = {
                    "valid": False,
                    "theta": -1234.,
                    "phi": -1234.,
                    "R": np.nan}
                continue
                

            reco_results[pol_key] = {
                "valid": True,
                "theta": theta,
                "phi": phi,
                "R": R}

        return reco_results

