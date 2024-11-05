import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.framework import waveform_utilities as wu
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from araproc.analysis import standard_reco as sr
from araproc.framework import data_visualization as dv
from araproc.analysis import snr


import logging
logging.getLogger().setLevel(logging.INFO)

# d = dataset.AnalysisDataset(
#     station_id=5,
#     path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
# )
# iter_start = 0
# iter_stop = 2000

d = dataset.AnalysisDataset(
    station_id = 2,
    path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
)
iter_start = 19690
iter_stop = 20900

# instantiate the standard reconstruction methods
reco = sr.StandardReco(d.station_id, 
                       excluded_channels=d.excluded_channels)

dd_vals = []
filt_vals = []
dd_corr_vals = []
filt_corr_vals = []

for e in range(iter_start, iter_stop, 1):

    useful_event = d.get_useful_event(e)
    
    process_event = useful_event.isCalpulserEvent()

    if process_event:
      
        # wavepacket_raw = d.get_wavepacket(useful_event, which_traces="calibrated")
        # wavepacket_interp = d.get_wavepacket(useful_event, which_traces="interpolated")
        wavepacket_dd = d.get_wavepacket(useful_event, which_traces="dedispersed")
        wavepacket_filt = d.get_wavepacket(useful_event, which_traces="filtered")    

        avg_snr_dd = snr.get_avg_snr(wavepacket_dd["waveforms"], excluded_channels=d.excluded_channels_and_hpol)
        avg_snr_filt = snr.get_avg_snr(wavepacket_filt["waveforms"], excluded_channels=d.excluded_channels_and_hpol)
        print(f"The Average SNR dd {avg_snr_dd:.1f}, avg snr filt is {avg_snr_filt:.1f}")
        dd_vals.append(avg_snr_dd)
        filt_vals.append(avg_snr_filt)

        # raw_t, raw_v = wu.tgraph_to_arrays(wavepacket_raw["waveforms"][1])
        # interp_t, interp_v = wu.tgraph_to_arrays(wavepacket_interp["waveforms"][1])
        # dd_t, dd_v = wu.tgraph_to_arrays(wavepacket_dd["waveforms"][1])
        # filt_t, filt_v = wu.tgraph_to_arrays(wavepacket_filt["waveforms"][1])

        # freqs_interp, spectrum_interp = wu.time2freq(interp_t, interp_v)
        # spectrum_interp = np.log10(np.abs(spectrum_interp))

        # freqs_dd, spectrum_dd = wu.time2freq(dd_t, dd_v)
        # spectrum_dd = np.log10(np.abs(spectrum_dd))        

        # freqs_filt, spectrum_filt = wu.time2freq(filt_t, filt_v)
        # spectrum_filt = np.log10(np.abs(spectrum_filt))



        # ## butterworth filter myself

        # # nyquist = 1./(2.*(filt_t[1]-filt_t[0])*1E-9) # interp speed in seconds
        # # freq_hipass = 100E6/nyquist # highpass filter at 90 MHz, in units of nyquist

        # # # b, a = signal.butter(5, freq_hipass*2*np.pi, "highpass", analog=True)
        # # b, a = signal.butter(5, freq_hipass, "highpass")
        # # response = signal.lfilter(b,a,dd_v)
        # # freqs_custm_filt, spectrum_custom_filt = wu.time2freq(dd_t, response)
        # # spectrum_custom_filt = np.log10(np.abs(spectrum_custom_filt))

        # fig, (ax, ax2) = plt.subplots(1,2, figsize=(10,5))
        
        # # ax.plot(raw_t, raw_v)
        # # ax.plot(interp_t, interp_v, label="Interpolated", color="C1")
        # ax.plot(dd_t, dd_v, label="DD", color="C2")
        # ax.plot(filt_t, filt_v, label="Filt", color="C3")
        # # ax.plot(dd_t, response, label="FiltCustom", color="C4", ls='--')
        # # ax.set_xlim([290,330])
        # # ax.set_ylim([-600,600])

        # ax.set_xlim([100,200])
        # ax.set_ylim([-200,200])

        # # ax2.plot(freqs_interp, spectrum_interp, color="C1")
        # ax2.plot(freqs_dd, spectrum_dd, color="C2")
        # ax2.plot(freqs_filt, spectrum_filt, color="C3")
        # # ax2.plot(freqs_custm_filt, spectrum_custom_filt, color="C4", ls='--')
        # ax2.set_xlim([0,0.15])
        
        # fig.savefig(f"event_{e}.png")


        reco_results_dd = reco.do_standard_reco(wavepacket_dd)
        reco_results_filt = reco.do_standard_reco(wavepacket_filt)

        dd_corr_vals.append(reco_results_dd["pulser_v"]["corr"])
        filt_corr_vals.append(reco_results_filt["pulser_v"]["corr"])



        # dv.plot_skymap(reco_results_dd["pulser_v"]["map"],
        #                ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_map_dd.png",
        #                )
        # dv.plot_skymap(reco_results_filt["pulser_v"]["map"],
        #                ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_map_filt.png",
        #                )

np.savez(f"values_{d.station_id}.npz", dd_vals=dd_vals, filt_vals=filt_vals, dd_corr_vals=dd_corr_vals, filt_corr_vals=filt_corr_vals)