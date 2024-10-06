import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.framework import waveform_utilities as wu

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('mpl.style')

# (eventNumber 213179, which is TTree index 20695 in the burn sample file)
d = dataset.AnalysisDataset(
    station_id = 2,
    path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
    path_to_pedestal_file="/cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/ara_build/share/araCalib/ATRI/araAtriStation2Pedestals.txt"
)
useful_event = d.get_useful_event(20695)
wave_bundle_raw = d.get_waveforms(useful_event, which_traces="interpolated")
wave_bundle_filt = d.get_waveforms(useful_event)


wave_raw_t, wave_raw_v = wu.tgraph_to_arrays(wave_bundle_raw[0])
wave_filt_t, wave_filt_v = wu.tgraph_to_arrays(wave_bundle_filt[0])

wave_raw_freqs, wave_raw_spectrum = wu.time2freq(wave_raw_t, wave_raw_v)
wave_raw_spectrum = np.log10(np.abs(wave_raw_spectrum))

wave_filt_freqs, wave_filt_spectrum = wu.time2freq(wave_filt_t, wave_filt_v)
wave_filt_spectrum = np.log10(np.abs(wave_filt_spectrum))

fig, ax = plt.subplots(1,1)
ax.plot(wave_raw_freqs, wave_raw_spectrum, label="Original")
ax.plot(wave_filt_freqs, wave_filt_spectrum, '--', label="Filtered")
# ax.set_xlim(0.25,0.55)
ax.set_xlabel("Frequency / GHz")
ax.set_ylabel("log10(Spectrum) / au ")
ax.legend()
plt.tight_layout()
fig.savefig("cw_demo.png", dpi=300)
plt.close(fig)
del fig, ax
