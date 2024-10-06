import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.analysis import standard_reco as sr
from araproc.framework import data_visualization as dv
from araproc.framework import waveform_utilities as wu

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('mpl.style')

d = dataset.AnalysisDataset(
    station_id=5,
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)

reco = sr.StandardReco(d.station_id)

useful_event = d.get_useful_event(54)

wave_bundle_raw = d.get_waveforms(useful_event, which_traces="interpolated")
wave_bundle_filt = d.get_waveforms(useful_event, which_traces="filtered")
reco_results = reco.do_standard_reco(wave_bundle_filt)

wave_raw_t, wave_raw_v = wu.tgraph_to_arrays(wave_bundle_raw[0])
wave_filt_t, wave_filt_v = wu.tgraph_to_arrays(wave_bundle_filt[0])

fig, ax = plt.subplots(1,1)
ax.plot(wave_raw_t, wave_raw_v, label="Original", lw=3)
ax.plot(wave_filt_t, wave_filt_v, '--', label="Filtered", lw=3)
ax.set_xlim(280,340)
ax.set_xlabel("Time / ns ")
ax.set_ylabel("Voltage / Volt ")
ax.legend()
plt.tight_layout()
fig.savefig("pulser_demo.png", dpi=300)
plt.close(fig)
del fig, ax


# plot the waveforms
dv.plot_waveform_bundle(wave_bundle_filt, 
            time_or_freq="time",
            ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_waves.png",
            )

# and also plot one of the skymaps (in this case, the vpol map the distance of the pulser)
dv.plot_skymap(reco_results["pulser_v"]["map"],
                ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_map.png",
                )
