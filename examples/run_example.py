import numpy as np

import logging
logging.getLogger().setLevel(logging.INFO)

import araproc
from araproc.framework import dataset 
from araproc.framework import analysis_event
from araproc.framework import data_visualization as dv
from araproc.framework import waveform_utilities as wu

d = dataset.AraDataset(
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)

useful_event = d.get_useful_event(15)
ana_event = analysis_event.AraAnalysisEvent(station_id = d.station_id,
                                            run_number = d.run_number,
                                            useful_event=useful_event,
                                            )
dv.plot_analysis_event(ana_event, 
                       set_to_visualize="interpolated",
                       time_or_freq="freq",
                       ouput_file_path=f"station_{ana_event.station_id}_run_{ana_event.run_number}_event_{ana_event.event_number}.png",
                        )

# waves = ana_event.waveform_sets["interpolated"]
# wave_x, wave_y = wu.tgraph_to_arrays(waves[0])

# from araproc.analysis import dedisperse as dd

# s = dd.load_arasim_phase_response_as_spline()
# dd.dedisperse_wave(wave_x, wave_y, s)



# freqs = np.linspace(0, 1.5, 2000)
# phases = s(freqs)
# phases_wrapped = np.angle(np.exp(1j*phases))


# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1,1)

# ax.plot(data["freq"]/1E3, data["phase"])
# ax.plot(freqs, phases_wrapped, '--')

# # ax.plot(data["freq"]/1E3, np.unwrap(data["phase"]))
# # ax.plot(freqs, np.unwrap(phases_wrapped), '--')
# fig.savefig("test.png")