import logging
logging.getLogger().setLevel(logging.INFO)

import araproc
from araproc.framework import dataset 
from araproc.framework import data_visualization as dv

d = dataset.AraDataset(
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)

for e in range(d.num_events):

    useful_event = d.get_useful_event(e)
    if useful_event.isCalpulserEvent():
        print(e)

        wave_bundle = d.get_waveforms(useful_event)

        dv.plot_waveform_bundle(wave_bundle, 
                            time_or_freq="time",
                            ouput_file_path=f"/scratch/brianclark/station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_tdomain.png",
                            )
