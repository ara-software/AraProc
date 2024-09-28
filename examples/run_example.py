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
                       ouput_file_path=f"station_{ana_event.station_id}_run_{ana_event.run_number}_event_{ana_event.event_number}.png",
                        )

waves = ana_event.waveform_sets["interpolated"]
wave_x, wave_y = wu.tgraph_to_arrays(waves[0])
print(wave_x)
