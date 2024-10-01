import logging
logging.getLogger().setLevel(logging.INFO)

import araproc
from araproc.framework import dataset 
from araproc.framework import data_visualization as dv

# good A5 data set with pulser
# d = dataset.AraDataset(
#     path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
#     path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
# )

# a nice A2 dataset with CW (go after eventNumber 213179, which is TTree index 20695 in the burn sample file)
d = dataset.AraDataset(
    # path_to_data_file="/data/wipac/ARA/2013/filtered/full2013Data/ARA02/root/run1548/event1548.root",
    path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
    path_to_pedestal_file="/data/user/brianclark/ARA/framework_dev/software_deps/cvmfs/ara.opensciencegrid.org/trunk/ara_build/share/araCalib/ATRI/araAtriStation2Pedestals.txt"
)

for e in range(20690, 20700, 1):

    useful_event = d.get_useful_event(e)
    if useful_event.eventNumber == 213179:
        wave_bundle = d.get_waveforms(useful_event)

        dv.plot_waveform_bundle(wave_bundle, 
                            time_or_freq="freq",
                            ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_fdomain.png",
                            )
