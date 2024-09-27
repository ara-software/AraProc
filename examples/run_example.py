import logging
logging.getLogger().setLevel(logging.INFO)

import araproc
from araproc.framework import dataset 

d = dataset.AraDataset(
    station_id=5,
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)
print(d.run_number)

e = d.get_calibrated_event(0)

# for i in range(d.num_events):
#     e = d.get_calibrated_event(i)
#     del e
