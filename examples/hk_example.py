import araproc # noqa: F401
from araproc.framework import dataset
from araproc.framework import housekeeping

import logging
logging.getLogger().setLevel(logging.DEBUG)

# for later use
# d = dataset.AnalysisDataset(
#     station_id=5,
#     path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
#     do_not_calibrate = False, # turn this flag to True if you only need raw atri events
# )
# iter_start = 0
# iter_stop = 200

shk = housekeeping.SensorHkWrapper(
    path_to_hk_file = "/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/sensorHk005626.root"
)
iter_start = 0
iter_stop = shk.num_data_points


for p in range(iter_start, iter_stop, 1):

    hk = shk.get_data_point(p)

    # how to get dda/tda voltages and currents
    dda_volt = []
    dda_curr = []
    tda_volt = []
    tda_curr = []
    for dda in range(4):
        dda_volt.append(hk.getDdaVoltage(dda))
        dda_curr.append(hk.getDdaCurrent(dda))
        tda_volt.append(hk.getTdaVoltage(dda))
        tda_curr.append(hk.getTdaCurrent(dda))
