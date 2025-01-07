import araproc # noqa: F401
from araproc.framework import dataset
from araproc.framework import housekeeping

import logging
logging.getLogger().setLevel(logging.DEBUG)

# an example of how to read a sensor Hk file
shk = housekeeping.SensorHkWrapper(
    path_to_hk_file = "/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/sensorHk005626.root"
)
iter_start = 0
iter_stop = shk.num_data_points

for p in range(iter_start, iter_stop, 1):

    hk = shk.get_sensor_data_point(p)

    # how to get the data point unix time
    unix_time = hk.unixTime

    # how to get dda/tda voltages and currents
    dda_volt = []
    dda_curr = []
    dda_temp = []
    tda_volt = []
    tda_curr = []
    tda_temp = []
    for dda in range(4):
        dda_volt.append(hk.getDdaVoltage(dda))
        dda_curr.append(hk.getDdaCurrent(dda))
        dda_temp.append(hk.getDdaTemp(dda))
        tda_volt.append(hk.getTdaVoltage(dda))
        tda_curr.append(hk.getTdaCurrent(dda))
        tda_temp.append(hk.getTdaTemp(dda))

