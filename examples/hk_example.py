import araproc # noqa: F401
from araproc.framework import dataset
from araproc.framework import housekeeping

###################################################
###################################################
###################################################
# an example of how to read a sensor Hk file
###################################################
###################################################
###################################################

shk = housekeeping.SensorHkWrapper(
    path_to_hk_file = "/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/sensorHk005626.root"
)
for p in range(shk.num_data_points):

    hk = shk.get_data_point(p)

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

###################################################
###################################################
###################################################
# an example of how to read a event Hk file
###################################################
###################################################
###################################################

ehk = housekeeping.EventHkWrapper(
    path_to_hk_file = "/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/eventHk005626.root"
)
for p in range(ehk.num_data_points):

    hk = ehk.get_data_point(p)

    # how to get the data point unix time
    unix_time = hk.unixTime

    # the L1 rate (threshold per channel)
    # nb: only the first four channels of each dda are actually meaningful
    rates = []
    thresh = []
    for dda in range(4):
        for chan in range(4):
            rates.append(hk.getSingleChannelRateHz(dda,chan))
            thresh.append(hk.getSingleChannelThreshold(dda,chan))

###################################################
###################################################
###################################################
# an example of how to read a data config file
###################################################
###################################################
###################################################

cfg = housekeeping.ConfigFileWrapper(
    path_to_config_file = "/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/configFile.run005626.dat"
)

print(cfg.cal_pulser_info)